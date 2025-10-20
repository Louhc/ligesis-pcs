use crate::pcs::prelude::*;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{DenseMultilinearExtension};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, marker::PhantomData, rand::Rng,
    sync::Arc, vec, vec::Vec, cmp::min,
};
use transcript::IOPTranscript;

mod rscode;
use rscode::ReedSolomon;
mod rand;
use rand::*;

/// LigeSIS Polynomial Commitment Scheme
pub struct LigeSISPCS<F: PrimeField> {
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct LigeSISProof<F: PrimeField> {
    pub f0: Vec<F>, // r^T * A
    pub f1: Vec<F>, // u0^T * A
    // pub mt_proofs: Vec<Vec<Byte32>>,
    pub cols: Vec<Vec<F>>,
}

impl<F: PrimeField> LigeSISPCS<F> {
    pub fn compute_value_from_proof(
        log_m0: usize,
        point: &Vec<F>,
        proof: &LigeroProof<F>,
    ) -> F {
        let u1 = get_tensor(&point[log_m0..].to_vec());
        (0..u1.len()).map(|i| proof.f1[i] * u1[i]).sum::<F>()
    }
}

fn transposition<F: Copy>( mat: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    (0..mat.len()).map(
        |i| mat.iter().map(|row| row[i]).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

fn decompose<F: PrimeField>( x: &F ) -> Vec<bool> {
    x.into_bigint().to_bits_be()
}

fn decompose_vector<F: PrimeField>( v: &Vec<F> ) -> Vec<bool> {
    let mut res = Vec::new();
    for x in v.iter() {
        res.append(&mut decompose(x));
    }
    res
}

fn mat_mul<F: PrimeField>( a: &Vec<Vec<F>>, b: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert!(m == b.len());
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| a[i][k] * b[k][j]).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

fn field_mat_mul_bool_mat<F: PrimeField>( a: &Vec<Vec<F>>, b: &Vec<Vec<bool>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert!(m == b.len());
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| if b[k][j] { a[i][k] } else { F::ZERO }).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

#[derive(Clone, Debug)]
struct LigeSISSRS<F: PrimeField> {
    mu: usize,
    log_m: usize,
    rs_len: usize,
    eta: usize,
    c: usize,
    mat_a: Vec<Vec<F>>,
}

#[derive(Clone)]
struct LigeSISProverParam<F: PrimeField> {
    mu: usize,
    log_m: usize,
    log_n: usize,
    rs_len: usize,
    eta: usize,
    c: usize,
    mat_a: Vec<Vec<F>>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
struct LigeSISVerifierParam<F: PrimeField> {
    mu: usize,
    log_m: usize,
    log_n: usize,
    rs_len: usize,
    eta: usize,
    c: usize,
    com_a: F,
}

impl<F: PrimeField> PolynomialCommitmentScheme<F> for LigeSISPCS<F> {
    // Parameters
    type ProverParam = LigeSISProverParam<F>;
    type VerifierParam = LigeSISVerifierParam<F>;
    type SRS = LigeSISSRS<F>; // (num of variables, length of RS code)
    // Polynomial and its associated types
    type Polynomial = Arc<DenseMultilinearExtension<F>>;
    type ProverCommitmentAdvice = (); // merkle tree structure
    type Point = Vec<F>;
    type Evaluation = F;
    // Commitments and proofs
    type Commitment = (); // merkle tree root
    type Proof = LigeroProof<F>; // merkle tree paths, columes of `E`
    type BatchProof = (); // 

    fn gen_srs_for_testing<R: Rng>(
        rng: &mut R, 
        log_size: usize
    ) -> Result<Self::SRS, PCSError> {
        // MultilinearUniversalParams::<E>::gen_srs_for_testing(rng, log_size)
        let mu = log_size;
        let log_m = mu / 2;
        let rs_len = (1 << log_m) * 2;
        let eta = 64;
        let c = 2;
        let mat_a = (0..c).map(
            |_| random_field_vector_from_rng((1 << log_m) * eta, rng)
        ).collect::<Vec<_>>();
        Ok(LigeSISSRS{
            mu, log_m, rs_len, eta, c, mat_a
        })
    }

    fn trim(
        srs: impl Borrow<Self::SRS>,
        _supported_degree: Option<usize>,
        _supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let &LigeSISSRS{mu, log_m, rs_len, eta, c, mat_a} = srs.borrow();
        let log_n = mu - log_m;
        let com_a = F::ZERO;
        let prover_param = LigeSISProverParam{
            mu, log_m, log_n, rs_len, eta, c, mat_a
        };
        let verifier_param = LigeSISVerifierParam{
            mu, log_m, log_n, rs_len, eta, c, com_a
        };
        Ok((prover_param, verifier_param))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        // trim parameters
        let &LigeSISProverParam{mu, log_m, log_n, rs_len, eta, c, mat_a} = prover_param.borrow();
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_f = reshape(&poly.evaluations, m, n);

        // encode `F`
        let rs = ReedSolomon::new(n, rs_len);
        let mat_f_prime = mat_f.iter().map(
            |row| rs.encode(row)
        ).collect::<Vec<_>>();

        let mat_b = transposition(&(transposition(&mat_f_prime)
            .iter()
            .map(|col| decompose_vector(col))
            .collect::<Vec<_>>()));

        let mat_h = field_mat_mul_bool_mat(&mat_a, &mat_b);
        let com_h = ();

        Ok((com_h, ()))
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        // trim parameters
        let &LigeSISProverParam{mu, log_m, log_n, rs_len, eta, c, mat_a} = prover_param.borrow();
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_a = reshape(&poly.evaluations, m, n);
        
        // encode `A`
        let rs = ReedSolomon::new(n, rs_len);
        let mat_e = mat_a.iter().map(
            |row| rs.encode(row)
        ).collect::<Vec<_>>();

        // generate `r` and compute the tensor vector `u0`
        let r = transcript.get_and_append_challenge_vectors(b"r", m)?;
        let u0 = get_tensor(&point[..log_m].to_vec());

        // compute `rA` and `u0A` and compute msg
        let f0: Vec<F> = (0..n).map(
            |j| (0..m).map(|i| r[i] * mat_a[i][j]).sum()
        ).collect::<Vec<_>>();
        let f1: Vec<F> = (0..n).map(
            |j| (0..m).map(|i| u0[i] * mat_a[i][j]).sum()
        ).collect::<Vec<_>>();
        // let msg: Vec<F> = { let mut f = f0.clone(); f.append(&mut f1.clone()); f };
        
        // get merkle tree on columes
        let mt = advice;
        
        // generate lambda indices and alpha
        let idx: Vec<usize> = transcript.get_and_append_challenge_indices(b"idx", min(128, m1 << 1), m1 << 1)?;        
        let _alpha: F = transcript.get_and_append_challenge(b"alpha")?;

        // trim all needed columes and compute merkle paths
        let cols = idx.iter().map(
            |&i| mat_e.iter().map(|row| row[i]).collect::<Vec<_>>()
        ).collect::<Vec<_>>();
        let mt_proofs = idx.iter().map(
            |&i| mt.prove(i)
        ).collect::<Vec<_>>();
        
        Ok(LigeroProof{f0, f1, cols, mt_proofs})
    }

    fn verify(
        verifier_param: &Self::VerifierParam,
        com: &Self::Commitment,
        point: &Self::Point,
        value: &F,
        proof: &Self::Proof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        // trim parameters
        let &(log_n, log_m0, rs_len) = verifier_param;
        let log_m1 = log_n - log_m0;
        let (n, m0, m1) = (1 << log_n, 1 << log_m0, 1 << log_m1);
        let f0 = proof.f0.clone();
        let f1 = proof.f1.clone();

        // generate the challenge and compuate the tensor vector
        let r = transcript.get_and_append_challenge_vectors(b"r", m0)?;
        let (u0, u1) = (get_tensor(&point[..log_m0].to_vec()), get_tensor(&point[log_m0..].to_vec()));

        // check if the final value is correctly computed
        if (0..m1).map(|i| f1[i] * u1[i]).sum::<F>() != *value {
            return Ok(false);
        }

        // choose lambda columes
        // generate a random value alpha and batch `f0 + alpha * f1`
        let idx: Vec<usize> = transcript.get_and_append_challenge_indices(b"idx", min(128, m1 << 1), m1 << 1)?;
        let alpha: F = transcript.get_and_append_challenge(b"alpha")?;
        let f = (0..m1).map(|i| f0[i] + alpha * f1[i]).collect();
        
        // encode `f`
        let rs = ReedSolomon::<F>::new(m1, rs_len);
        let enc = rs.encode(&f);
        let enc_i = idx.iter().map(
            |&i| enc[i]
        ).collect::<Vec<_>>();

        // check if `Enc(f)` and `(r + alpha * u0)^T E` meet at lambda points
        let cmp_i = (0..idx.len()).map(
            |i| (0..m0).map(|j| proof.cols[i][j] * (r[j] + alpha * u0[j])).sum::<F>()
        ).collect::<Vec<_>>();
        if cmp_i != enc_i {
            return Ok(false);
        }

        // check merkle paths
        for i in 0..idx.len() {
            if !MerkleTree::verify(com, idx[i], 
                &compute_sha256_row(&proof.cols[i]), 
                &proof.mt_proofs[i]) {
                return Ok(false);
            }
        }

        return Ok(true);
    }
}

pub fn reshape<F: PrimeField>( a: &Vec<F>, n: usize, m: usize ) -> Vec<Vec<F>> {
    (0..n).map(
        |i| (0..m).map(
            |j| if i * m + j < a.len() { a[i * m + j] } else { F::ZERO }
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn get_tensor<F: PrimeField>( r: &Vec<F> ) -> Vec<F> {
    let mut res = vec![F::ONE];
    for i in 0..r.len() {
        let mut tmp = res.iter().map(|&x| x * r[i]).collect::<Vec<_>>();
        res = res.iter().map(|&x| x * (F::ONE - r[i])).collect::<Vec<_>>();
        res.append(&mut tmp);
    }
    res
}

#[cfg(test)]
mod tests;