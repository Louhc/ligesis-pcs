use crate::pcs::{deepfold::{self, DeepFoldCommitment, DeepFoldPCS, DeepFoldProverCommitmentAdvice, DeepFoldProverParam, DeepFoldSRS, DeepFoldVerifierParam}, prelude::*};
use arithmetic::math::Math;
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, marker::PhantomData, rand::Rng,
    sync::Arc, vec, vec::Vec, cmp::min,
};
use transcript::IOPTranscript;

#[cfg(test)]
mod tests;

mod rand;
use rand::*;

/// LigeSIS Polynomial Commitment Scheme
pub struct LigeSISPCS<F: PrimeField> {
    #[doc(hidden)]
    phantom: PhantomData<F>,
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

#[derive(Clone, Debug)]
pub struct LigeSISSRS<F: PrimeField> {
    mu: usize,
    log_m: usize,
    rs_len: usize,
    eta: usize,
    c: usize,
    mat_a: Vec<Vec<F>>,
    deepfold_srs: DeepFoldSRS<F>,
}

#[derive(Clone)]
pub struct LigeSISProverParam<F: PrimeField> {
    mu: usize,
    log_m: usize,
    log_n: usize,
    rs_len: usize,
    eta: usize,
    c: usize,
    mat_a: Vec<Vec<F>>,
    deepfold_prover_param: DeepFoldProverParam<F>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct LigeSISVerifierParam<F: PrimeField> {
    mu: usize,
    log_m: usize,
    log_n: usize,
    rs_len: usize,
    eta: usize,
    c: usize,
    com_a: F,
    deepfold_verifier_param: DeepFoldVerifierParam<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct LigeSISProof<F: PrimeField> {
    pub f0: Vec<F>, // r^T * A
    pub f1: Vec<F>, // u0^T * A
    // pub mt_proofs: Vec<Vec<Byte32>>,
    pub cols: Vec<Vec<F>>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
struct LigeSISProverCommitmentAdvice<F: PrimeField> {
    com_h_advice: DeepFoldProverCommitmentAdvice<F>
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
struct LigeSISCommitment<F: PrimeField> {
    com_h: DeepFoldCommitment<F>,
}

impl<F: PrimeField> PolynomialCommitmentScheme<F> for LigeSISPCS<F> {
    // Parameters
    type ProverParam = LigeSISProverParam<F>;
    type VerifierParam = LigeSISVerifierParam<F>;
    type SRS = LigeSISSRS<F>;
    // Polynomial and its associated types
    type Polynomial = Arc<DenseMultilinearExtension<F>>;
    type ProverCommitmentAdvice = LigeSISProverCommitmentAdvice<F>;
    type Point = Vec<F>;
    type Evaluation = F;
    // Commitments and proofs
    type Commitment = LigeSISCommitment<F>;
    type Proof = LigeSISProof<F>;
    type BatchProof = ();

    fn gen_srs_for_testing<R: Rng>(
        rng: &mut R, 
        log_size: usize
    ) -> Result<Self::SRS, PCSError> {
        // MultilinearUniversalParams::<E>::gen_srs_for_testing(rng, log_size)
        let mu = log_size;
        let log_m = mu / 2;
        let rs_len = (1 << log_m) * 2;
        let eta = 64;
        let log_c = 4;
        let c = 1 << log_c;
        let mat_a = (0..c).map(
            |_| random_field_vector_from_rng((1 << log_m) * eta, rng)
        ).collect::<Vec<_>>();
        let deepfold_srs = DeepFoldPCS::<F>::gen_srs_for_testing(rng, log_c + mu - log_m + 1)?;
        Ok(LigeSISSRS{
            mu, log_m, rs_len, eta, c, mat_a, deepfold_srs,
        })
    }

    fn trim(
        srs: impl Borrow<Self::SRS>,
        _supported_degree: Option<usize>,
        _supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let LigeSISSRS{mu, log_m, rs_len, eta, c, mat_a, deepfold_srs} = srs.borrow().clone();
        let log_n = mu - log_m;
        let com_a = F::ZERO;
        let (deepfold_prover_param, deepfold_verifier_param) = DeepFoldPCS::<F>::trim(deepfold_srs, Some(deepfold_srs.mu), Some(deepfold_srs.mu))?;
        let prover_param = LigeSISProverParam{
            mu, log_m, log_n, rs_len, eta, c, mat_a, deepfold_prover_param,
        };
        let verifier_param = LigeSISVerifierParam{
            mu, log_m, log_n, rs_len, eta, c, com_a, deepfold_verifier_param, 
        };
        Ok((prover_param, verifier_param))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        // trim parameters
        let LigeSISProverParam{mu, log_m, log_n, rs_len, eta, c, mat_a, deepfold_prover_param} = prover_param.borrow().clone();
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
        let (com_h, com_h_advice) = DeepFoldPCS::<F>::commit(deepfold_prover_param, &evals_to_arcpoly(&mat_h.concat()), transcript)?;

        Ok((LigeSISCommitment{com_h}, LigeSISProverCommitmentAdvice{com_h_advice}))
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        let eta = F::MODULUS_BIT_SIZE;
        let s_lambda = 128;
        let &LigeSISProverParam{mu, log_m, log_n, rs_len, eta, c, mat_a, deepfold_prover_param} = prover_param.borrow();
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_f = reshape(&poly.evaluations, m, n);
        let rs = ReedSolomon::new(n, rs_len);
        let mat_f_prime = mat_a.iter().map(
            |row| rs.encode(row)
        ).collect::<Vec<_>>();
        let mat_b = transposition(&(transposition(&mat_f_prime)
            .iter()
            .map(|col| decompose_vector(col))
            .collect::<Vec<_>>()));

        let (z1, z2) = (point[..log_m].to_vec(), point[log_m..].to_vec());
        let eq_z1 = get_tensor(&z1);
        let a = (0..m).map(
            |i| (0..n).map(|j| eq_z1[i] * mat_f[i][j]).sum()
        ).collect::<Vec<F>>();

        let I = transcript.get_and_append_challenge_indices(b"I", 128, 2 * n)?;
        let mat_b_trans = transposition(&mat_b);
        let mat_b_I = transposition(&I.iter().map(|&i| mat_b_trans[i]).collect::<Vec<_>>());

        let (com_b_I, com_b_I_advice) = DeepFoldPCS::commit(deepfold_prover_param, &evals_to_arcpoly(&bool_vec_to_field_vec(&mat_b_I.concat())), transcript)?;
        
        let alpha1 = transcript.get_and_append_challenge_vectors(b"alpha1", (m * mu * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;


        Ok(LigeSISProof { f0: (), f1: (), cols: () })
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

fn evals_to_arcpoly<F: PrimeField>( a: &Vec<F> ) -> Arc<DenseMultilinearExtension<F>> {
    Arc::new(DenseMultilinearExtension::<F>::from_evaluations_vec(a.len().ilog2() as usize, a.clone()))
}

fn otimes<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> Vec<F> {
    a.iter().map(
        |x| b.iter().map(
            |y| x * y
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

fn bool_vec_to_field_vec<F: PrimeField>( a: &Vec<bool> ) -> Vec<F> {
    (0..a.len()).map(
        |i| if a[i] {F::ONE} else {F::ZERO}
    ).collect::<Vec<_>>()
}