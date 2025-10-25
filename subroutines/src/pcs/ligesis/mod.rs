use crate::{pcs::{deepfold::{self, DeepFoldCommitment, DeepFoldPCS, DeepFoldProof, DeepFoldProverCommitmentAdvice, DeepFoldProverParam, DeepFoldSRS, DeepFoldVerifierCommitmentAdvice, DeepFoldVerifierParam}, prelude::*}, rand::random_field_vector_from_rng, IOPProof, PolyIOP, SumCheck};
use arithmetic::{math::Math, VirtualPolynomial, VPAuxInfo};
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, marker::PhantomData, rand::Rng,
    sync::Arc, vec, vec::Vec, cmp::min,
};
use serde::Serialize;
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
        proof: &LigeSISProof<F>,
    ) -> F {
        F::ZERO
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
    com_mat_a_advice: DeepFoldProverCommitmentAdvice<F>,
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
    com_mat_a: DeepFoldCommitment<F>,
    com_mat_a_v_advice: DeepFoldVerifierCommitmentAdvice<F>,
    deepfold_verifier_param: DeepFoldVerifierParam<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct LigeSISProof<F: PrimeField> {
    pub com_a: DeepFoldCommitment<F>,
    pub com_bI: DeepFoldCommitment<F>,
    pub bI_check_proof: IOPProof<F>,
    pub alpha2_a_bI_r2_check_proof: IOPProof<F>,
    pub v_bI_r2_check_proof: IOPProof<F>,
    pub lookup_proof: (),
    pub com_a_proof: DeepFoldProof<F>,
    pub rs_a_proof: (),
    pub com_h_proof: DeepFoldProof<F>,
    pub com_mat_a_proof: DeepFoldProof<F>,
    pub com_bI_proof0: DeepFoldProof<F>,
    pub com_bI_proof1: DeepFoldProof<F>,
    pub com_bI_proof2: DeepFoldProof<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct LigeSISProverCommitmentAdvice<F: PrimeField> {
    pub mat_h: Vec<Vec<F>>,
    pub com_mat_h_advice: DeepFoldProverCommitmentAdvice<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct LigeSISVerifierCommitmentAdvice<F: PrimeField> {
    pub com_mat_h_v_advice: DeepFoldVerifierCommitmentAdvice<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct LigeSISCommitment<F: PrimeField> {
    pub com_mat_h: DeepFoldCommitment<F>,
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
    type VerifierCommitmentAdvice = LigeSISVerifierCommitmentAdvice<F>;

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
        let mat_a = random_field_vector_from_rng(c * 2 * (1 << (mu - log_m)), rng)
            .chunks(2 * (1 << (mu - log_m)))
            .map(|chunk| chunk.to_vec())
            .collect::<Vec<_>>();
        let deepfold_srs = DeepFoldPCS::<F>::gen_srs_for_testing(rng, log_c + mu - log_m + 1)?;
        Ok(LigeSISSRS{
            mu, log_m, rs_len, eta, c, mat_a, deepfold_srs,
        })
    }

    fn setup(
        srs: impl Borrow<Self::SRS>,
        _supported_degree: Option<usize>,
        _supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let LigeSISSRS{mu, log_m, rs_len, eta, c, mat_a, deepfold_srs} = srs.borrow().clone();
        let log_n = mu - log_m;
        let n = 1 << log_n;
        let (deepfold_prover_param, deepfold_verifier_param) = DeepFoldPCS::<F>::setup(deepfold_srs, Some(deepfold_srs.mu.clone()), Some(deepfold_srs.mu.clone()))?;
        
        let mut transcript = IOPTranscript::new(b"setup");
        let mut transcript_clone = transcript.clone();
        let (com_mat_a, com_mat_a_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&mat_a.concat()), &mut transcript)?;
        let com_mat_a_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_mat_a, &mut transcript_clone)?;

        let prover_param = LigeSISProverParam{
            mu, log_m, log_n, rs_len, eta, c, mat_a, com_mat_a_advice, deepfold_prover_param,
        };
        let verifier_param = LigeSISVerifierParam{
            mu, log_m, log_n, rs_len, eta, c, com_mat_a, com_mat_a_v_advice, deepfold_verifier_param, 
        };
        Ok((prover_param, verifier_param))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        // trim parameters
        let LigeSISProverParam{mu, log_m, log_n, rs_len, eta, c, mat_a, com_mat_a_advice, deepfold_prover_param} = prover_param.borrow().clone();
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_f = reshape(&poly.evaluations, m, n);
        // encode `F`
        let rs = ReedSolomon::new(n, rs_len);
        let mat_f_prime = mat_f.iter().map(
            |row| rs.encode(row)
        ).collect::<Vec<_>>();
        // decompose `F`
        let mat_b = transposition(&(transposition(&mat_f_prime)
            .iter()
            .map(|col| decompose_vector(col))
            .collect::<Vec<_>>()));
        // compute `H`
        let mat_h = field_mat_mul_bool_mat(&mat_a, &mat_b);

        // compute com(H)
        let (com_mat_h, com_mat_h_advice) = DeepFoldPCS::commit(deepfold_prover_param, &evals_to_arcpoly(&mat_h.concat()), transcript)?;

        Ok((LigeSISCommitment{com_mat_h}, LigeSISProverCommitmentAdvice{mat_h, com_mat_h_advice}))
    }

    fn verifier_receive_commit(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::VerifierCommitmentAdvice, PCSError> {
        let com_mat_h_v_advice = DeepFoldPCS::verifier_receive_commit(
            &verifier_param.deepfold_verifier_param, 
            &commitment.com_mat_h,
            transcript)?;
        Ok(LigeSISVerifierCommitmentAdvice { com_mat_h_v_advice })
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
        let LigeSISProverParam{mu, log_m, log_n, rs_len, eta, c, mat_a, com_mat_a_advice, deepfold_prover_param} = prover_param.borrow().clone();
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
        let LigeSISProverCommitmentAdvice{mat_h, com_mat_h_advice} = advice.clone();

        // Step 1
        let (z1, z2): (Vec<F>, Vec<F>) = (point[..log_m].to_vec(), point[log_m..].to_vec());
        let eq_z1 = get_tensor(&z1);
        
        // Step 2
        let a = (0..m).map(
            |i| (0..n).map(|j| eq_z1[i] * mat_f[i][j]).sum()
        ).collect::<Vec<F>>();
        let (com_a, com_a_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&a), transcript)?;

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", 128, 2 * n)?;
        
        // Step 4
        let mat_b_trans = transposition(&mat_b);
        let mat_bI = transposition(&I.iter().map(|&i| mat_b_trans[i].clone()).collect::<Vec<_>>());
        let bI_field = bool_vec_to_field_vec(&mat_bI.concat());
        let (com_bI, com_bI_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), transcript)?;
        
        // Step 5
        let alpha1 = transcript.get_and_append_challenge_vectors(b"alpha1", (m * mu * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;

        // Step 6
        let mut bI_check = VirtualPolynomial::new_from_mle(&evals_to_arcpoly(&bI_field), F::ONE);
        bI_check.mul_by_mle(evals_to_arcpoly(&bI_field.iter().map(|&x| x - F::ONE).collect::<Vec<F>>()), F::ONE).unwrap();
        bI_check.mul_by_mle(evals_to_arcpoly(&alpha1), F::ONE).unwrap();
        let bI_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(bI_check, transcript).unwrap();
        let r1 = bI_check_proof.point.clone();
        let com_bI_proof0 = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), &com_bI_advice, &r1, transcript)?;

        // Step 7
        let eq_alpha2_a_bI = mat_mul(&vec![get_tensor(&alpha2)], &field_mat_mul_bool_mat(&mat_a, &mat_bI));
        let v = otimes(&get_tensor(&alpha2), &(0..eta).map(|i| F::from(1u64 << i)).collect::<Vec<_>>());
        let v_bI = field_mat_mul_bool_mat(&vec![v.clone()], &mat_bI);
        let eq_alpha2_h = mat_mul(&vec![get_tensor(&alpha2)], &mat_h);
        let rs_a = rs.encode(&a);

        /*
        Lookup Argument for (I, eq_alpha2_a_bI, v_bI) in ([2n], eq_alpha2_h, rs_a)
        */
        let lookup_proof = ();
        let r2 = vec![F::ZERO; s_lambda];
        let r3 = vec![F::ZERO; 2 * n];

        // Step 8
        let alpha2_a = mat_mul(&vec![get_tensor(&alpha2)], &mat_a)[0].clone();
        let bI_r2 = field_mat_mul_bool_mat(&vec![get_tensor(&r2)], &transposition(&mat_bI))[0].clone();
        let mut alpha2_a_bI_r2_check = VirtualPolynomial::new_from_mle(&evals_to_arcpoly(&alpha2_a), F::ONE);
        alpha2_a_bI_r2_check.mul_by_mle(evals_to_arcpoly(&bI_r2), F::ONE);
        let alpha2_a_bI_r2_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(alpha2_a_bI_r2_check, transcript).unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();

        // Step 9
        let mut v_bI_r2_check = VirtualPolynomial::new_from_mle(&evals_to_arcpoly(&v), F::ONE);
        v_bI_r2_check.mul_by_mle(evals_to_arcpoly(&bI_r2), F::ONE);
        let v_bI_r2_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(v_bI_r2_check, transcript).unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();

        // Step 10
        let com_a_proof = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&a), &com_a_advice, &z2, transcript)?;
        // rs(a)(r3) = y5 to be done
        let rs_a_proof = ();
        let com_h_proof = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&mat_h.concat()), &com_mat_h_advice, &vec![r3.clone(), alpha2.clone()].concat(), transcript)?;
        let com_mat_a_proof = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&mat_a.concat()), &com_mat_a_advice, &vec![r4.clone(), alpha2.clone()].concat(), transcript)?;
        let com_bI_proof1 = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), &com_bI_advice, &vec![r2.clone(), r4.clone()].concat(), transcript)?;
        let com_bI_proof2 = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), &com_bI_advice, &vec![r2.clone(), r5.clone()].concat(), transcript)?;

        Ok(LigeSISProof {
            com_a,
            com_bI,
            bI_check_proof,
            alpha2_a_bI_r2_check_proof,
            v_bI_r2_check_proof,
            lookup_proof,
            com_a_proof,
            rs_a_proof,
            com_h_proof,
            com_mat_a_proof,
            com_bI_proof0,
            com_bI_proof1,
            com_bI_proof2,
        })
    }

    fn verify(
        verifier_param: &Self::VerifierParam,
        com: &Self::Commitment,
        point: &Self::Point,
        value: &F,
        advice: &Self::VerifierCommitmentAdvice,
        proof: &Self::Proof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        // trim parameters
        let eta = F::MODULUS_BIT_SIZE;
        let s_lambda = 128;
        let LigeSISVerifierParam{mu, log_m, log_n, rs_len, eta, c, com_mat_a, com_mat_a_v_advice, deepfold_verifier_param} = verifier_param.borrow().clone();
        let LigeSISCommitment{ com_mat_h } = com.clone();
        let LigeSISVerifierCommitmentAdvice{com_mat_h_v_advice} = advice.clone();
        let (m, n) = (1 << log_m, 1 << log_n);
        let LigeSISProof {
            com_a,
            com_bI,
            bI_check_proof,
            alpha2_a_bI_r2_check_proof,
            v_bI_r2_check_proof,
            lookup_proof,
            com_a_proof,
            rs_a_proof,
            com_h_proof,
            com_mat_a_proof,
            com_bI_proof0,
            com_bI_proof1,
            com_bI_proof2,
        } = proof.clone();

        // Step 1
        let (z1, z2): (Vec<F>, Vec<F>) = (point[..log_m].to_vec(), point[log_m..].to_vec());

        // Step 2
        let com_a_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_a, transcript)?;

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", 128, 2 * n)?;

        // Step 4
        let com_bI_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_bI, transcript)?;

        // Step 5
        let alpha1 = transcript.get_and_append_challenge_vectors(b"alpha1", (m * mu * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;

        // Step 6
        let bI_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&bI_check_proof);
        let _sum_check_claim0 = <PolyIOP<F> as SumCheck<F>>::verify(bI_check_sum, &bI_check_proof, &VPAuxInfo{
            max_degree: 1, 
            num_variables: m * eta * s_lambda, 
            phantom: PhantomData::<F>::default()
        }, transcript).unwrap();
        let r1 = bI_check_proof.point.clone();
        if !DeepFoldPCS::verify(
                &deepfold_verifier_param, &com_bI, &r1, 
                &DeepFoldPCS::compute_value_from_proof(&r1, &com_bI_proof0), 
                &com_bI_v_advice, &com_bI_proof0, transcript).unwrap() {
            return Ok(false);
        }

        // Step 7
        let r2 = vec![F::ZERO; s_lambda];
        let r3 = vec![F::ZERO; 2 * n];

        // Step 8
        let alpha2_a_bI_r2_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&alpha2_a_bI_r2_check_proof);
        let _sum_check_claim1 = <PolyIOP<F> as SumCheck<F>>::verify(alpha2_a_bI_r2_check_sum, &alpha2_a_bI_r2_check_proof, &VPAuxInfo{
            max_degree: 1, 
            num_variables: m * eta, 
            phantom: PhantomData::<F>::default()
        }, transcript).unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();
        
        // Step 9
        let v_bI_r2_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&v_bI_r2_check_proof);
        let _sum_check_claim2 = <PolyIOP<F> as SumCheck<F>>::verify(v_bI_r2_check_sum, &v_bI_r2_check_proof, &VPAuxInfo{
            max_degree: 1, 
            num_variables: m * eta, 
            phantom: PhantomData::<F>::default()
        }, transcript).unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();

        // Step 10
        let coms = vec![
            com_a, com_mat_h, com_mat_a, com_bI.clone(), com_bI,
        ];
        let points = vec![
            z2, vec![r3.clone(), alpha2.clone()].concat(), 
            vec![r4.clone(), alpha2.clone()].concat(), 
            vec![r2.clone(), r4.clone()].concat(),
            vec![r2.clone(), r5.clone()].concat(),
        ];
        let proofs = vec![
            com_a_proof, com_h_proof, com_mat_a_proof, com_bI_proof1, com_bI_proof2,
        ];
        let values = points.iter().zip(proofs.iter())
            .map(|(point, proof)| DeepFoldPCS::compute_value_from_proof(point, proof))
            .collect::<Vec<_>>();
        let advices = vec![
            com_a_v_advice, com_mat_h_v_advice, com_mat_a_v_advice, com_bI_v_advice.clone(), com_bI_v_advice,
        ];
        
        for ((((com, point), value), proof), advice) 
            in coms.iter().zip(points.iter()).zip(values.iter()).zip(proofs.iter()).zip(advices.iter()) {
                if DeepFoldPCS::verify(&deepfold_verifier_param, com, point, value, advice, proof, transcript)? {
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
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| if b[k][j] { a[i][k] } else { F::ZERO }).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

fn bool_mat_mul_field_mat<F: PrimeField>( a: &Vec<Vec<bool>>, b: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert!(m == b.len());
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| if a[i][k] { b[k][j] } else { F::ZERO }).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

fn evals_to_arcpoly<F: PrimeField>( a: &Vec<F> ) -> Arc<DenseMultilinearExtension<F>> {
    Arc::new(DenseMultilinearExtension::<F>::from_evaluations_vec(a.len().ilog2() as usize, a.clone()))
}

fn otimes<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> Vec<F> {
    a.iter().map(
        |x| b.iter().map(
            |y| (*x) * (*y)
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>().concat()
}

fn bool_vec_to_field_vec<F: PrimeField>( a: &Vec<bool> ) -> Vec<F> {
    (0..a.len()).map(
        |i| if a[i] {F::ONE} else {F::ZERO}
    ).collect::<Vec<_>>()
}