use crate::{pcs::prelude::*, IOPProof, PolyIOP, SumCheck};
use arithmetic::{math::Math, VirtualPolynomial, VPAuxInfo};
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, marker::PhantomData, rand::Rng,
    sync::Arc, vec, vec::Vec, cmp::min, cmp::max,
};
use serde::Serialize;
use transcript::IOPTranscript;

mod types;
use types::*;

#[cfg(test)]
mod tests;

/// LigeSIS Polynomial Commitment Scheme
pub struct LigeSISPCS<F: PrimeField> {
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField> LigeSISPCS<F> {
    pub fn compute_value_from_proof(
        log_n: usize,
        point: &Vec<F>,
        proof: &LigeSISProof<F>,
    ) -> F {
        DeepFoldPCS::compute_value_from_proof(
            &point[..log_n].to_vec(),
            &proof.com_a_proof,
        )
    }
}

#[derive(Clone, Debug)]
pub struct LigeSISSRS<F: PrimeField> {
    lambda: usize,
    eta: usize,
    mu: usize,
    log_m: usize,
    rs_len: usize,
    c: usize,
    mat_a: Vec<Vec<F>>,
    deepfold_srs: DeepFoldSRS<F>,
}

#[derive(Clone)]
pub struct LigeSISProverParam<F: PrimeField> {
    eta: usize,
    s_lambda: usize,
    mu: usize,
    log_m: usize,
    log_n: usize,
    rs_len: usize,
    c: usize,
    mat_a: Vec<Vec<F>>,
    com_mat_a_advice: DeepFoldProverCommitmentAdvice<F>,
    deepfold_prover_param: DeepFoldProverParam<F>,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct LigeSISVerifierParam<F: PrimeField> {
    eta: usize,
    s_lambda: usize,
    mu: usize,
    log_m: usize,
    log_n: usize,
    rs_len: usize,
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
    pub com_rs_a: DeepFoldCommitment<F>,
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
        let eta = F::ONE.into_bigint().to_bits_be().len();
        let lambda = 128usize;
        let mu = log_size;
        let log_m = if log_size < 8 {0} else { (log_size - 8) / 2};
        let rs_len = (1 << (mu - log_m)) * 2;
        let log_c = 3;
        let c = 1 << log_c;
        let log_eta = eta.ilog2() as usize;
        let log_n = mu - log_m;
        let log_s_lambda = lambda.ilog2() as usize;

        let mat_a = random_field_vector_from_rng(c * (1 << log_m) * eta, rng)
            .chunks(eta * (1 << log_m))
            .map(|chunk| chunk.to_vec())
            .collect::<Vec<_>>();

        let deepfold_srs = DeepFoldPCS::<F>::gen_srs_for_testing(rng, max(max(log_c, log_s_lambda) + log_m + log_eta, log_c + 1 + log_n))?;
        Ok(LigeSISSRS{
            eta, lambda, mu, log_m, rs_len, c, mat_a, deepfold_srs,
        })
    }

    fn setup(
        srs: impl Borrow<Self::SRS>,
        _supported_degree: Option<usize>,
        _supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let LigeSISSRS{eta, lambda, mu, log_m, rs_len, c, mat_a, deepfold_srs} = srs.borrow().clone();
        let log_n = mu - log_m;
        let n = 1 << log_n;
        let s_lambda = min(lambda, rs_len);
        let (deepfold_prover_param, deepfold_verifier_param) = DeepFoldPCS::<F>::setup(deepfold_srs, Some(deepfold_srs.max_mu.clone()), Some(deepfold_srs.max_mu.clone()))?;
        
        let mut transcript = IOPTranscript::new(b"setup");
        let mut transcript_clone = transcript.clone();
        let (com_mat_a, com_mat_a_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&mat_a.concat()), &mut transcript)?;
        
        let com_mat_a_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_mat_a, &mut transcript_clone)?;
        let prover_param = LigeSISProverParam{
            eta, s_lambda, mu, log_m, log_n, rs_len, c, mat_a, com_mat_a_advice, deepfold_prover_param,
        };
        let verifier_param = LigeSISVerifierParam{
            eta, s_lambda, mu, log_m, log_n, rs_len, c, com_mat_a, com_mat_a_v_advice, deepfold_verifier_param, 
        };
        Ok((prover_param, verifier_param))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        // trim parameters
        let LigeSISProverParam{eta, s_lambda, mu, log_m, log_n, rs_len, c, mat_a, com_mat_a_advice, deepfold_prover_param} = prover_param.borrow().clone();
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_f = reshape(&poly.evaluations, m, n);
        // encode `F`
        let start = std::time::Instant::now();
        let rs = ReedSolomon::new(n, rs_len);
        let mat_f_prime = mat_f.iter().map(
            |row| rs.encode(row)
        ).collect::<Vec<_>>();
        println!("RS(F): {} s", start.elapsed().as_secs_f64());

        // decompose `F`
        let start = std::time::Instant::now();
        let mat_b = decompose_mat_by_col(&mat_f_prime);
        println!("Decompose(F): {} s", start.elapsed().as_secs_f64());

        // compute `H`
        let start = std::time::Instant::now();
        let mat_h = field_mat_mul_bool_mat(&mat_a, &mat_b);
        println!("SIS_Hash(B): {} s", start.elapsed().as_secs_f64());

        // compute com(H)         
        let start = std::time::Instant::now();
        let (com_mat_h, com_mat_h_advice) = DeepFoldPCS::commit(deepfold_prover_param, &evals_to_arcpoly(&mat_h.concat()), transcript)?;
        println!("DeepFold.Commit(H): {} s", start.elapsed().as_secs_f64());
        
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
        Ok(LigeSISVerifierCommitmentAdvice {com_mat_h_v_advice})
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        // let start = std::time::Instant::now();
        let LigeSISProverParam{eta, s_lambda, mu, log_m, log_n, rs_len, c, mat_a, com_mat_a_advice, deepfold_prover_param} = prover_param.borrow().clone();
        let (m, n) = (1 << log_m, 1 << log_n);

        assert_eq!(mu, log_m + log_n);
        assert_eq!(poly.num_vars, mu);

        let mat_f = reshape(&poly.evaluations, m, n);
        let rs = ReedSolomon::new(n, rs_len);
        let mat_f_prime = mat_f.iter().map(
            |row| rs.encode(row)
        ).collect::<Vec<_>>();
        let mat_b = transposition(&(transposition(&mat_f_prime)
            .iter()
            .map(|col| decompose_vector(col))
            .collect::<Vec<_>>()));
        let LigeSISProverCommitmentAdvice{mat_h, com_mat_h_advice} = advice.clone();

        // Step 1
        let (z1, z2): (Vec<F>, Vec<F>) = (point[log_n..].to_vec(), point[..log_n].to_vec());
        let eq_z1 = get_tensor(&z1);
        
        // Step 2
        let a = (0..n).map(
            |j| (0..m).map(|i| eq_z1[i] * mat_f[i][j]).sum()
        ).collect::<Vec<F>>();
        let (com_a, com_a_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&a), transcript)?;

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;
        
        // Step 4
        let mat_b_trans = transposition(&mat_b);
        let mat_bI = transposition(&I.iter().map(|&i| mat_b_trans[i].clone()).collect::<Vec<_>>());
        let bI_field = bool_vec_to_field_vec(&mat_bI.concat());
        let (com_bI, com_bI_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), transcript)?;
        
        // Step 5
        let alpha1 = transcript.get_and_append_challenge_vectors(b"alpha1", (m * eta * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
        let alpha3 = transcript.get_and_append_challenge_vectors(b"alpha3", rs_len.ilog2() as usize)?;

        // Step 6
        let mut bI_check = VirtualPolynomial::new(bI_field.len().ilog2() as usize);
        bI_check.add_mle_list([
            evals_to_arcpoly(&bI_field),
            evals_to_arcpoly(&bI_field.iter().map(|&x| x - F::ONE).collect::<Vec<F>>()),
            evals_to_arcpoly(&get_tensor(&alpha1)),
        ], F::ONE).unwrap();
        let bI_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(bI_check, transcript).unwrap();
        let r1 = bI_check_proof.point.clone();
        let com_bI_proof0 = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), &com_bI_advice, &r1, transcript)?;
        
        // Step 7 Check rs_a
        let rs_a = rs.encode(&a);
        let (com_rs_a, com_rs_a_advice) = DeepFoldPCS::commit(&deepfold_prover_param, &evals_to_arcpoly(&rs_a), transcript)?;
        
        // let 

        // let mut rs_a_check = VirtualPolynomial::new(bI_field.len().ilog2() as usize);
        // rs_a_check.add_mle_list([
        //     evals_to_arcpoly(&bI_field),
        //     evals_to_arcpoly(&bI_field.iter().map(|&x| x - F::ONE).collect::<Vec<F>>()),
        //     evals_to_arcpoly(&get_tensor(&alpha1)),
        // ], F::ONE).unwrap();

        // Step 8 Lookup Argument
        let eq_alpha2_a_bI = mat_mul(&vec![get_tensor(&alpha2)], &field_mat_mul_bool_mat(&mat_a, &mat_bI));
        let v = otimes(&get_tensor(&z1), &(0..eta).map(|i| F::from(2u64).pow([i as u64])).collect::<Vec<_>>());
        let v_bI = field_mat_mul_bool_mat(&vec![v.clone()], &mat_bI);
        let eq_alpha2_h = mat_mul(&vec![get_tensor(&alpha2)], &mat_h);

        /*
        Lookup Argument for (I, eq_alpha2_a_bI, v_bI) in ([2n], eq_alpha2_h, rs_a)
        */
        let lookup_proof = ();
        let r2 = vec![F::ZERO; s_lambda.ilog2() as usize];
        let r3 = vec![F::ZERO; 1 + log_n];

        // Step 9
        let alpha2_a = mat_mul(&vec![get_tensor(&alpha2)], &mat_a)[0].clone();
        let bI_r2 = field_mat_mul_bool_mat(&vec![get_tensor(&r2)], &transposition(&mat_bI))[0].clone();
        let mut alpha2_a_bI_r2_check = VirtualPolynomial::new(mat_bI.len().ilog2() as usize);
        alpha2_a_bI_r2_check.add_mle_list([
            evals_to_arcpoly(&alpha2_a),
            evals_to_arcpoly(&bI_r2),
        ], F::ONE).unwrap();
        let alpha2_a_bI_r2_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(alpha2_a_bI_r2_check, transcript).unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();
        
        // Step 10
        let mut v_bI_r2_check = VirtualPolynomial::new(v.len().ilog2() as usize);
        v_bI_r2_check.add_mle_list([
            evals_to_arcpoly(&v),
            evals_to_arcpoly(&bI_r2),
        ], F::ONE).unwrap();
        let v_bI_r2_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(v_bI_r2_check, transcript).unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();

        // Step 11
        let com_a_proof = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&a), &com_a_advice, &z2, transcript)?;
        // rs(a)(r3) = y5 to be done
        let rs_a_proof = ();
        let com_h_proof = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&mat_h.concat()), &com_mat_h_advice, &vec![r3.clone(), alpha2.clone()].concat(), transcript)?;
        let com_mat_a_proof = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&mat_a.concat()), &com_mat_a_advice, &vec![r4.clone(), alpha2.clone()].concat(), transcript)?;
        let com_bI_proof1 = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), &com_bI_advice, &vec![r2.clone(), r4.clone()].concat(), transcript)?;
        let com_bI_proof2 = DeepFoldPCS::open(&deepfold_prover_param, &evals_to_arcpoly(&bI_field), &com_bI_advice, &vec![r2.clone(), r5.clone()].concat(), transcript)?;
        // println!("{} ms", start.elapsed().as_millis());
        // assert!(false);
        

        Ok(LigeSISProof {
            com_a,
            com_bI,
            com_rs_a,
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
        let LigeSISVerifierParam{eta, s_lambda, mu, log_m, log_n, rs_len, c, com_mat_a, com_mat_a_v_advice, deepfold_verifier_param} = verifier_param.borrow().clone();
        let LigeSISCommitment{ com_mat_h } = com.clone();
        let LigeSISVerifierCommitmentAdvice{com_mat_h_v_advice} = advice.clone();
        let (m, n) = (1 << log_m, 1 << log_n);
        let LigeSISProof {
            com_a,
            com_bI,
            com_rs_a,
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
        let (z1, z2): (Vec<F>, Vec<F>) = (point[log_n..].to_vec(), point[..log_n].to_vec());

        // Step 2
        let com_a_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_a, transcript)?;

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;

        // Step 4
        let com_bI_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_bI, transcript)?;

        // Step 5
        let alpha1 = transcript.get_and_append_challenge_vectors(b"alpha1", (m * eta * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
        let alpha3 = transcript.get_and_append_challenge_vectors(b"alpha3", rs_len.ilog2() as usize)?;
        // println!("verify : {}", alpha2[0]);

        // Step 6
        let bI_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&bI_check_proof);
        let _sum_check_claim0 = <PolyIOP<F> as SumCheck<F>>::verify(bI_check_sum, &bI_check_proof, &VPAuxInfo{
            max_degree: 3, 
            num_variables: (m * eta * s_lambda).ilog2() as usize, 
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
        let com_rs_a_v_advice = DeepFoldPCS::verifier_receive_commit(&deepfold_verifier_param, &com_rs_a, transcript)?;

        // Step 8
        let r2 = vec![F::ZERO; s_lambda.ilog2() as usize];
        let r3 = vec![F::ZERO; 1 + log_n];

        // Step 9
        let alpha2_a_bI_r2_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&alpha2_a_bI_r2_check_proof);
        let _sum_check_claim1 = <PolyIOP<F> as SumCheck<F>>::verify(alpha2_a_bI_r2_check_sum, &alpha2_a_bI_r2_check_proof, &VPAuxInfo{
            max_degree: 2, 
            num_variables: (m * eta).ilog2() as usize, 
            phantom: PhantomData::<F>::default()
        }, transcript).unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();
        
        // Step 10
        let v_bI_r2_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&v_bI_r2_check_proof);
        let _sum_check_claim2 = <PolyIOP<F> as SumCheck<F>>::verify(v_bI_r2_check_sum, &v_bI_r2_check_proof, &VPAuxInfo{
            max_degree: 2, 
            num_variables: (m * eta).ilog2() as usize, 
            phantom: PhantomData::<F>::default()
        }, transcript).unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();

        // Step 11
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
        if DeepFoldPCS::compute_value_from_proof(&points[0], &proofs[0]) != *value {
            return Ok(false);
        }
        
        for ((((com, point), value), proof), advice) 
            in coms.iter().zip(points.iter()).zip(values.iter()).zip(proofs.iter()).zip(advices.iter()) {
                if !DeepFoldPCS::verify(&deepfold_verifier_param, com, point, value, advice, proof, transcript)? {
                    return Ok(false);
                }
        }
        
        return Ok(true);
    }
}
