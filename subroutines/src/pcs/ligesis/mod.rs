use crate::{pcs::prelude::*, IOPProof, PolyIOP, SumCheck};
use arithmetic::{math::Math, VPAuxInfo, VirtualPolynomial};
use ark_ff::{BigInteger, PrimeField};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow,
    cmp::{max, min},
    marker::PhantomData,
    rand::Rng,
    sync::Arc,
    vec,
    vec::Vec,
};
use serde::Serialize;
use transcript::IOPTranscript;

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

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
    pub fn compute_value_from_proof(log_n: usize, point: &Vec<F>, proof: &LigeSISProof<F>) -> F {
        proof.deepfold_batched_proof.evals[0]
        // DeepFoldPCS::compute_value_from_proof(
        //     &point[..log_n].to_vec(),
        //     &proof.com_a_proof,
        // )
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
    rs: ReedSolomon<F>,
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
    g: F,
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
    pub rs_a_check_proof: IOPProof<F>,
    pub mat_g_check_proofs: Vec<IOPProof<F>>,
    pub lookup_proof: (),
    pub deepfold_batched_proof: DeepFoldBatchedProof<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct LigeSISProverCommitmentAdvice<F: PrimeField> {
    pub mat_f_prime: Vec<Vec<F>>,
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

    fn gen_srs_for_testing<R: Rng>(rng: &mut R, log_size: usize) -> Result<Self::SRS, PCSError> {
        let eta = F::ONE.into_bigint().to_bits_be().len();
        let lambda = 128usize;
        let mu = log_size;
        let log_m = if log_size < 8 { 0 } else { (log_size - 8) / 2 };
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

        let deepfold_srs = DeepFoldPCS::<F>::gen_srs_for_testing(
            rng,
            max(
                max(log_c, log_s_lambda) + log_m + log_eta,
                log_c + 1 + log_n,
            ),
        )?;
        Ok(LigeSISSRS {
            eta,
            lambda,
            mu,
            log_m,
            rs_len,
            c,
            mat_a,
            deepfold_srs,
        })
    }

    fn setup(
        srs: impl Borrow<Self::SRS>,
        _supported_degree: Option<usize>,
        _supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let LigeSISSRS {
            eta,
            lambda,
            mu,
            log_m,
            rs_len,
            c,
            mat_a,
            deepfold_srs,
        } = srs.borrow().clone();
        let log_n = mu - log_m;
        let n = 1 << log_n;
        let s_lambda = min(lambda, rs_len);
        println!("deepfold mu = {}", deepfold_srs.max_mu);
        let (deepfold_prover_param, deepfold_verifier_param) = DeepFoldPCS::<F>::setup(
            deepfold_srs,
            Some(deepfold_srs.max_mu.clone()),
            Some(deepfold_srs.max_mu.clone()),
        )?;

        let mut transcript = IOPTranscript::new(b"setup");
        let mut transcript_clone = transcript.clone();
        let (com_mat_a, com_mat_a_advice) = DeepFoldPCS::commit(
            &deepfold_prover_param,
            &evals_to_arcpoly(&mat_a.concat()),
            &mut transcript,
        )?;

        let rs = ReedSolomon::<F>::new(n, rs_len);
        let g = rs.get_generator();

        let com_mat_a_v_advice = DeepFoldPCS::verifier_receive_commit(
            &deepfold_verifier_param,
            &com_mat_a,
            &mut transcript_clone,
        )?;
        let prover_param = LigeSISProverParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            rs,
            mat_a,
            com_mat_a_advice,
            deepfold_prover_param,
        };
        let verifier_param = LigeSISVerifierParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            g,
            com_mat_a,
            com_mat_a_v_advice,
            deepfold_verifier_param,
        };
        Ok((prover_param, verifier_param))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        // trim parameters
        let &LigeSISProverParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            ref rs,
            ref mat_a,
            ref com_mat_a_advice,
            ref deepfold_prover_param,
        } = prover_param.borrow();
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_f = reshape(&poly.evaluations, m, n);

        // encode `F`
        let start = std::time::Instant::now();
        let mat_f_prime = mat_f.iter().map(|row| rs.encode(row)).collect::<Vec<_>>();
        println!("RS(F): {} s", start.elapsed().as_secs_f64());

        // compute `H`
        let start = std::time::Instant::now();
        let mat_a_prime = mat_a
            .iter()
            .map(|row| {
                (0..eta * m / 8)
                    .map(|i| get_mat_a_byte_bucket(&row[i * 8..(i + 1) * 8].to_vec()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mat_h = mat_a_prime
            .iter()
            .map(|row| {
                (0..n * 2)
                    .map(|j| {
                        (0..m)
                            .map(|i| {
                                mat_f_prime[i][j]
                                    .into_bigint()
                                    .to_bytes_le()
                                    .iter()
                                    .enumerate()
                                    .map(|(k, x)| row[eta * i / 8 + k][*x as usize])
                                    .sum::<F>()
                            })
                            .sum::<F>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        println!("SIS_Hash(B): {} s", start.elapsed().as_secs_f64());

        // compute com(H)
        let (com_mat_h, com_mat_h_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&mat_h.concat()),
            &com_mat_a_advice.alpha0,
            transcript,
        )?;

        Ok((
            LigeSISCommitment { com_mat_h },
            LigeSISProverCommitmentAdvice {
                mat_f_prime,
                mat_h,
                com_mat_h_advice,
            },
        ))
    }

    fn d_commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<(Option<Self::Commitment>, Self::ProverCommitmentAdvice), PCSError> {
        // trim parameters
        let num_party = Net::n_parties();
        let num_party_vars = Net::n_parties().log_2() as usize;

        let &LigeSISProverParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            ref rs,
            ref mat_a,
            ref com_mat_a_advice,
            ref deepfold_prover_param,
        } = prover_param.borrow();
        let log_m = log_m - num_party_vars;
        let (m, n) = (1 << log_m, 1 << log_n);
        let mat_f = reshape(&poly.evaluations, m, n);
        // encode `F`
        let mat_f_prime = mat_f.iter().map(|row| rs.encode(row)).collect::<Vec<_>>();

        // decompose `F`
        let mat_b_trans = transposition(&mat_f_prime)
            .iter()
            .map(|col| decompose_vector(col))
            .collect::<Vec<_>>();
        let mat_h_i = field_mat_mul_trans_bool_mat(&mat_a, &mat_b_trans);

        let all_mat_h = Net::send_to_master(&mat_h_i);

        if Net::am_master() {
            let all_mat_h = all_mat_h.unwrap();
            let mat_h = (0..c)
                .map(|i| {
                    (0..2 * n)
                        .map(|j| (0..num_party).map(|k| all_mat_h[k][i][j]).sum::<F>())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let (com_mat_h, com_mat_h_advice) = DeepFoldPCS::commit(
                deepfold_prover_param,
                &evals_to_arcpoly(&mat_h.concat()),
                transcript,
            )?;
            Ok((
                Some(LigeSISCommitment { com_mat_h }),
                LigeSISProverCommitmentAdvice {
                    mat_f_prime,
                    mat_h,
                    com_mat_h_advice,
                },
            ))
        } else {
            Ok((
                None,
                LigeSISProverCommitmentAdvice {
                    mat_f_prime,
                    mat_h: mat_h_i,
                    com_mat_h_advice: DeepFoldProverCommitmentAdvice::default(),
                },
            ))
        }
    }

    fn verifier_receive_commit(
        verifier_param: &Self::VerifierParam,
        commitment: &Self::Commitment,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::VerifierCommitmentAdvice, PCSError> {
        let com_mat_h_v_advice = DeepFoldPCS::verifier_receive_commit_given_alpha(
            &verifier_param.deepfold_verifier_param,
            &commitment.com_mat_h,
            &verifier_param.com_mat_a_v_advice.alpha0,
            transcript,
        )?;
        Ok(LigeSISVerifierCommitmentAdvice { com_mat_h_v_advice })
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        let &LigeSISProverParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            ref rs,
            ref mat_a,
            ref com_mat_a_advice,
            ref deepfold_prover_param,
        } = prover_param.borrow();
        let (m, n) = (1 << log_m, 1 << log_n);
        let log_rs_len = rs_len.ilog2() as usize;
        let deepfold_alpha0 = com_mat_a_advice.alpha0;

        assert_eq!(mu, log_m + log_n);
        assert_eq!(poly.num_vars, mu);

        let mat_f = reshape(&poly.evaluations, m, n);

        let LigeSISProverCommitmentAdvice {
            mat_f_prime,
            mat_h,
            com_mat_h_advice,
        } = advice;

        // Step 1
        let (z1, z2): (Vec<F>, Vec<F>) = (point[log_n..].to_vec(), point[..log_n].to_vec());
        let eq_z1 = get_tensor(&z1);

        // Step 2
        let start = std::time::Instant::now();
        let a = (0..n)
            .map(|j| (0..m).map(|i| eq_z1[i] * mat_f[i][j]).sum())
            .collect::<Vec<F>>();
        let (com_a, com_a_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&a),
            &deepfold_alpha0,
            transcript,
        )?;
        println!("DeepFold.Commit(a): {} s", start.elapsed().as_secs_f64());

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;

        // Step 4
        // let mat_b_trans = transposition(&mat_b);
        let start = std::time::Instant::now();
        let mat_f_prime_trans = transposition(&mat_f_prime);
        let mat_bI = transposition(
            &I.iter()
                .map(|&i| decompose_vector(&mat_f_prime_trans[i]))
                .collect::<Vec<_>>(),
        );
        let bI_field = bool_vec_to_field_vec(&mat_bI.concat());
        let (com_bI, com_bI_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&bI_field),
            &deepfold_alpha0,
            transcript,
        )?;
        println!("DeepFold.Commit(B_I): {} s", start.elapsed().as_secs_f64());

        // Step 5
        let alpha1 = transcript
            .get_and_append_challenge_vectors(b"alpha1", (m * eta * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
        let alpha3 = transcript.get_and_append_challenge_vectors(b"alpha3", log_rs_len)?;

        // Step 6
        let start = std::time::Instant::now();
        let mut bI_check = VirtualPolynomial::new(bI_field.len().ilog2() as usize);
        bI_check
            .add_mle_list(
                [
                    evals_to_arcpoly(&bI_field),
                    evals_to_arcpoly(&bI_field.iter().map(|&x| x - F::ONE).collect::<Vec<F>>()),
                    evals_to_arcpoly(&get_tensor(&alpha1)),
                ],
                F::ONE,
            )
            .unwrap();
        let bI_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(bI_check, transcript).unwrap();
        let r1 = bI_check_proof.point.clone();

        // Step 7 Check rs_a
        let rs_a = rs.encode(&a);
        let g = rs.get_generator();
        let (com_rs_a, com_rs_a_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&rs_a),
            &deepfold_alpha0,
            transcript,
        )?;

        // Step 7.1 check eq_alpha3^T * G * a = eq_alpha3^T * rs_a
        //  \sum_i alpha3_mat_g(i) * a(i) = eq_alpha3^T * rs_a
        // reduce to alpha3_mat_g(r6), a(r6)
        let alpha3_mat_g = compute_alpha_mat_g(log_rs_len as usize, log_n, &g, &alpha3);
        let mut rs_a_check = VirtualPolynomial::new(log_n);
        rs_a_check
            .add_mle_list(
                [
                    evals_to_arcpoly(&alpha3_mat_g[log_rs_len][..n].to_vec()),
                    evals_to_arcpoly(&a),
                ],
                F::ONE,
            )
            .unwrap();
        let rs_a_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(rs_a_check, transcript).unwrap();
        let r6 = rs_a_check_proof.point.clone();

        // Step 7.2 check alpha3_mat_g(r6)
        let mut cur_p = vec![r6.clone(), vec![F::ZERO; log_rs_len - log_n]].concat();
        let mut mat_g_check_proofs = Vec::new();
        for i in (2..=log_rs_len).rev() {
            let (x, b) = (cur_p[..i - 1].to_vec(), cur_p[i - 1]);
            let gi = g.pow([1u64 << (log_rs_len - i)]);
            let w = (0..1 << (i - 1))
                .map(|z| {
                    F::ONE - alpha3[log_rs_len - i]
                        + alpha3[log_rs_len - i]
                            * (gi.pow([z]) * (F::ONE - b) + gi.pow([z + (1 << (i - 1))]) * b)
                })
                .collect::<Vec<_>>();
            let mut mat_g_check = VirtualPolynomial::new(i - 1);
            mat_g_check
                .add_mle_list(
                    [
                        evals_to_arcpoly(&get_tensor(&x)),
                        evals_to_arcpoly(&alpha3_mat_g[i - 1]),
                        evals_to_arcpoly(&w),
                    ],
                    F::ONE,
                )
                .unwrap();
            let mat_g_check_proof =
                <PolyIOP<F> as SumCheck<F>>::prove(mat_g_check, transcript).unwrap();

            cur_p = mat_g_check_proof.point.clone();
            mat_g_check_proofs.push(mat_g_check_proof);
        }

        // Step 8 Lookup Argument
        let eq_alpha2_a_bI = mat_mul(
            &vec![get_tensor(&alpha2)],
            &field_mat_mul_bool_mat(&mat_a, &mat_bI),
        )
        .concat();
        let v = otimes(
            &get_tensor(&z1),
            &(0..eta)
                .map(|i| F::from(2u64).pow([i as u64]))
                .collect::<Vec<_>>(),
        );
        let v_bI = field_mat_mul_bool_mat(&vec![v.clone()], &mat_bI).concat();
        let eq_alpha2_h = mat_mul(&vec![get_tensor(&alpha2)], &mat_h).concat();

        // Lookup Argument for (I, eq_alpha2_a_bI, v_bI) in ([2n], eq_alpha2_h, rs_a)
        assert!(
            (0..s_lambda).all(|k| eq_alpha2_a_bI[k] == eq_alpha2_h[I[k]] && v_bI[k] == rs_a[I[k]])
        );
        let lookup_proof = ();
        let r2 = vec![F::ONE; s_lambda.ilog2() as usize];
        let r3 = vec![F::ONE; 1 + log_n];

        // Step 9
        let alpha2_a = mat_mul(&vec![get_tensor(&alpha2)], &mat_a)[0].clone();
        let bI_r2 =
            field_mat_mul_bool_mat(&vec![get_tensor(&r2)], &transposition(&mat_bI))[0].clone();
        let mut alpha2_a_bI_r2_check = VirtualPolynomial::new(mat_bI.len().ilog2() as usize);
        alpha2_a_bI_r2_check
            .add_mle_list(
                [evals_to_arcpoly(&alpha2_a), evals_to_arcpoly(&bI_r2)],
                F::ONE,
            )
            .unwrap();
        let alpha2_a_bI_r2_check_proof =
            <PolyIOP<F> as SumCheck<F>>::prove(alpha2_a_bI_r2_check, transcript).unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();

        // Step 10
        let mut v_bI_r2_check = VirtualPolynomial::new(v.len().ilog2() as usize);
        v_bI_r2_check
            .add_mle_list([evals_to_arcpoly(&v), evals_to_arcpoly(&bI_r2)], F::ONE)
            .unwrap();
        let v_bI_r2_check_proof =
            <PolyIOP<F> as SumCheck<F>>::prove(v_bI_r2_check, transcript).unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();
        println!("Sumchecks: {} s", start.elapsed().as_secs_f64());

        // Step 11
        let start = std::time::Instant::now();
        let polys = [
            &a,
            &a,
            &rs_a,
            &rs_a,
            &mat_h.concat(),
            &mat_a.concat(),
            &bI_field,
            &bI_field,
            &bI_field,
        ]
        .map(|p| evals_to_arcpoly(p))
        .to_vec();
        let advices = [
            &com_a_advice,
            &com_a_advice,
            &com_rs_a_advice,
            &com_rs_a_advice,
            &com_mat_h_advice,
            &com_mat_a_advice,
            &com_bI_advice,
            &com_bI_advice,
            &com_bI_advice,
        ]
        .map(|a| a.clone());
        let points = [
            &z2,
            &r6,
            &r3,
            &alpha3,
            &vec![r3.clone(), alpha2.clone()].concat(),
            &vec![r4.clone(), alpha2.clone()].concat(),
            &r1,
            &vec![r2.clone(), r4.clone()].concat(),
            &vec![r2.clone(), r5.clone()].concat(),
        ]
        .map(|p| p.clone());
        let evals = [F::ZERO; 9];

        let deepfold_batched_proof = DeepFoldPCS::multi_open(
            deepfold_prover_param,
            polys,
            &advices,
            &points,
            &evals,
            transcript,
        )?;

        println!("DeepFold.Open: {} s", start.elapsed().as_secs_f64());

        Ok(LigeSISProof {
            // commitments
            com_a,
            com_bI,
            com_rs_a,
            // sumcheck proofs
            bI_check_proof,
            alpha2_a_bI_r2_check_proof,
            v_bI_r2_check_proof,
            rs_a_check_proof,
            mat_g_check_proofs,
            // lookup proof
            lookup_proof,
            // commitment proofs
            deepfold_batched_proof,
        })
    }

    fn d_open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Option<Self::Proof>, PCSError> {
        let &LigeSISProverParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            ref rs,
            ref mat_a,
            ref com_mat_a_advice,
            ref deepfold_prover_param,
        } = prover_param.borrow();
        let (m, n) = (1 << log_m, 1 << log_n);
        let log_rs_len = rs_len.ilog2() as usize;
        let deepfold_alpha0 = com_mat_a_advice.alpha0;

        assert_eq!(mu, log_m + log_n);
        assert_eq!(poly.num_vars, mu);

        let mat_f = reshape(&poly.evaluations, m, n);

        let LigeSISProverCommitmentAdvice {
            mat_f_prime,
            mat_h,
            com_mat_h_advice,
        } = advice;

        // Step 1
        let (z1, z2): (Vec<F>, Vec<F>) = (point[log_n..].to_vec(), point[..log_n].to_vec());
        let eq_z1 = get_tensor(&z1);

        // Step 2
        let start = std::time::Instant::now();
        let a = (0..n)
            .map(|j| (0..m).map(|i| eq_z1[i] * mat_f[i][j]).sum())
            .collect::<Vec<F>>();
        let (com_a, com_a_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&a),
            &deepfold_alpha0,
            transcript,
        )?;
        println!("DeepFold.Commit(a): {} s", start.elapsed().as_secs_f64());

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;

        // Step 4
        // let mat_b_trans = transposition(&mat_b);
        let start = std::time::Instant::now();
        let mat_f_prime_trans = transposition(&mat_f_prime);
        let mat_bI = transposition(
            &I.iter()
                .map(|&i| decompose_vector(&mat_f_prime_trans[i]))
                .collect::<Vec<_>>(),
        );
        let bI_field = bool_vec_to_field_vec(&mat_bI.concat());
        let (com_bI, com_bI_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&bI_field),
            &deepfold_alpha0,
            transcript,
        )?;
        println!("DeepFold.Commit(B_I): {} s", start.elapsed().as_secs_f64());

        // Step 5
        let alpha1 = transcript
            .get_and_append_challenge_vectors(b"alpha1", (m * eta * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
        let alpha3 = transcript.get_and_append_challenge_vectors(b"alpha3", log_rs_len)?;

        // Step 6
        let start = std::time::Instant::now();
        let mut bI_check = VirtualPolynomial::new(bI_field.len().ilog2() as usize);
        bI_check
            .add_mle_list(
                [
                    evals_to_arcpoly(&bI_field),
                    evals_to_arcpoly(&bI_field.iter().map(|&x| x - F::ONE).collect::<Vec<F>>()),
                    evals_to_arcpoly(&get_tensor(&alpha1)),
                ],
                F::ONE,
            )
            .unwrap();
        let bI_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(bI_check, transcript).unwrap();
        let r1 = bI_check_proof.point.clone();

        // Step 7 Check rs_a
        let rs_a = rs.encode(&a);
        let g = rs.get_generator();
        let (com_rs_a, com_rs_a_advice) = DeepFoldPCS::commit_given_alpha(
            deepfold_prover_param,
            &evals_to_arcpoly(&rs_a),
            &deepfold_alpha0,
            transcript,
        )?;

        // Step 7.1 check eq_alpha3^T * G * a = eq_alpha3^T * rs_a
        //  \sum_i alpha3_mat_g(i) * a(i) = eq_alpha3^T * rs_a
        // reduce to alpha3_mat_g(r6), a(r6)
        let alpha3_mat_g = compute_alpha_mat_g(log_rs_len as usize, log_n, &g, &alpha3);
        let mut rs_a_check = VirtualPolynomial::new(log_n);
        rs_a_check
            .add_mle_list(
                [
                    evals_to_arcpoly(&alpha3_mat_g[log_rs_len][..n].to_vec()),
                    evals_to_arcpoly(&a),
                ],
                F::ONE,
            )
            .unwrap();
        let rs_a_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(rs_a_check, transcript).unwrap();
        let r6 = rs_a_check_proof.point.clone();

        // Step 7.2 check alpha3_mat_g(r6)
        let mut cur_p = vec![r6.clone(), vec![F::ZERO; log_rs_len - log_n]].concat();
        let mut mat_g_check_proofs = Vec::new();
        for i in (2..=log_rs_len).rev() {
            let (x, b) = (cur_p[..i - 1].to_vec(), cur_p[i - 1]);
            let gi = g.pow([1u64 << (log_rs_len - i)]);
            let w = (0..1 << (i - 1))
                .map(|z| {
                    F::ONE - alpha3[log_rs_len - i]
                        + alpha3[log_rs_len - i]
                            * (gi.pow([z]) * (F::ONE - b) + gi.pow([z + (1 << (i - 1))]) * b)
                })
                .collect::<Vec<_>>();
            let mut mat_g_check = VirtualPolynomial::new(i - 1);
            mat_g_check
                .add_mle_list(
                    [
                        evals_to_arcpoly(&get_tensor(&x)),
                        evals_to_arcpoly(&alpha3_mat_g[i - 1]),
                        evals_to_arcpoly(&w),
                    ],
                    F::ONE,
                )
                .unwrap();
            let mat_g_check_proof =
                <PolyIOP<F> as SumCheck<F>>::prove(mat_g_check, transcript).unwrap();

            cur_p = mat_g_check_proof.point.clone();
            mat_g_check_proofs.push(mat_g_check_proof);
        }

        // Step 8 Lookup Argument
        let eq_alpha2_a_bI = mat_mul(
            &vec![get_tensor(&alpha2)],
            &field_mat_mul_bool_mat(&mat_a, &mat_bI),
        )
        .concat();
        let v = otimes(
            &get_tensor(&z1),
            &(0..eta)
                .map(|i| F::from(2u64).pow([i as u64]))
                .collect::<Vec<_>>(),
        );
        let v_bI = field_mat_mul_bool_mat(&vec![v.clone()], &mat_bI).concat();
        let eq_alpha2_h = mat_mul(&vec![get_tensor(&alpha2)], &mat_h).concat();

        // Lookup Argument for (I, eq_alpha2_a_bI, v_bI) in ([2n], eq_alpha2_h, rs_a)
        assert!(
            (0..s_lambda).all(|k| eq_alpha2_a_bI[k] == eq_alpha2_h[I[k]] && v_bI[k] == rs_a[I[k]])
        );
        let lookup_proof = ();
        let r2 = vec![F::ONE; s_lambda.ilog2() as usize];
        let r3 = vec![F::ONE; 1 + log_n];

        // Step 9
        let alpha2_a = mat_mul(&vec![get_tensor(&alpha2)], &mat_a)[0].clone();
        let bI_r2 =
            field_mat_mul_bool_mat(&vec![get_tensor(&r2)], &transposition(&mat_bI))[0].clone();
        let mut alpha2_a_bI_r2_check = VirtualPolynomial::new(mat_bI.len().ilog2() as usize);
        alpha2_a_bI_r2_check
            .add_mle_list(
                [evals_to_arcpoly(&alpha2_a), evals_to_arcpoly(&bI_r2)],
                F::ONE,
            )
            .unwrap();
        let alpha2_a_bI_r2_check_proof =
            <PolyIOP<F> as SumCheck<F>>::prove(alpha2_a_bI_r2_check, transcript).unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();

        // Step 10
        let mut v_bI_r2_check = VirtualPolynomial::new(v.len().ilog2() as usize);
        v_bI_r2_check
            .add_mle_list([evals_to_arcpoly(&v), evals_to_arcpoly(&bI_r2)], F::ONE)
            .unwrap();
        let v_bI_r2_check_proof =
            <PolyIOP<F> as SumCheck<F>>::prove(v_bI_r2_check, transcript).unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();
        println!("Sumchecks: {} s", start.elapsed().as_secs_f64());

        // Step 11
        let start = std::time::Instant::now();
        let polys = [
            &a,
            &a,
            &rs_a,
            &rs_a,
            &mat_h.concat(),
            &mat_a.concat(),
            &bI_field,
            &bI_field,
            &bI_field,
        ]
        .map(|p| evals_to_arcpoly(p))
        .to_vec();
        let advices = [
            &com_a_advice,
            &com_a_advice,
            &com_rs_a_advice,
            &com_rs_a_advice,
            &com_mat_h_advice,
            &com_mat_a_advice,
            &com_bI_advice,
            &com_bI_advice,
            &com_bI_advice,
        ]
        .map(|a| a.clone());
        let points = [
            &z2,
            &r6,
            &r3,
            &alpha3,
            &vec![r3.clone(), alpha2.clone()].concat(),
            &vec![r4.clone(), alpha2.clone()].concat(),
            &r1,
            &vec![r2.clone(), r4.clone()].concat(),
            &vec![r2.clone(), r5.clone()].concat(),
        ]
        .map(|p| p.clone());
        let evals = [F::ZERO; 9];

        let deepfold_batched_proof = DeepFoldPCS::multi_open(
            deepfold_prover_param,
            polys,
            &advices,
            &points,
            &evals,
            transcript,
        )?;

        println!("DeepFold.Open: {} s", start.elapsed().as_secs_f64());

        if Net::am_master() {
            Ok(Some(LigeSISProof {
                com_a,
                com_bI,
                com_rs_a,
                bI_check_proof,
                alpha2_a_bI_r2_check_proof,
                v_bI_r2_check_proof,
                rs_a_check_proof,
                mat_g_check_proofs,
                lookup_proof,
                deepfold_batched_proof,
            }))
        } else {
            Ok(None)
        }
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
        let LigeSISVerifierParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            rs_len,
            c,
            g,
            com_mat_a,
            com_mat_a_v_advice,
            deepfold_verifier_param,
        } = verifier_param.borrow().clone();
        let LigeSISCommitment { com_mat_h } = com.clone();
        let LigeSISVerifierCommitmentAdvice { com_mat_h_v_advice } = advice.clone();
        let (m, n) = (1 << log_m, 1 << log_n);
        let log_rs_len = rs_len.ilog2() as usize;
        let deepfold_alpha0 = com_mat_a_v_advice.alpha0;
        let LigeSISProof {
            com_a,
            com_bI,
            com_rs_a,
            bI_check_proof,
            alpha2_a_bI_r2_check_proof,
            v_bI_r2_check_proof,
            rs_a_check_proof,
            mat_g_check_proofs,
            lookup_proof,
            deepfold_batched_proof,
        } = proof.clone();
        let values = deepfold_batched_proof.evals.clone();

        // Step 1
        let (z1, z2): (Vec<F>, Vec<F>) = (point[log_n..].to_vec(), point[..log_n].to_vec());

        // Step 2
        let com_a_v_advice = DeepFoldPCS::verifier_receive_commit_given_alpha(
            &deepfold_verifier_param,
            &com_a,
            &deepfold_alpha0,
            transcript,
        )?;

        // Step 3
        let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;

        // Step 4
        let com_bI_v_advice = DeepFoldPCS::verifier_receive_commit_given_alpha(
            &deepfold_verifier_param,
            &com_bI,
            &deepfold_alpha0,
            transcript,
        )?;

        // Step 5
        let alpha1 = transcript
            .get_and_append_challenge_vectors(b"alpha1", (m * eta * s_lambda).ilog2() as usize)?;
        let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
        let alpha3 =
            transcript.get_and_append_challenge_vectors(b"alpha3", rs_len.ilog2() as usize)?;
        // println!("verify : {}", alpha2[0]);

        // Step 6
        let bI_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&bI_check_proof);
        let bI_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
            bI_check_sum,
            &bI_check_proof,
            &VPAuxInfo {
                max_degree: 3,
                num_variables: (m * eta * s_lambda).ilog2() as usize,
                phantom: PhantomData::<F>::default(),
            },
            transcript,
        )
        .unwrap();
        let r1 = bI_check_proof.point.clone();
        let bI_r1 = values[6];
        if bI_check_claim.expected_evaluation
            != bI_r1 * (bI_r1 - F::ONE) * eval_mle_eq(&alpha1, &r1)
        {
            return Ok(false);
        }

        // Step 7
        let com_rs_a_v_advice = DeepFoldPCS::verifier_receive_commit_given_alpha(
            &deepfold_verifier_param,
            &com_rs_a,
            &deepfold_alpha0,
            transcript,
        )?;

        // Step 7.1
        let rs_a_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&rs_a_check_proof);
        let rs_a_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
            rs_a_check_sum,
            &rs_a_check_proof,
            &VPAuxInfo {
                max_degree: 2,
                num_variables: log_n,
                phantom: PhantomData::<F>::default(),
            },
            transcript,
        )
        .unwrap();
        let r6 = rs_a_check_proof.point.clone();
        let a_r6 = values[1];
        if rs_a_check_sum != values[3] {
            return Ok(false);
        }

        // Step 7.2
        let mut cur_p = vec![r6.clone(), vec![F::ZERO; log_rs_len - log_n]].concat();
        let mut expected_eval = rs_a_check_claim.expected_evaluation / a_r6;
        for i in (2..=log_rs_len).rev() {
            let (x, b) = (cur_p[..i - 1].to_vec(), cur_p[i - 1]);
            let mat_g_check_proof = &mat_g_check_proofs[log_rs_len - i];
            let mat_g_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(mat_g_check_proof);
            let mat_g_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
                mat_g_check_sum,
                mat_g_check_proof,
                &VPAuxInfo {
                    max_degree: 3,
                    num_variables: i - 1,
                    phantom: PhantomData::<F>::default(),
                },
                transcript,
            )
            .unwrap();
            let mut r = mat_g_check_proof.point.clone();
            r.push(b);
            let v = (0..i)
                .map(|k| F::ONE - r[k] + r[k] * g.pow([rs_len as u64 >> (i - k)]))
                .product::<F>();
            expected_eval = mat_g_check_claim.expected_evaluation
                / eval_mle_eq(&mat_g_check_proof.point, &x)
                / (F::ONE - alpha3[log_rs_len - i] + alpha3[log_rs_len - i] * v);
            cur_p = mat_g_check_proof.point.clone();
        }
        if expected_eval
            != F::ONE - alpha3[log_rs_len - 1]
                + alpha3[log_rs_len - 1]
                    * (F::ONE - cur_p[0] + cur_p[0] * g.pow([rs_len as u64 >> 1]))
        {
            return Ok(false);
        }

        // Step 8
        let r2 = vec![F::ONE; s_lambda.ilog2() as usize];
        let r3 = vec![F::ONE; 1 + log_n];

        // Step 9
        let alpha2_a_bI_r2_check_sum =
            <PolyIOP<F> as SumCheck<F>>::extract_sum(&alpha2_a_bI_r2_check_proof);
        let alpha2_a_bI_r2_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
            alpha2_a_bI_r2_check_sum,
            &alpha2_a_bI_r2_check_proof,
            &VPAuxInfo {
                max_degree: 2,
                num_variables: (m * eta).ilog2() as usize,
                phantom: PhantomData::<F>::default(),
            },
            transcript,
        )
        .unwrap();
        let r4 = alpha2_a_bI_r2_check_proof.point.clone();
        let bI_r4_r2 = values[7];
        let mat_a_alpha2_r_4 = values[5];
        if alpha2_a_bI_r2_check_claim.expected_evaluation != bI_r4_r2 * mat_a_alpha2_r_4 {
            return Ok(false);
        }

        // Step 10
        let v_bI_r2_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&v_bI_r2_check_proof);
        let v_bI_r2_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
            v_bI_r2_check_sum,
            &v_bI_r2_check_proof,
            &VPAuxInfo {
                max_degree: 2,
                num_variables: (m * eta).ilog2() as usize,
                phantom: PhantomData::<F>::default(),
            },
            transcript,
        )
        .unwrap();
        let r5 = v_bI_r2_check_proof.point.clone();
        let bI_r5_r2 = values[8];
        let (r5_0, r5_1) = (
            r5[..eta.ilog2() as usize].to_vec(),
            r5[eta.ilog2() as usize..].to_vec(),
        );
        let t_r5_0 = get_tensor(&r5_0);
        let powers = {
            let mut res = vec![F::ONE];
            for i in 1..eta {
                res.push(res[i - 1] + res[i - 1])
            }
            res
        };
        let r5_0_p = (0..eta).map(|i| t_r5_0[i] * powers[i]).sum::<F>();
        if v_bI_r2_check_claim.expected_evaluation != bI_r5_r2 * eval_mle_eq(&r5_1, &z1) * r5_0_p {
            return Ok(false);
        }

        // Step 11
        let coms = [
            &com_a, &com_a, &com_rs_a, &com_rs_a, &com_mat_h, &com_mat_a, &com_bI, &com_bI, &com_bI,
        ]
        .map(|c| c.clone());
        let points = [
            z2,
            r6,
            r3.clone(),
            alpha3.clone(),
            vec![r3.clone(), alpha2.clone()].concat(),
            vec![r4.clone(), alpha2.clone()].concat(),
            r1,
            vec![r2.clone(), r4.clone()].concat(),
            vec![r2.clone(), r5.clone()].concat(),
        ];
        let advices = [
            &com_a_v_advice,
            &com_a_v_advice,
            &com_rs_a_v_advice,
            &com_rs_a_v_advice,
            &com_mat_h_v_advice,
            &com_mat_a_v_advice,
            &com_bI_v_advice,
            &com_bI_v_advice,
            &com_bI_v_advice,
        ]
        .map(|a| a.clone());
        if values[0] != *value {
            return Ok(false);
        }
        if !DeepFoldPCS::batch_verify(
            &deepfold_verifier_param,
            &coms,
            &points,
            &advices,
            &deepfold_batched_proof,
            transcript,
        )? {
            return Ok(false);
        }

        return Ok(true);
    }
}
