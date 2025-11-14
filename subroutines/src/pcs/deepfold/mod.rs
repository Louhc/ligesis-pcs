use core::num;
use std::time::Instant;

use crate::{pcs::{deepfold, prelude::*}, IOPProof, PolyIOP, SumCheck};
use arithmetic::{VirtualPolynomial, VPAuxInfo};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, marker::PhantomData, rand::Rng,
    sync::Arc, vec, vec::Vec, cmp::max,
    collections::HashSet,
};
use itertools::concat;
use transcript::IOPTranscript;

#[cfg(test)]
mod tests;

/// DeepFold Polynomial Commitment Scheme
pub struct DeepFoldPCS<F: PrimeField> {
    #[doc(hidden)]
    phantom: PhantomData<F>,
}


#[derive(Clone, Debug, Copy)]
pub struct DeepFoldSRS<F: PrimeField> {
    pub max_mu: usize,
    pub l0: GeneralEvaluationDomain<F>,
    pub s: usize,
}

#[derive(Clone)]
pub struct DeepFoldProverParam<F: PrimeField> {
    pub max_mu: usize,
    pub l0: GeneralEvaluationDomain<F>,
    pub s: usize,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct DeepFoldVerifierParam<F: PrimeField> {
    pub max_mu: usize,
    pub len_l0: usize,
    pub g: F,
    pub s: usize,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct DeepFoldProof<F: PrimeField> {
    pub linear_polys: Vec<Vec<(F, F)>>,
    pub mt_roots: Vec<Byte32>,
    pub f_mu: F,
    pub mt_proofs: Vec<Vec<(F, Vec<Byte32>, usize)>>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
pub struct DeepFoldBatchedProof<F: PrimeField> {
    pub deepfold_proof: DeepFoldProof<F>,
    pub sum_check_proof: IOPProof<F>,
    pub mt_proofs_for_mt0: Vec<Vec<(F, Vec<Byte32>, usize)>>,
    pub evals: Vec<F>,
    pub sum_check_evals: Vec<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct DeepFoldProverCommitmentAdvice<F: PrimeField> {
    pub f0: Vec<F>,
    pub mt0: MerkleTree,
    pub v0: Vec<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct DeepFoldCommitment {
    pub mu: usize,
    pub rt0: Byte32,
}


impl<F: PrimeField> DeepFoldPCS<F> {
    pub fn compute_value_from_proof(
        point: &Vec<F>,
        proof: &DeepFoldProof<F>,
    ) -> F {
        eval_linear_poly(&proof.linear_polys[0][0], &point[0])
    }
}

impl<F: PrimeField> PolynomialCommitmentScheme<F> for DeepFoldPCS<F> {
    // Parameters
    type ProverParam = DeepFoldProverParam<F>; 
    type VerifierParam = DeepFoldVerifierParam<F>;
    type SRS = DeepFoldSRS<F>;
    // Polynomial and its associated types
    type Polynomial = Arc<DenseMultilinearExtension<F>>;
    type ProverCommitmentAdvice = DeepFoldProverCommitmentAdvice<F>;
    type Point = Vec<F>;
    type Evaluation = F;
    // Commitments and proofs
    type Commitment = DeepFoldCommitment; // merkle tree root
    type Proof = DeepFoldProof<F>; // merkle tree paths, columes of `E`
    type BatchProof = DeepFoldBatchedProof<F>; // 

    fn gen_srs_for_testing<R: Rng>(
        _rng: &mut R, 
        log_size: usize
    ) -> Result<Self::SRS, PCSError> {
        let max_mu = log_size;
        let len_l0 = (1 << max_mu) * 2;
        let l0 = GeneralEvaluationDomain::<F>::new(len_l0).unwrap();
        let s = 10;
        Ok(DeepFoldSRS{max_mu, l0, s})
    }

    fn setup(
        srs: impl Borrow<Self::SRS>,
        _supported_degree: Option<usize>,
        _supported_num_vars: Option<usize>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let srs = srs.borrow();
        Ok((
            DeepFoldProverParam{max_mu: srs.max_mu, l0: srs.l0, s: srs.s},
            DeepFoldVerifierParam{max_mu: srs.max_mu, len_l0: srs.l0.size(), g: srs.l0.element(1), s: srs.s},
        ))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        let &Self::ProverParam{max_mu, l0, s} = prover_param.borrow();
        // let mu = poly.num_vars;
        let mu = max_mu;
        let poly = resize_poly(&poly, mu);
        assert!(mu <= max_mu);

        let f0 = evals_to_coeffs(mu, &poly.evaluations);
        let v0 = l0.fft(&f0);

        let mt0 = MerkleTree::new(&v0.iter().map(|&x| compute_sha256_row(&[x])).collect());

        let rt0 = mt0.root();
        Ok((
            DeepFoldCommitment{mu, rt0},
            DeepFoldProverCommitmentAdvice{f0, mt0, v0},
        ))
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        let &Self::ProverParam{max_mu, l0, s} = prover_param.borrow();
        
        // let mu = poly.num_vars;
        let mu = max_mu;
        let poly = resize_poly(&poly, mu);
        let point = resize_point(&point, mu);
        
        assert!(mu <= max_mu);

        let Self::ProverCommitmentAdvice{f0, mt0, v0} = advice.clone();
        let mut a = vec![Vec::new()];
        let mut f_tilde = vec![poly.evaluations.clone()];
        let mut f = vec![f0];
        let mut alpha = vec![F::ZERO];
        let mut linear_polys = Vec::new();
        let mut l = vec![l0];
        l.append(&mut (1..mu + 1)
            .map(|i| GeneralEvaluationDomain::<F>::new(l0.size() >> i).unwrap())
            .collect::<Vec<_>>());
        let mut v = vec![v0];
        let mut mt_roots = vec![mt0.root().clone()];
        let mut mt = vec![mt0];
        let mut mt_proofs = Vec::new();
        let mut f_mu = F::ZERO;
        let mut r = vec![F::ZERO];

        // Step 1
        a[0].push(point.clone());
        
        // Step 2
        for i in 1..mu + 1 {
            // Step 2.a
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);
            a[i - 1].push(get_alpha_powers::<F>(alpha[i], mu - i + 1));
            let (f0, f1) = split_even_odd(&f_tilde[i - 1]);
            let (fe, fo) = split_even_odd(&f[i - 1]);
            // Step 2.b
            if i == mu {
                linear_polys.push(vec![(f_tilde[i - 1][0], f_tilde[i - 1][1])]);
            } else {
                linear_polys.push(
                    a[i - 1].iter().map(
                        |w| {
                            assert!(!w.is_empty());
                            let w_tensor = get_tensor(&w[1..].to_vec());
                            (inner_product(&w_tensor, &f0), inner_product(&w_tensor, &f1))
                        }
                    ).collect::<Vec<_>>()
                );

                a.push(
                    a[i - 1].iter().map(
                        |w| w[1..].to_vec()
                    ).collect::<Vec<_>>()
                );
            }
            // Step 2.c
            let ri = transcript.get_and_append_challenge(b"r")?;
            r.push(ri);
            // Step 2.d
            f.push(vector_add(
                &fe,
                &scalar_vector_product(ri, &fo)
            ));
            f_tilde.push(vector_add(
                &scalar_vector_product(F::ONE - ri, &f0),
                &scalar_vector_product(ri, &f1)
            ));
            // Step 2.e
            v.push(l[i].fft(&f[i]));
            if i == mu {
                f_mu = v[i][0];
            } else {
                let mti = MerkleTree::new(&v[i].iter().map(|&x| compute_sha256_row(&[x])).collect());
                mt_roots.push(mti.root().clone());
                mt.push(mti);    
            }
        }
        // Step 4
        for t in 0..s {
            // Step 4.a
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, l[0].size())?[0];
            // Step 4.b
            mt_proofs.push(Vec::new());
            for i in 0..mu {
                let offset = l[i + 1].size();
                let beta0 = if beta >= offset {beta - offset} else {beta + offset};
                mt_proofs[t].push((v[i][beta], mt[i].prove(beta), beta));
                mt_proofs[t].push((v[i][beta0], mt[i].prove(beta0), beta0));
                if beta >= offset {
                    beta -= offset;
                }
            }
        }
        Ok(DeepFoldProof{
            linear_polys,
            mt_roots,
            f_mu,
            mt_proofs,
        })
    }

    fn multi_open(
            prover_param: impl Borrow<Self::ProverParam>,
            polynomials: Vec<Self::Polynomial>,
            advices: &[Self::ProverCommitmentAdvice],
            points: &[Self::Point],
            _evals: &[Self::Evaluation],
            transcript: &mut IOPTranscript<F>,
        ) -> Result<Self::BatchProof, PCSError> {
        let &Self::ProverParam{max_mu, l0, s} = prover_param.borrow();
        let num_poly = polynomials.len();
        let mu = max_mu;
        let polynomials = polynomials.iter().map(|poly| resize_poly(&poly, mu)).collect::<Vec<_>>();
        let points = points.iter().map(|point| resize_point(&point, mu)).collect::<Vec<_>>();
        let mt0_list = advices.iter().map(|advice| &advice.mt0).collect::<Vec<_>>();

        // SumCheck Phase
        let start = std::time::Instant::now();
        let r = transcript.get_and_append_challenge(b"batched_sumcheck")?;
        let mut sum_check = VirtualPolynomial::new(max_mu);
        for i in 0..num_poly {
            sum_check.add_mle_list([
                evals_to_arcpoly(&polynomials[i].evaluations),
                evals_to_arcpoly(&get_tensor(&points[i])),
            ], r.pow([i as u64])).unwrap();
        }
        let sum_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(sum_check, transcript).unwrap();
        let point = sum_check_proof.point.clone();
        let sum_check_evals = polynomials.iter().map(
            |poly| eval_mle_poly(&poly.evaluations, &point)
        ).collect::<Vec<_>>();
        println!("DeepFoldPCS sumcheck : {} ms", start.elapsed().as_millis());

        // Batched Open Phase
        let start = std::time::Instant::now();
        let gamma = transcript.get_and_append_challenge_vectors(b"gamma", num_poly)?;
        let poly = evals_to_arcpoly(&(0..1 << max_mu).map(
            |i| (0..num_poly).map(
                |j| gamma[j] * polynomials[j].evaluations[i]
            ).sum::<F>()
        ).collect::<Vec<_>>());
        let f0 = evals_to_coeffs(mu, &poly.evaluations);
        let v0 = l0.fft(&f0);
        let mt0 = MerkleTree::new(&v0.iter().map(|&x| compute_sha256_row(&[x])).collect());
        let deepfold_prover_advice = DeepFoldProverCommitmentAdvice{f0, mt0, v0};
        println!("DeepFoldPCS before multi_open : {} ms", start.elapsed().as_millis());
        let start = std::time::Instant::now();
        let deepfold_proof = Self::open(prover_param, &poly, &deepfold_prover_advice, &point, transcript)?;
        println!("DeepFoldPCS multi_open : {} ms", start.elapsed().as_millis());

        // Additional checks for mt0
        let start = std::time::Instant::now();
        let mut mt_proofs_for_mt0 = Vec::new();
        for t in 0..s {
            mt_proofs_for_mt0.push(Vec::new());
            for k in 0..num_poly {
                mt_proofs_for_mt0[t].push((
                    advices[k].v0[deepfold_proof.mt_proofs[t][0].2], 
                    mt0_list[k].prove(deepfold_proof.mt_proofs[t][0].2), 
                    deepfold_proof.mt_proofs[t][0].2,
                ));
            }
        }

        let evals = polynomials.iter().zip(points.iter()).map(
            |(poly, point)| eval_mle_poly(&poly.evaluations, point)
        ).collect::<Vec<_>>();
        println!("DeepFoldPCS merkle tree : {} ms", start.elapsed().as_millis());
        

        Ok(Self::BatchProof{
            deepfold_proof,
            sum_check_proof,
            mt_proofs_for_mt0,
            evals,
            sum_check_evals,
        })
    }

    fn verify(
        verifier_param: &Self::VerifierParam,
        com: &Self::Commitment,
        point: &Self::Point,
        value: &F,
        proof: &Self::Proof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        let Self::VerifierParam{max_mu, len_l0, g, s} = verifier_param.clone();
        let Self::Commitment{mu, rt0} = com.clone();
        let mu = max_mu;
        let point = resize_point(&point, mu);
        assert!(mu <= max_mu);
        let Self::Proof{linear_polys, mt_roots, f_mu, mt_proofs} = proof.clone();
        
        if rt0 != mt_roots[0] {
            return Ok(false);
        }

        let mut alpha = vec![F::ZERO];
        let mut r = vec![F::ZERO];
        
        for _ in 1..mu + 1 {
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);
            r.push(transcript.get_and_append_challenge(b"r")?);
        }

        if eval_linear_poly(&linear_polys[0][0], &point[0]) != *value 
            || eval_linear_poly(&linear_polys[mu - 1][0], &r[mu]) != f_mu {
            return Ok(false);
        }

        for i in 1..mu {
            for j in 0..linear_polys[i - 1].len() {
                let k = if i < mu - 1 { j } else { 0 };
                let w1 = if j == 0 { point[i] } else {
                    alpha[j].pow([1 << (i + 1 - j) as u64])
                };
                if eval_linear_poly(&linear_polys[i - 1][j], &r[i])
                    != eval_linear_poly(&linear_polys[i][k], &w1) {
                    return Ok(false);
                }
            }
        }

        for t in 0..s {
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, len_l0)?[0];
            let mut beta_point = g.pow([beta as u64]);
            for i in 0..mu {
                let offset = len_l0 >> (i + 1);
                let beta0 = if beta >= offset {beta - offset} else {beta + offset};
                if !MerkleTree::verify(
                    &mt_roots[i],
                    beta,
                    &compute_sha256_row(&[mt_proofs[t][i * 2].0]),
                    &mt_proofs[t][i * 2].1,
                ) {
                    return Ok(false);
                }
                if !MerkleTree::verify(
                    &mt_roots[i],
                    beta0,
                    &compute_sha256_row(&[mt_proofs[t][i * 2 + 1].0]),
                    &mt_proofs[t][i * 2 + 1].1,
                ) {
                    return Ok(false);
                }

                let next_beta = if beta >= offset {beta - offset} else {beta};
                let val = if i < mu - 1 {mt_proofs[t][i * 2 + 2].0} else {f_mu};
                
                if !is_collinear( (beta_point, mt_proofs[t][i * 2].0), 
                                  (-beta_point, mt_proofs[t][i * 2 + 1].0),
                                  (r[i + 1], val) ) {
                    return Ok(false);
                }

                beta = next_beta;
                beta_point *= beta_point;
            }
        }

        Ok(true)
    }

    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        commitments: &[Self::Commitment],
        points: &[Self::Point],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        let Self::VerifierParam{max_mu, len_l0, g, s} = verifier_param.clone();
        let mu = max_mu;
        let num_poly = commitments.len();
        let points = points.iter().map(|point| resize_point(&point, mu)).collect::<Vec<_>>();
        assert!(points.len() == num_poly);
        let Self::BatchProof{deepfold_proof, sum_check_proof, mt_proofs_for_mt0, evals, sum_check_evals} = batch_proof.clone();

        // Sumcheck Phase
        let r = transcript.get_and_append_challenge(b"batched_sumcheck")?;
        let sum_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&sum_check_proof);
        if sum_check_sum != (0..num_poly).map(|k| r.pow([k as u64]) * evals[k]).sum::<F>() {
            return Ok(false);
        }
        let sum_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
            sum_check_sum, 
            &sum_check_proof, 
            &VPAuxInfo{max_degree: 2, num_variables: mu, phantom: PhantomData::<F>::default()}, 
            transcript).unwrap();
        let point = sum_check_proof.point.clone();
        if sum_check_claim.expected_evaluation != 
            (0..num_poly).map(|k| r.pow([k as u64]) * eval_mle_eq(&point, &points[k]) * sum_check_evals[k]).sum::<F>() {
            return Ok(false);
        }
        
        // Batched Open Phase
        let gamma = transcript.get_and_append_challenge_vectors(b"gamma", num_poly)?;
        let com = DeepFoldCommitment{
            mu,
            rt0: deepfold_proof.mt_roots[0].clone(),
        };
        let value = DeepFoldPCS::compute_value_from_proof(&point, &deepfold_proof);
        if value != (0..num_poly).map(|k| gamma[k] * sum_check_evals[k]).sum::<F>() {
            return Ok(false);
        }
        if !Self::verify(verifier_param, &com, &point, &value, &deepfold_proof, transcript)? {
            return Ok(false);
        }

        // Additional checks for mt0
        for t in 0..s {
            let mut sum = F::ZERO;
            for k in 0..num_poly {
                if !MerkleTree::verify(
                    &commitments[k].rt0,
                    deepfold_proof.mt_proofs[t][0].2,
                    &compute_sha256_row(&[mt_proofs_for_mt0[t][k].0]),
                    &mt_proofs_for_mt0[t][k].1,
                ) {
                    return Ok(false);
                }
                sum += gamma[k] * mt_proofs_for_mt0[t][k].0;
            }
            if sum != deepfold_proof.mt_proofs[t][0].0 {
                return Ok(false);
            }
        }

        Ok(true)
    }
}
