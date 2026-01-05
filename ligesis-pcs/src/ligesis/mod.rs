mod commit;
mod open;
mod verify;

pub use commit::{ligesis_commit, ligesis_d_commit};
pub use open::{ligesis_open, ligesis_d_open};
pub use verify::ligesis_verify;

use crate::{
    deepfold::*, errors::PCSError, rand::*, rscode::*, utils::*,
    PolynomialCommitmentScheme,
    ext_sumcheck::ExtSumCheckProof,
    types::{HasQuadraticExtension, FieldExtension},
};
use ark_ff::{BigInteger, PrimeField};
use ark_poly::DenseMultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow,
    cmp::{max, min},
    marker::PhantomData,
    rand::Rng,
    sync::Arc,
    vec::Vec,
};
use transcript::IOPTranscript;

pub use crate::FGoldilocks;

/// Compute the SIS hash matrix H = A' * B where B is the byte decomposition of F'.
///
/// Parameters:
/// - `mat_a`: The A matrix of shape `c x (eta * m_rows)`
/// - `mat_f_prime`: The RS-encoded F' matrix of shape `m_rows x cols`
/// - `eta`: The eta parameter (bit length)
/// - `m_rows`: Number of rows in mat_f_prime (can be m or m/num_party for distributed)
///
/// Returns `mat_h` of shape `c x cols`
fn compute_sis_hash<F: PrimeField>(
    mat_a: &[Vec<F>],
    mat_f_prime: &[Vec<F>],
    eta: usize,
    m_rows: usize,
) -> Vec<Vec<F>> {
    let c = mat_a.len();
    let cols = mat_f_prime[0].len();

    // Precompute byte buckets for mat_a
    let mat_a_prime: Vec<Vec<Vec<F>>> = mat_a
        .iter()
        .map(|row| {
            (0..eta * m_rows / 8)
                .map(|i| get_mat_a_byte_bucket(&row[i * 8..(i + 1) * 8].to_vec()))
                .collect()
        })
        .collect();

    // Compute H using SIS hash
    mat_a_prime
        .iter()
        .map(|row| {
            (0..cols)
                .map(|j| {
                    (0..m_rows)
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
                .collect()
        })
        .collect()
}

// TODO: Lookup code is currently incorrect, commented out for now
// #[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
// /// proof of lookup
// pub struct LigeSISLookupProof<F: PrimeField> {
//     pub sumcheck_proof: IOPProof<F>,
// }

#[cfg(test)]
mod tests;

/// LigeSIS Polynomial Commitment Scheme
pub struct LigeSISPCS<F: PrimeField> {
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField + HasQuadraticExtension> LigeSISPCS<F> {
    pub fn compute_value_from_proof(_log_n: usize, _point: &Vec<F>, proof: &LigeSISProof<F>) -> F {
        F::ext_real(&proof.deepfold_batched_proof.evals[0])
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, Default)]
pub struct LigeSISSRS<F: PrimeField> {
    lambda: usize,
    eta: usize,
    pub mu: usize,
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
    c: usize,
    rs: ReedSolomon<F>,
    mat_a: Vec<Vec<F>>,
    mat_a_pad: Arc<DenseMultilinearExtension<F>>,
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
    com_mat_a: DeepFoldCommitment,
    deepfold_verifier_param: DeepFoldVerifierParam<F>,
}

/// Extension field SumCheck proof (direct, no reduction needed)
/// Used with direct extension field opening in DeepFold (128-bit soundness)
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct ExtSumCheckWithReductionProof<F: PrimeField + HasQuadraticExtension> {
    /// Extension field SumCheck proof (contains the extension field point)
    pub ext_proof: ExtSumCheckProof<F::Extension>,
}

/// Proof for Ligesis - provides 128-bit soundness
/// Uses extension field SumCheck and direct extension field opening in DeepFold
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct LigeSISProof<F: PrimeField + HasQuadraticExtension> {
    pub com_a: <DeepFoldPCS<F> as PolynomialCommitmentScheme<F>>::Commitment,
    pub com_bI: <DeepFoldPCS<F> as PolynomialCommitmentScheme<F>>::Commitment,
    pub com_rs_a: <DeepFoldPCS<F> as PolynomialCommitmentScheme<F>>::Commitment,

    /// Extension field SumCheck proofs (no reduction needed)
    pub bI_check_proof: ExtSumCheckWithReductionProof<F>,
    pub alpha2_a_bI_r2_check_proof: ExtSumCheckWithReductionProof<F>,
    pub v_bI_r2_check_proof: ExtSumCheckWithReductionProof<F>,
    pub rs_a_check_proof: ExtSumCheckWithReductionProof<F>,
    pub mat_g_check_proofs: Vec<ExtSumCheckWithReductionProof<F>>,

    /// Extension field DeepFold batch proof (direct extension field opening)
    pub deepfold_batched_proof: DeepFoldExtBatchedProof<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Debug, PartialEq, Eq, Default)]
pub struct LigeSISProverCommitmentAdvice<F: PrimeField> {
    pub mat_f_prime: Vec<Vec<F>>,
    pub mat_h: Vec<Vec<F>>,
    pub mat_h_pad: Arc<DenseMultilinearExtension<F>>,
    pub com_mat_h_advice: <DeepFoldPCS<F> as PolynomialCommitmentScheme<F>>::ProverCommitmentAdvice,
}

impl<F: PrimeField> Clone for LigeSISProverCommitmentAdvice<F> {
    fn clone(&self) -> Self {
        LigeSISProverCommitmentAdvice {
            mat_f_prime: self.mat_f_prime.clone(),
            mat_h: self.mat_h.clone(),
            mat_h_pad: Arc::clone(&self.mat_h_pad),
            com_mat_h_advice: self.com_mat_h_advice.clone(),
        }
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct LigeSISCommitment<F: PrimeField> {
    pub num_vars: usize,
    pub com_mat_h: <DeepFoldPCS<F> as PolynomialCommitmentScheme<F>>::Commitment,
}

impl<F: PrimeField + HasQuadraticExtension> PolynomialCommitmentScheme<F> for LigeSISPCS<F> {
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

    fn gen_srs_for_testing<R: Rng>(rng: &mut R, log_size: usize) -> Result<Self::SRS, PCSError> {
        let eta = F::ONE.into_bigint().to_bits_be().len();
        let lambda = 128usize;
        let mu = log_size;
        let log_m = if log_size < 4 { 0 } else { (log_size - 8) / 2 };
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
        let (deepfold_prover_param, deepfold_verifier_param) = DeepFoldPCS::<F>::setup(
            deepfold_srs,
        )?;

        let mat_a_pad = evals_to_arcpoly(&resize_eval(&mat_a.concat(), deepfold_srs.max_mu));
        let (com_mat_a, com_mat_a_advice) =
            DeepFoldPCS::commit(&deepfold_prover_param, &mat_a_pad)?;

        let rs = ReedSolomon::<F>::new(n, rs_len);
        let g = rs.get_generator();

        let prover_param = LigeSISProverParam {
            eta,
            s_lambda,
            mu,
            log_m,
            log_n,
            c,
            rs,
            mat_a,
            mat_a_pad,
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
            deepfold_verifier_param,
        };
        Ok((prover_param, verifier_param))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        ligesis_commit(prover_param.borrow(), poly)
    }

    fn d_commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<(Option<Self::Commitment>, Self::ProverCommitmentAdvice), PCSError> {
        ligesis_d_commit(prover_param.borrow(), poly)
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        ligesis_open(prover_param.borrow(), poly, advice, point, transcript)
    }

    fn d_open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Option<Self::Proof>, PCSError> {
        ligesis_d_open(prover_param.borrow(), poly, advice, point, transcript)
    }

    fn verify(
        verifier_param: &Self::VerifierParam,
        com: &Self::Commitment,
        point: &Self::Point,
        value: &F,
        proof: &Self::Proof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        ligesis_verify(verifier_param, com, point, value, proof, transcript)
    }
}

