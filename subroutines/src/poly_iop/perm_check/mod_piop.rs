// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

//! Main module for the Permutation Check protocol

use crate::{
    poly_iop::{errors::PolyIOPErrors, PolyIOP},
    MultiRationalSumcheck, MultiRationalSumcheckProof, PolynomialCommitmentScheme,
};
use arithmetic::{eq_eval, math::Math, products_except_self};
use ark_ec::pairing::Pairing;
use ark_ff::{batch_inversion, PrimeField};
use ark_poly::DenseMultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer, One, Zero};
use itertools::{izip, Itertools};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::sync::Arc;
use transcript::IOPTranscript;
use util::compute_leaves;

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

use super::multi_rational_sumcheck::MultiRationalSumcheckSubClaim;

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct PermutationCheckProof<F, PCS>
where
    F: PrimeField,
    PCS: PolynomialCommitmentScheme<F>,
{
    pub proofs: Vec<MultiRationalSumcheckProof<F>>,
    pub h_comms: Vec<PCS::Commitment>,
}

/// A permutation subclaim consists of
/// - the SubClaim from the ProductCheck
/// - Challenges beta and gamma
#[derive(Clone, Debug, Default, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct PermutationCheckSubClaim<F>
where
    F: PrimeField,
{
    pub subclaims: Vec<(MultiRationalSumcheckSubClaim<F>, usize)>,
    /// Challenges beta and gamma
    pub challenges: (F, F),
}

pub mod util;

/// A PermutationCheck w.r.t. `(fs, gs, perms)`
/// proves that (g1, ..., gk) is a permutation of (f1, ..., fk) under
/// permutation `(p1, ..., pk)`
/// It is derived from ProductCheck.
///
/// A Permutation Check IOP takes the following steps:
///
/// Inputs:
/// - fs = (f1, ..., fk)
/// - gs = (g1, ..., gk)
/// - permutation oracles = (p1, ..., pk)
pub trait PermutationCheck<F, PCS>
where
    F: PrimeField,
    PCS: PolynomialCommitmentScheme<F>,
{
    type PermutationCheckSubClaim;
    type PermutationProof: CanonicalSerialize + CanonicalDeserialize;

    type MultilinearExtension;
    type Transcript;

    /// Initialize the system with a transcript
    ///
    /// This function is optional -- in the case where a PermutationCheck is
    /// an building block for a more complex protocol, the transcript
    /// may be initialized by this complex protocol, and passed to the
    /// PermutationCheck prover/verifier.
    fn init_transcript() -> Self::Transcript;

    /// Inputs:
    /// - fs = (f1, ..., fk)
    /// - gs = (g1, ..., gk)
    /// - permutation oracles = (p1, ..., pk)
    /// Outputs:
    /// - a permutation check proof proving that gs is a permutation of fs under
    ///   permutation
    ///
    /// Cost: O(N)
    #[allow(clippy::type_complexity)]
    fn prove(
        prover_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gxs: &[Self::MultilinearExtension],
        perms: &[Self::MultilinearExtension],
        transcript: &mut IOPTranscript<F>,
    ) -> Result<
        (
            Self::PermutationProof,
            Vec<PCS::ProverCommitmentAdvice>,
            Vec<Vec<F>>,
            Vec<Arc<DenseMultilinearExtension<F>>>,
        ),
        PolyIOPErrors,
    >;

    fn d_prove_prepare(
        prover_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gxs: &[Self::MultilinearExtension],
        perms: &[Self::MultilinearExtension],
        transcript: &mut IOPTranscript<F>,
    ) -> Result<
        (
            Vec<(
                Vec<Arc<DenseMultilinearExtension<F>>>,
                Arc<DenseMultilinearExtension<F>>,
                F,
                Option<PCS::Commitment>,
                PCS::ProverCommitmentAdvice,
            )>,
            Vec<F>,
        ),
        PolyIOPErrors,
    >;

    fn d_prove(
        prover_param: &PCS::ProverParam,
        to_prove: Vec<(
            Vec<Arc<DenseMultilinearExtension<F>>>,
            Arc<DenseMultilinearExtension<F>>,
            F,
            Option<PCS::Commitment>,
            PCS::ProverCommitmentAdvice,
        )>,
        claims: Vec<F>,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<
        (
            Option<(Self::PermutationProof, Vec<Vec<F>>)>,
            Vec<PCS::ProverCommitmentAdvice>,
            Vec<Arc<DenseMultilinearExtension<F>>>,
        ),
        PolyIOPErrors,
    >;

    /// Verify that (g1, ..., gk) is a permutation of
    /// (f1, ..., fk) over the permutation oracles (perm1, ..., permk)
    fn verify(
        proof: &Self::PermutationProof,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors>;

    fn check_openings(
        subclaim: &Self::PermutationCheckSubClaim,
        f_openings: &[F],
        g_openings: &[F],
        h_openings: &[F],
        perm_openings: &[F],
    ) -> Result<(), PolyIOPErrors>;
}

impl<F, PCS> PermutationCheck<F, PCS> for PolyIOP<F>
where
    F: PrimeField,
    PCS: PolynomialCommitmentScheme<F, Polynomial = Arc<DenseMultilinearExtension<F>>>,
{
    type PermutationCheckSubClaim = PermutationCheckSubClaim<F>;
    type PermutationProof = PermutationCheckProof<F, PCS>;
    type MultilinearExtension = Arc<DenseMultilinearExtension<F>>;
    type Transcript = IOPTranscript<F>;

    fn init_transcript() -> Self::Transcript {
        IOPTranscript::<F>::new(b"Initializing PermutationCheck transcript")
    }

    // Strictly speaking the list of points is redundant as it is present in the
    // proofs, but we try to keep the interface uniform
    fn prove(
        prover_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gxs: &[Self::MultilinearExtension],
        perms: &[Self::MultilinearExtension],
        transcript: &mut IOPTranscript<F>,
    ) -> Result<
        (
            Self::PermutationProof,
            Vec<PCS::ProverCommitmentAdvice>,
            Vec<Vec<F>>,
            Vec<Arc<DenseMultilinearExtension<F>>>,
        ),
        PolyIOPErrors,
    > {
        let start = start_timer!(|| "Permutation check prove");
        if fxs.is_empty() {
            return Err(PolyIOPErrors::InvalidParameters("fxs is empty".to_string()));
        }
        if (fxs.len() != gxs.len()) || (fxs.len() != perms.len()) {
            return Err(PolyIOPErrors::InvalidProof(format!(
                "fxs.len() = {}, gxs.len() = {}, perms.len() = {}",
                fxs.len(),
                gxs.len(),
                perms.len(),
            )));
        }

        // generate challenge `beta` and `gamma` from current transcript
        let beta = transcript.get_and_append_challenge(b"beta")?;
        let gamma = transcript.get_and_append_challenge(b"gamma")?;
        let leaves = compute_leaves::<F, false>(&beta, &gamma, fxs, gxs, perms)?;

        let leaves_len = leaves.len();

        let to_prove = leaves
            .into_iter()
            // .into_par_iter()
            .map(|leave| {
                let half_len = leave.len() / 2;
                let nv = leave[0].len().log_2();
                let (g_polys, inv_evals): (Vec<_>, Vec<_>) = leave
                    .into_par_iter()
                    .map(|evals| {
                        let mut inv_evals = evals.clone();
                        batch_inversion(&mut inv_evals);

                        (
                            Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, evals)),
                            inv_evals,
                        )
                    })
                    .unzip();
                let h_evals = (0..inv_evals[0].len())
                    .into_par_iter()
                    .map(|i| {
                        inv_evals[..half_len]
                            .iter()
                            .map(|eval| eval[i])
                            .sum::<F>()
                            - inv_evals[half_len..]
                                .iter()
                                .map(|eval| eval[i])
                                .sum::<F>()
                    })
                    .collect::<Vec<_>>();
                let claim = if leaves_len == 1 {
                    F::zero()
                } else {
                    h_evals.iter().sum::<F>()
                };

                let h_poly = Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, h_evals));
                let (h_comm, h_advice) = PCS::commit(prover_param, &h_poly).unwrap();

                (g_polys, h_poly, claim, h_comm, h_advice)
            })
            .collect::<Vec<_>>();

        let (proofs, points, comms, advices, polys): (Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
            to_prove
                .into_iter()
                .map(|(g_polys, h_poly, claim, h_comm, h_advice)| {
                    let mut f_values = vec![F::one(); g_polys.len()];
                    f_values[g_polys.len() / 2..].fill(-F::one());
                    let (proof, point) = <Self as MultiRationalSumcheck<F>>::prove(
                        &f_values,
                        g_polys,
                        Arc::new(DenseMultilinearExtension::clone(&h_poly)),
                        claim,
                        transcript,
                    )
                    .unwrap();
                    (proof, point, h_comm, h_advice, h_poly)
                })
                .multiunzip();

        end_timer!(start);

        Ok((
            Self::PermutationProof {
                proofs,
                h_comms: comms,
            },
            advices,
            points,
            polys,
        ))
    }

    fn d_prove_prepare(
        prover_param: &PCS::ProverParam,
        fxs: &[Self::MultilinearExtension],
        gxs: &[Self::MultilinearExtension],
        perms: &[Self::MultilinearExtension],
        transcript: &mut IOPTranscript<F>,
    ) -> Result<
        (
            Vec<(
                Vec<Arc<DenseMultilinearExtension<F>>>,
                Arc<DenseMultilinearExtension<F>>,
                F,
                Option<PCS::Commitment>,
                PCS::ProverCommitmentAdvice,
            )>,
            Vec<F>,
        ),
        PolyIOPErrors,
    > {
        let start = start_timer!(|| "Permutation check prove");
        if fxs.is_empty() {
            return Err(PolyIOPErrors::InvalidParameters("fxs is empty".to_string()));
        }
        if (fxs.len() != gxs.len()) || (fxs.len() != perms.len()) {
            return Err(PolyIOPErrors::InvalidProof(format!(
                "fxs.len() = {}, gxs.len() = {}, perms.len() = {}",
                fxs.len(),
                gxs.len(),
                perms.len(),
            )));
        }

        let (beta, gamma) = if Net::am_master() {
            let beta = transcript.get_and_append_challenge(b"beta")?;
            let gamma = transcript.get_and_append_challenge(b"gamma")?;
            Net::recv_from_master_uniform(Some((beta, gamma)))
        } else {
            Net::recv_from_master_uniform(None)
        };

        let leaves = compute_leaves::<F, true>(&beta, &gamma, fxs, gxs, perms)?;

        let leaves_len = leaves.len();
        let to_prove = leaves
            .into_iter()
            .map(|leave| {
                let half_len = leave.len() / 2;
                let nv = leave[0].len().log_2();
                let (g_polys, inv_evals): (Vec<_>, Vec<_>) = leave
                    .into_par_iter()
                    .map(|evals| {
                        let mut inv_evals = evals.clone();
                        batch_inversion(&mut inv_evals);

                        (
                            Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, evals)),
                            inv_evals,
                        )
                    })
                    .unzip();
                let h_evals = (0..inv_evals[0].len())
                    .into_par_iter()
                    .map(|i| {
                        inv_evals[..half_len].iter().map(|eval| eval[i]).sum::<F>()
                            - inv_evals[half_len..].iter().map(|eval| eval[i]).sum::<F>()
                    })
                    .collect::<Vec<_>>();
                let claim = if leaves_len == 1 {
                    F::zero()
                } else {
                    h_evals.iter().sum::<F>()
                };
                let h_poly = Arc::new(DenseMultilinearExtension::from_evaluations_vec(nv, h_evals));
                let (h_comm, h_advice) = PCS::d_commit(prover_param, &h_poly, transcript).unwrap();

                (g_polys, h_poly, claim, h_comm, h_advice)
            })
            .collect::<Vec<_>>();

        let mut claims = to_prove
            .iter()
            .map(|(_, _, claim, ..)| *claim)
            .collect::<Vec<_>>();
        let all_claims = Net::send_to_master(&claims);
        if Net::am_master() {
            let all_claims = all_claims.unwrap();
            claims = (0..all_claims[0].len())
                .map(|i| all_claims.iter().map(|claims| claims[i]).sum::<F>())
                .collect::<Vec<_>>();
        }

        end_timer!(start);

        Ok((to_prove, claims))
    }

    fn d_prove(
        _prover_param: &PCS::ProverParam,
        to_prove: Vec<(
            Vec<Arc<DenseMultilinearExtension<F>>>,
            Arc<DenseMultilinearExtension<F>>,
            F,
            Option<PCS::Commitment>,
            PCS::ProverCommitmentAdvice,
        )>,
        claims: Vec<F>,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<
        (
            Option<(Self::PermutationProof, Vec<Vec<F>>)>,
            Vec<PCS::ProverCommitmentAdvice>,
            Vec<Arc<DenseMultilinearExtension<F>>>,
        ),
        PolyIOPErrors,
    > {
        let start = start_timer!(|| "Permutation check prove");

        if !Net::am_master() {
            let (advices, polys): (Vec<_>, Vec<_>) = to_prove
                .into_iter()
                .map(|(g_polys, h_poly, _, _, h_advice)| {
                    let mut f_values = vec![F::one(); g_polys.len()];
                    f_values[g_polys.len() / 2..].fill(-F::one());

                    <Self as MultiRationalSumcheck<F>>::d_prove(
                        &f_values,
                        g_polys,
                        Arc::new(DenseMultilinearExtension::clone(&h_poly)),
                        F::zero(),
                        transcript,
                    )
                    .unwrap();
                    (h_advice, h_poly)
                })
                .unzip();

            end_timer!(start);
            return Ok((None, advices, polys));
        }

        let (proofs, points, comms, advices, polys): (Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
            to_prove
                .into_iter()
                .zip(claims)
                .map(|((g_polys, h_poly, _, h_comm, h_advice), claim)| {
                    let mut f_values = vec![F::one(); g_polys.len()];
                    f_values[g_polys.len() / 2..].fill(-F::one());

                    let (proof, point) = <Self as MultiRationalSumcheck<F>>::d_prove(
                        &f_values,
                        g_polys,
                        Arc::new(DenseMultilinearExtension::clone(&h_poly)),
                        claim,
                        transcript,
                    )
                    .unwrap()
                    .unwrap();
                    (proof, point, h_comm, h_advice, h_poly)
                })
                .multiunzip();

        end_timer!(start);

        let comms = comms
            .into_iter()
            .map(|comm| comm.unwrap())
            .collect::<Vec<_>>();
        Ok((
            Some((
                Self::PermutationProof {
                    proofs,
                    h_comms: comms,
                },
                points,
            )),
            advices,
            polys,
        ))
    }

    fn verify(
        proof: &Self::PermutationProof,
        transcript: &mut Self::Transcript,
    ) -> Result<Self::PermutationCheckSubClaim, PolyIOPErrors> {
        let start = start_timer!(|| "Permutation check verify");

        let beta = transcript.get_and_append_challenge(b"beta")?;
        let gamma = transcript.get_and_append_challenge(b"gamma")?;

        let mut subclaims = Vec::with_capacity(proof.proofs.len());
        let mut claimed_sum = F::zero();
        for proof in proof.proofs.iter() {
            claimed_sum += proof.claimed_sum;
            let subclaim = <Self as MultiRationalSumcheck<F>>::verify(proof, transcript)?;
            subclaims.push((subclaim, proof.num_polys / 2));
        }

        if claimed_sum != F::zero() {
            return Err(PolyIOPErrors::InvalidProof(format!(
                "Claimed sums do not add to zero",
            )));
        }

        end_timer!(start);
        Ok(PermutationCheckSubClaim {
            subclaims,
            challenges: (beta, gamma),
        })
    }

    fn check_openings(
        subclaim: &Self::PermutationCheckSubClaim,
        f_openings: &[F],
        g_openings: &[F],
        h_openings: &[F],
        perm_openings: &[F],
    ) -> Result<(), PolyIOPErrors> {
        let (beta, gamma) = subclaim.challenges;

        let mut shift = 0;
        let mut offset = 0;
        for (subclaim_idx, (subclaim, len)) in subclaim.subclaims.iter().enumerate() {
            let num_vars = subclaim.sumcheck_point.len();

            let sid: F = (0..num_vars)
                .map(|i| F::from_u64(i.pow2() as u64).unwrap() * subclaim.sumcheck_point[i])
                .sum::<F>()
                + F::from_u64(shift as u64).unwrap();

            let eq_eval = eq_eval(&subclaim.sumcheck_point, &subclaim.zerocheck_r)?;
            let g_evals = f_openings[offset..offset + len]
                .iter()
                .enumerate()
                .map(|(i, f)| *f + beta * (sid + F::from((i * (1 << num_vars)) as u64)) + gamma)
                .chain(
                    g_openings[offset..offset + len]
                        .iter()
                        .zip(perm_openings[offset..offset + len].iter())
                        .map(|(g, perm)| *g + beta * perm + gamma),
                )
                .collect::<Vec<_>>();
            let g_products = products_except_self(&g_evals);
            let sum = h_openings[subclaim_idx]
                + subclaim.coeff
                    * eq_eval
                    * (g_products[0] * g_evals[0] * h_openings[subclaim_idx]
                        - g_products[..*len].iter().sum::<F>()
                        + g_products[*len..].iter().sum::<F>());

            if sum != subclaim.sumcheck_expected_evaluation {
                return Err(PolyIOPErrors::InvalidVerifier("wrong subclaim".to_string()));
            }

            shift += len * num_vars.pow2();
            offset += len;
        }

        Ok(())
    }
}
