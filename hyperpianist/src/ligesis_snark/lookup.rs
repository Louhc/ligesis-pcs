// use super::*;

// use crate::{
//     errors::HyperPlonkErrors,
//     lookup::HyperPlonkLookupPlugin,
//     prelude::HyperPlonkParams,
//     structs::{HyperPlonkIndex, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey},
//     utils::{
//         build_f, eval_f, prover_sanity_check, PcsDynamicAccumulator, PcsDynamicOpenings,
//         PcsDynamicVerifier,
//     },
//     witness::WitnessColumn,
//     HyperPlonkSNARK,
// };
// use arithmetic::{evaluate_opt, math::Math, VPAuxInfo};
// use ark_ec::pairing::Pairing;
// use ark_ff::PrimeField;
// use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
// use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
// use ark_std::{end_timer, log2, start_timer, Zero};
// use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};
// use itertools::izip;
// use lazy_static::lazy_static;
// #[cfg(feature = "parallel")]
// use rayon::iter::ParallelIterator;
// use rayon::{iter::IntoParallelRefIterator, ThreadPoolBuilder};
// use std::{iter::zip, marker::PhantomData, mem::take, sync::Arc};
// use subroutines::{
//     pcs::prelude::HashBasedPCS,
//     poly_iop::{prelude::CombinedCheck, PolyIOP},
//     BatchProof,
// };
// use transcript::IOPTranscript;

// pub trait HashBasedPlonkLookupPlugin<F, PCS>
// where
//     F: PrimeField,
//     PCS: HashBasedPCS<F>,
// {
//     type Ops: Sync;
//     type Preprocessing: Clone + Sync;
//     type Transcript;
//     type Proof: Send + Sync + CanonicalSerialize + CanonicalDeserialize;

//     fn preprocess() -> Self::Preprocessing;
//     fn construct_witnesses(ops: &Self::Ops) -> Vec<Arc<DenseMultilinearExtension<E::ScalarField>>>;
//     fn num_witness_columns() -> Vec<usize>;
//     fn max_num_variables() -> usize;
//     fn prove(
//         preprocessing: &Self::Preprocessing,
//         pcs_param: &PCS::ProverParam,
//         ops: &Self::Ops,
//         transcript: &mut Self::Transcript,
//     ) -> (Self::Proof, HyperPlonkLookupProverOpeningPoints<E, PCS>);
//     fn d_prove(
//         preprocessing: &Self::Preprocessing,
//         pcs_param: &PCS::ProverParam,
//         ops: &Self::Ops,
//         transcript: &mut Self::Transcript,
//     ) -> (
//         Option<Self::Proof>,
//         HyperPlonkLookupProverOpeningPoints<E, PCS>,
//     );
//     fn num_regular_openings(proof: &Self::Proof) -> usize;
//     fn verify(
//         proof: &Self::Proof,
//         witness_openings: &[E::ScalarField],
//         regular_openings: &[E::ScalarField],
//         transcript: &mut Self::Transcript,
//     ) -> Result<HyperPlonkLookupVerifierOpeningPoints<E, PCS>, PolyIOPErrors>;

// }