use super::*;

use crate::{
    errors::HyperPlonkErrors,
    lookup::HyperPlonkLookupPlugin,
    prelude::HyperPlonkParams,
    structs::{HyperPlonkIndex, HyperPlonkProof, HyperPlonkProvingKey, HyperPlonkVerifyingKey},
    utils::{
        build_f, eval_f, prover_sanity_check, PcsDynamicAccumulator, PcsDynamicOpenings,
        PcsDynamicVerifier,
    },
    witness::WitnessColumn,
    HyperPlonkSNARK,
};
use arithmetic::{evaluate_opt, math::Math, VPAuxInfo};
use ark_ec::pairing::Pairing;
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_std::{end_timer, log2, start_timer, Zero};
use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};
use itertools::izip;
use lazy_static::lazy_static;
#[cfg(feature = "parallel")]
use rayon::iter::ParallelIterator;
use rayon::{iter::IntoParallelRefIterator, ThreadPoolBuilder};
use std::{iter::zip, marker::PhantomData, mem::take, sync::Arc};
use subroutines::{
    pcs::prelude::HashBasedPCS,
    poly_iop::{prelude::CombinedCheck, PolyIOP},
    BatchProof,
};
use transcript::IOPTranscript;

pub mod structs;
pub use structs::*;
pub mod lookup;
pub use lookup::*;

impl<F, PCS, Lookup> HashBasedHyperPlonkSNARK<F, PCS, Lookup> for PolyIOP<F>
where
    F: PrimeField,
    PCS: HashBasedPCS<
        F,
        Polynomial = Arc<DenseMultilinearExtension<F>>,
        Point = Vec<F>,
        Evaluation = F,
        BatchProof = (),
    >,
{
    type Index = HyperPlonkIndex<F>;
    type ProvingKey = ();
    type VerifyingKey = ();
    type Proof = ();

    fn preprocess(
        index: &Self::Index,
        pcs_srs: &PCS::SRS,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), HyperPlonkErrors> {
        
        Ok((
            (),
            (),
        ))
    }

    fn d_preprocess(
        index: &Self::Index,
        pcs_srs: &PCS::SRS,
    ) -> Result<(Self::ProvingKey, Option<Self::VerifyingKey>), HyperPlonkErrors> {
        unimplemented!()
    }

    fn prove(
        pk: &Self::ProvingKey,
        pub_input: &[F],
        witnesses: &[WitnessColumn<F>],
        // ops: &Lookup::Ops,
    ) -> Result<Self::Proof, HyperPlonkErrors> {
        unimplemented!()
    }

    fn d_prove(
        pk: &Self::ProvingKey,
        pub_input: &[F],
        witnesses: &[WitnessColumn<F>],
        // ops: &Lookup::Ops,
    ) -> Result<Option<Self::Proof>, HyperPlonkErrors> {
        unimplemented!()
    }

    fn verify(
        vk: &Self::VerifyingKey,
        pub_input: &[F],
        proof: &Self::Proof,
    ) -> Result<bool, HyperPlonkErrors> {
        unimplemented!()
    }
}
