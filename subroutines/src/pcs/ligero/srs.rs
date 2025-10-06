use crate::pcs::{
    hashpcs::HashBasedSRS, multilinear_kzg::util::{eq_eval, eq_extension}, prelude::PCSError, StructuredReferenceString
};
use ark_ec::{pairing::Pairing, scalar_mul::fixed_base::FixedBase, AffineRepr, CurveGroup};
use ark_ff::{Field, PrimeField, Zero};
use ark_poly::DenseMultilinearExtension;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    collections::LinkedList, end_timer, format, rand::Rng, start_timer, string::ToString, vec::Vec,
    UniformRand,
};
use core::iter::FromIterator;
use std::marker::PhantomData;

/// Evaluations over {0,1}^n for G1 or G2
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct Evaluations<C: AffineRepr> {
    /// The evaluations.
    pub evals: Vec<C>,
}

/// Universal Parameter
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct LigeroUniversalParams<F: PrimeField> {
    pub num_vars: usize,
    pub log_n: usize,
    pub code_len: usize,

    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField> std::default::Default for LigeroUniversalParams<F> {
    fn default() -> Self {
        Self {
            num_vars: 0,
            log_n: 0,
            code_len: 0,
            phantom: PhantomData::default(),
        }
    }
}


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct LigeroProverParam<F: PrimeField> {
    pub num_vars: usize,
    pub log_n: usize,
    pub code_len: usize,

    #[doc(hidden)]
    phantom: PhantomData<F>,
}


#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug)]
pub struct LigeroVerifierParam<F: PrimeField> {
    pub num_vars: usize,
    pub log_n: usize,
    pub code_len: usize,

    #[doc(hidden)]
    phantom: PhantomData<F>,
}

impl<F: PrimeField> HashBasedSRS<F> for LigeroUniversalParams<F> {
    type ProverParam = LigeroProverParam<F>;
    type VerifierParam = LigeroVerifierParam<F>;

    /// Extract the prover parameters from the public parameters.
    fn extract_prover_param(&self) -> Self::ProverParam {
        LigeroProverParam{
            num_vars: self.num_vars,
            log_n: self.log_n,
            code_len: self.code_len,
            phantom: self.phantom,
        }
    }

    /// Extract the verifier parameters from the public parameters.
    fn extract_verifier_param(&self) -> Self::VerifierParam {
        LigeroVerifierParam{
            num_vars: self.num_vars,
            log_n: self.log_n,
            code_len: self.code_len,
            phantom: self.phantom,
        }
    }

    /// Trim the universal parameters to specialize the public parameters
    /// for multilinear polynomials.
    fn trim(&self) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        Ok((self.extract_prover_param(), self.extract_verifier_param()))
    }

    /// Build SRS for testing.
    /// WARNING: THIS FUNCTION IS FOR TESTING PURPOSE ONLY.
    /// THE OUTPUT SRS SHOULD NOT BE USED IN PRODUCTION.
    fn gen_srs_for_testing(rng: &mut impl Rng, num_vars: usize) -> Result<Self, PCSError> {
        let log_n = num_vars / 2;
        let code_len = (1 << (num_vars - log_n)) * rng.gen_range(2..4);
        Ok(LigeroUniversalParams { num_vars, log_n, code_len, phantom: PhantomData::default() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::test_rng;
    type F = ark_bls12_381::Fr;

    #[test]
    fn test_srs_gen() -> Result<(), PCSError> {
        let mut rng = test_rng();
        for nv in 4..10 {
            let _ = LigeroUniversalParams::<F>::gen_srs_for_testing(&mut rng, nv)?;
        }

        Ok(())
    }
}
