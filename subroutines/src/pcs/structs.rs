// Copyright (c) 2023 Espresso Systems (espressosys.com)
// This file is part of the HyperPlonk library.

// You should have received a copy of the MIT License
// along with the HyperPlonk library. If not, see <https://mit-license.org/>.

use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use derivative::Derivative;
use crate::{PolynomialCommitmentScheme, IOPProof};

// #[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
// #[derivative(
//     Default(bound = ""),
//     Hash(bound = ""),
//     Clone(bound = ""),
//     Copy(bound = ""),
//     Debug(bound = ""),
//     PartialEq(bound = ""),
//     Eq(bound = "")
// )]
// /// A commitment is an Affine point.
// pub struct Commitment<E: Pairing>(
//     /// the actual commitment is an affine point.
//     pub E::G1Affine,
// );

#[derive(Clone, Debug, Default, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct BatchProof<F, PCS>
where
    F: PrimeField,
    PCS: PolynomialCommitmentScheme<F>,
{
    /// A sum check proof proving tilde g's sum
    pub(crate) sum_check_proof: IOPProof<F>,
    /// f_i(point_i)
    pub f_i_eval_at_point_i: Vec<F>,
    /// proof for g'(a_2)
    pub(crate) g_prime_proof: PCS::Proof,
}
