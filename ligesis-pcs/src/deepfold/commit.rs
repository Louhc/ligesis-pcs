//! DeepFold commit functions
//!
//! This module contains the commit implementations for DeepFold PCS:
//! - `deepfold_commit`: Standard commit
//! - `deepfold_d_commit`: Distributed commit

use crate::{errors::PCSError, hash::*, utils::*};
use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::{end_timer, start_timer, sync::Arc, vec::Vec};
use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

use super::{
    DeepFoldCommitment, DeepFoldProverCommitmentAdvice, DeepFoldProverParam,
    utils::{build_merkle_tree, compute_leaf_hashes},
};

/// Standard DeepFold commit
pub fn deepfold_commit<F: PrimeField>(
    prover_param: &DeepFoldProverParam<F>,
    poly: &Arc<ark_poly::DenseMultilinearExtension<F>>,
) -> Result<(DeepFoldCommitment, DeepFoldProverCommitmentAdvice<F>), PCSError> {
    let DeepFoldProverParam { max_mu, l0, s: _ } = prover_param;
    let mu = poly.num_vars;
    assert!(mu <= *max_mu);

    let f0 = evals_to_coeffs(mu, &poly.evaluations);
    let v0 = l0.fft(&f0);

    let mt0 = build_merkle_tree(&v0);

    let rt0 = mt0.root();
    Ok((
        DeepFoldCommitment { mu, rt0 },
        DeepFoldProverCommitmentAdvice {
            f0,
            mt0,
            v0,
            f_tilde: poly.evaluations.clone(),
            upper_tree: None,
        },
    ))
}

/// Distributed commit: each party has local polynomial evaluations
/// Each party builds local subtree, master builds upper tree from collected roots
/// Returns (Option<Commitment>, Advice) - commitment is Some only for master
pub fn deepfold_d_commit<F: PrimeField>(
    prover_param: &DeepFoldProverParam<F>,
    poly: &Arc<ark_poly::DenseMultilinearExtension<F>>,
) -> Result<(Option<DeepFoldCommitment>, DeepFoldProverCommitmentAdvice<F>), PCSError> {
    let DeepFoldProverParam { max_mu, l0, s: _ } = prover_param;
    let num_party = Net::n_parties();
    let num_party_vars = num_party.ilog2() as usize;

    // Each party has local evaluations of size 2^local_mu
    let local_mu = poly.num_vars;
    let mu = local_mu + num_party_vars;
    assert!(mu <= *max_mu);

    // Step 1: Gather all evaluations to master
    let timer = start_timer!(|| "DCommit.GatherEvals");
    let all_evals_opt = Net::send_to_master(&poly.evaluations);
    end_timer!(timer);

    // Step 2: Master computes full f0, v0, and distributes leaf hashes
    let (f0, v0, f_tilde, local_leaves, leaf_size): (Vec<F>, Vec<F>, Vec<F>, Vec<Byte32>, usize) =
        if Net::am_master() {
            let all_evals: Vec<Vec<F>> = all_evals_opt.unwrap();
            let full_evals: Vec<F> = all_evals.into_iter().flatten().collect();

            // Compute full coefficients and FFT
            let timer = start_timer!(|| "DCommit.FFT");
            let f0 = evals_to_coeffs(mu, &full_evals);
            let v0 = l0.fft(&f0);
            end_timer!(timer);

            // Compute leaf hashes from v0
            let timer = start_timer!(|| "DCommit.LeafHashes");
            let (all_leaves, leaf_size) = compute_leaf_hashes(&v0);
            end_timer!(timer);

            // Split leaf hashes into chunks for each party
            let chunk_size = all_leaves.len() / num_party;
            let leaf_chunks: Vec<Vec<Byte32>> = (0..num_party)
                .map(|i| all_leaves[i * chunk_size..(i + 1) * chunk_size].to_vec())
                .collect();

            // Distribute leaf hash chunks to each party
            let local_leaves = Net::recv_from_master(Some(leaf_chunks));

            // Also broadcast leaf_size
            Net::recv_from_master_uniform(Some(leaf_size));

            (f0, v0, full_evals, local_leaves, leaf_size)
        } else {
            // Workers receive their portion of leaf hashes
            let local_leaves: Vec<Byte32> = Net::recv_from_master(None);
            let leaf_size: usize = Net::recv_from_master_uniform(None);
            (vec![], vec![], vec![], local_leaves, leaf_size)
        };

    // Step 3: Each party builds local subtree
    let timer = start_timer!(|| "DCommit.DMerkle");
    let local_mt0 = MerkleTree::with_leaf_size(&local_leaves, leaf_size);
    let local_root = local_mt0.root();

    // Gather all local roots to master to build upper tree
    let all_roots_opt = Net::send_to_master(&local_root);
    end_timer!(timer);

    if Net::am_master() {
        let all_roots: Vec<Byte32> = all_roots_opt.unwrap();

        // Build upper tree from all party roots
        let upper_tree = MerkleTree::with_leaf_size(&all_roots, leaf_size);
        let rt0 = upper_tree.root();

        Ok((
            Some(DeepFoldCommitment { mu, rt0 }),
            DeepFoldProverCommitmentAdvice {
                f0,
                mt0: local_mt0,
                v0,
                f_tilde,
                upper_tree: Some(upper_tree),
            },
        ))
    } else {
        Ok((
            None,
            DeepFoldProverCommitmentAdvice {
                f0: vec![],
                mt0: local_mt0,
                v0: vec![],
                f_tilde: vec![],
                upper_tree: None,
            },
        ))
    }
}
