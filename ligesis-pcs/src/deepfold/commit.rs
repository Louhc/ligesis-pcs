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
    let all_evals_opt = Net::send_to_master(&poly.evaluations);

    // Step 2: Master computes full f0, v0
    let (f0, v0, f_tilde): (Vec<F>, Vec<F>, Vec<F>) = if Net::am_master() {
        let all_evals: Vec<Vec<F>> = all_evals_opt.unwrap();
        let full_evals: Vec<F> = all_evals.into_iter().flatten().collect();

        // Compute full coefficients and FFT
        let timer = start_timer!(|| "DDeepFold.Commit.FFT");
        let f0 = evals_to_coeffs(mu, &full_evals);
        let v0 = l0.fft(&f0);
        end_timer!(timer);

        (f0, v0, full_evals)
    } else {
        (vec![], vec![], vec![])
    };

    // Step 3: Master distributes v0 data to workers for parallel leaf hash computation
    // Leaf structure: leaf i contains v0[i], v0[i+step], ..., v0[i+(leaf_size-1)*step]
    // where step = v0.len() / leaf_size
    // We distribute so each worker computes a contiguous range of leaves
    let (local_v0_data, leaf_size, leaves_per_party): (Vec<F>, usize, usize) = if Net::am_master() {
        use super::utils::LEAF_SIZE;
        let leaf_size = LEAF_SIZE.min(v0.len());
        let step = v0.len() / leaf_size;  // number of leaves
        let leaves_per_party = step / num_party;

        // Reorganize v0 for distribution: for each worker, collect the elements needed for their leaves
        // Worker k computes leaves [k*leaves_per_party, (k+1)*leaves_per_party)
        // For leaf i, elements are at positions: i, i+step, i+2*step, ..., i+(leaf_size-1)*step
        let v0_chunks: Vec<Vec<F>> = (0..num_party)
            .map(|k| {
                let start_leaf = k * leaves_per_party;
                let end_leaf = (k + 1) * leaves_per_party;
                // Collect elements for this worker's leaves
                let mut chunk = Vec::with_capacity(leaves_per_party * leaf_size);
                for leaf_idx in start_leaf..end_leaf {
                    for j in 0..leaf_size {
                        chunk.push(v0[leaf_idx + j * step]);
                    }
                }
                chunk
            })
            .collect();

        let local_data: Vec<F> = Net::recv_from_master(Some(v0_chunks));
        Net::recv_from_master_uniform(Some((leaf_size, leaves_per_party)));
        (local_data, leaf_size, leaves_per_party)
    } else {
        let local_data: Vec<F> = Net::recv_from_master(None);
        let (leaf_size, leaves_per_party): (usize, usize) = Net::recv_from_master_uniform(None);
        (local_data, leaf_size, leaves_per_party)
    };

    // Step 4: Each party computes leaf hashes and builds local subtree
    let timer = start_timer!(|| "DDeepFold.Commit.LeafHash");
    // Compute leaf hashes from local v0 data (already organized by leaves)
    let local_leaves: Vec<Byte32> = (0..leaves_per_party)
        .map(|i| {
            let leaf: Vec<F> = (0..leaf_size)
                .map(|j| local_v0_data[i * leaf_size + j])
                .collect();
            compute_sha256_row(&leaf)
        })
        .collect();
    end_timer!(timer);

    let timer = start_timer!(|| "DDeepFold.Commit.Merkle");
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
