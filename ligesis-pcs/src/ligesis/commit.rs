//! Ligesis commit functions
//!
//! This module contains the commit implementations for Ligesis PCS:
//! - `ligesis_commit`: Standard commit
//! - `ligesis_d_commit`: Distributed commit

use crate::{
    deepfold::{*, multi_commit, split_polynomial},
    errors::PCSError, rscode::*, utils::*,
    PolynomialCommitmentScheme,
};
use arithmetic::math::Math;
use ark_ff::PrimeField;
use ark_poly::DenseMultilinearExtension;
use ark_std::{
    end_timer,
    start_timer,
    sync::Arc,
    vec::Vec,
};

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

use super::{
    LigeSISProverParam, LigeSISProverCommitmentAdvice, LigeSISCommitment,
    compute_sis_hash,
};

/// Standard Ligesis commit
#[allow(non_snake_case)]
pub fn ligesis_commit<F: PrimeField>(
    prover_param: &LigeSISProverParam<F>,
    poly: &Arc<DenseMultilinearExtension<F>>,
) -> Result<(LigeSISCommitment<F>, LigeSISProverCommitmentAdvice<F>), PCSError> {
    let &LigeSISProverParam {
        eta,
        s_lambda: _,
        mu,
        log_m,
        log_n,
        c: _,
        ref rs,
        ref mat_a,
        ref mat_a_pad,
        ref com_mat_a_advice,
        ref deepfold_prover_param,
    } = prover_param;
    let _ = (mat_a_pad, com_mat_a_advice);  // Suppress unused warnings
    let (m, n) = (1 << log_m, 1 << log_n);

    // Record original num_vars and pad if needed
    let num_vars = poly.num_vars;
    let poly_evals = if num_vars < mu {
        resize_eval(&poly.evaluations, mu)
    } else {
        poly.evaluations.clone()
    };
    let mat_f = reshape(&poly_evals, m, n);

    // encode `F`
    let timer = start_timer!(|| format!("Ligesis.Commit.RSEncode({}x{})", m, n));
    let mat_f_prime = mat_f.iter().map(|row| rs.encode(row)).collect::<Vec<_>>();
    end_timer!(timer);

    // compute `H`
    let timer = start_timer!(|| "Ligesis.Commit.SISHash");
    let mat_h = compute_sis_hash(mat_a, &mat_f_prime, eta, m);
    end_timer!(timer);

    // compute com(H) using multi_commit for potential splitting
    let timer = start_timer!(|| "Ligesis.Commit.DeepFold");
    let mat_h_poly = evals_to_arcpoly(&mat_h.concat());
    let (com_mat_h_list, com_mat_h_advices, _num_chunks_log2) =
        multi_commit(deepfold_prover_param, &mat_h_poly)?;

    // Get the chunk polynomials for later use in open
    let (mat_h_chunks, _) = split_polynomial(&mat_h_poly, deepfold_prover_param.max_mu);
    end_timer!(timer);

    Ok((
        LigeSISCommitment { num_vars, com_mat_h_list },
        LigeSISProverCommitmentAdvice {
            mat_f_prime,
            mat_h,
            mat_h_chunks,
            com_mat_h_advices,
        },
    ))
}

/// Distributed Ligesis commit
#[allow(non_snake_case)]
pub fn ligesis_d_commit<F: PrimeField>(
    prover_param: &LigeSISProverParam<F>,
    poly: &Arc<DenseMultilinearExtension<F>>,
) -> Result<(Option<LigeSISCommitment<F>>, LigeSISProverCommitmentAdvice<F>), PCSError> {
    let num_party = Net::n_parties();
    let num_party_vars = Net::n_parties().log_2() as usize;
    let party_id = Net::party_id();

    let &LigeSISProverParam {
        eta,
        s_lambda: _,
        mu: _,
        log_m,
        log_n,
        c,
        ref rs,
        ref mat_a,
        ref mat_a_pad,
        ref com_mat_a_advice,
        ref deepfold_prover_param,
    } = prover_param;
    let _ = (mat_a_pad, com_mat_a_advice);  // Suppress unused warnings
    let (m, n) = (1 << log_m, 1 << log_n);

    // Record original num_vars (for distributed case, actual num_vars = poly.num_vars + num_party_vars)
    let num_vars = poly.num_vars + num_party_vars;
    let mat_f = reshape(&poly.evaluations, m / num_party, n);

    // encode `F`
    let timer = start_timer!(|| format!("DLigesis.Commit.RSEncode({}x{})", m / num_party, n));
    let mat_f_prime = mat_f.iter().map(|row| rs.encode(row)).collect::<Vec<_>>();
    end_timer!(timer);

    // compute `H`
    let mat_a_k: Vec<Vec<F>> = mat_a.iter().map(
        |row| row[party_id * eta * m / num_party..(party_id + 1) * eta * m / num_party].to_vec()
    ).collect();

    let timer = start_timer!(|| format!("DLigesis.Commit.SISHash({}x{}x{})", c, m / num_party, n * 2));
    let mat_h_i = compute_sis_hash(&mat_a_k, &mat_f_prime, eta, m / num_party);
    end_timer!(timer);

    let all_mat_h = Net::send_to_master(&mat_h_i);

    // Master computes full mat_h and distributes portions to workers for d_commit
    let mat_h = if Net::am_master() {
        let all_mat_h = all_mat_h.ok_or(PCSError::UnexpectedNone("all_mat_h".into()))?;
        (0..c)
            .map(|i| {
                (0..2 * n)
                    .map(|j| (0..num_party).map(|k| all_mat_h[k][i][j]).sum::<F>())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    } else {
        vec![]
    };

    // Split mat_h into chunks on master
    let mat_h_poly = if Net::am_master() {
        evals_to_arcpoly(&mat_h.concat())
    } else {
        evals_to_arcpoly(&vec![F::ZERO])  // Dummy for non-master
    };

    let (mat_h_chunks_master, num_chunks_log2) = if Net::am_master() {
        split_polynomial(&mat_h_poly, deepfold_prover_param.max_mu)
    } else {
        (vec![], 0)
    };
    let num_chunks = if Net::am_master() { mat_h_chunks_master.len() } else { 0 };

    // Broadcast number of chunks to all parties
    let num_chunks: usize = if Net::am_master() {
        Net::recv_from_master_uniform(Some(num_chunks));
        num_chunks
    } else {
        Net::recv_from_master_uniform(None)
    };

    // Each party gets their portion of each chunk for d_commit
    let local_chunk_size = (1 << deepfold_prover_param.max_mu) / num_party;

    let timer = start_timer!(|| "DLigesis.Commit.DeepFold");
    let mut com_mat_h_list = Vec::with_capacity(num_chunks);
    let mut com_mat_h_advices = Vec::with_capacity(num_chunks);
    let mut mat_h_chunks_local = Vec::with_capacity(num_chunks);

    for chunk_idx in 0..num_chunks {
        // Distribute this chunk's data
        let local_chunk: Vec<F> = if Net::am_master() {
            let chunk_full = &mat_h_chunks_master[chunk_idx].evaluations;
            let chunks_for_parties: Vec<Vec<F>> = (0..num_party)
                .map(|k| chunk_full[k * local_chunk_size..(k + 1) * local_chunk_size].to_vec())
                .collect();
            Net::recv_from_master(Some(chunks_for_parties))
        } else {
            Net::recv_from_master(None)
        };

        let local_chunk_poly = evals_to_arcpoly(&local_chunk);
        let (com_opt, advice) = DeepFoldPCS::d_commit(deepfold_prover_param, &local_chunk_poly)?;

        if Net::am_master() {
            com_mat_h_list.push(com_opt.unwrap());
        }
        com_mat_h_advices.push(advice);
        mat_h_chunks_local.push(local_chunk_poly);
    }
    end_timer!(timer);

    if Net::am_master() {
        Ok((
            Some(LigeSISCommitment { num_vars, com_mat_h_list }),
            LigeSISProverCommitmentAdvice {
                mat_f_prime,
                mat_h,
                mat_h_chunks: mat_h_chunks_local,
                com_mat_h_advices,
            },
        ))
    } else {
        Ok((
            None,
            LigeSISProverCommitmentAdvice {
                mat_f_prime,
                mat_h: mat_h_i,
                mat_h_chunks: mat_h_chunks_local,
                com_mat_h_advices,
            },
        ))
    }
}
