//! Ligesis commit functions
//!
//! This module contains the commit implementations for Ligesis PCS:
//! - `ligesis_commit`: Standard commit
//! - `ligesis_d_commit`: Distributed commit

use crate::{
    deepfold::*, errors::PCSError, rscode::*, utils::*,
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
    let timer = start_timer!(|| "Commit.RS");
    let mat_f_prime = mat_f.iter().map(|row| rs.encode(row)).collect::<Vec<_>>();
    end_timer!(timer);

    // compute `H`
    let timer = start_timer!(|| "Commit.SISHash");
    let mat_h = compute_sis_hash(mat_a, &mat_f_prime, eta, m);
    end_timer!(timer);

    // compute com(H)
    let mat_h_pad =
        evals_to_arcpoly(&resize_eval(&mat_h.concat(), deepfold_prover_param.max_mu));
    let (com_mat_h, com_mat_h_advice) = DeepFoldPCS::commit(deepfold_prover_param, &mat_h_pad)?;

    Ok((
        LigeSISCommitment { num_vars, com_mat_h },
        LigeSISProverCommitmentAdvice {
            mat_f_prime,
            mat_h,
            mat_h_pad,
            com_mat_h_advice,
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
    let mat_f_prime = mat_f.iter().map(|row| rs.encode(row)).collect::<Vec<_>>();

    // compute `H`
    let mat_a_k: Vec<Vec<F>> = mat_a.iter().map(
        |row| row[party_id * eta * m / num_party..(party_id + 1) * eta * m / num_party].to_vec()
    ).collect();

    let timer = start_timer!(|| format!("Commit.DistributedSIS({}x{}x{})", c, m / num_party * eta, n * 2));
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

    // Distribute mat_h_pad portions for distributed DeepFold commit
    let mat_h_full = if Net::am_master() {
        resize_eval(&mat_h.concat(), deepfold_prover_param.max_mu)
    } else {
        vec![]
    };

    // Each party gets their portion of mat_h_pad for d_commit
    let local_mat_h_size = (1 << deepfold_prover_param.max_mu) / num_party;
    let local_mat_h: Vec<F> = if Net::am_master() {
        let chunks: Vec<Vec<F>> = (0..num_party)
            .map(|k| mat_h_full[k * local_mat_h_size..(k + 1) * local_mat_h_size].to_vec())
            .collect();
        Net::recv_from_master(Some(chunks))
    } else {
        Net::recv_from_master(None)
    };

    // All parties call d_commit
    let mat_h_pad = evals_to_arcpoly(&local_mat_h);
    let (com_mat_h_opt, com_mat_h_advice) =
        DeepFoldPCS::d_commit(deepfold_prover_param, &mat_h_pad)?;

    if Net::am_master() {
        Ok((
            Some(LigeSISCommitment { num_vars, com_mat_h: com_mat_h_opt.unwrap() }),
            LigeSISProverCommitmentAdvice {
                mat_f_prime,
                mat_h,
                mat_h_pad,
                com_mat_h_advice,
            },
        ))
    } else {
        Ok((
            None,
            LigeSISProverCommitmentAdvice {
                mat_f_prime,
                mat_h: mat_h_i,
                mat_h_pad,
                com_mat_h_advice,
            },
        ))
    }
}
