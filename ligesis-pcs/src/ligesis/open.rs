//! Ligesis open functions
//!
//! This module contains the open implementations for Ligesis PCS:
//! - `ligesis_open`: Extension field SumCheck with 128-bit security
//! - `ligesis_d_open`: Distributed open with 128-bit security

use crate::{
    deepfold::{*, split_point_for_chunks, expand_for_batch_open_ext},
    errors::PCSError, rscode::*, utils::*,
    ext_sumcheck::ExtSumCheckBuilder,
    types::{HasQuadraticExtension, FieldExtension},
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
use transcript::IOPTranscript;

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

use super::{
    LigeSISProverParam, LigeSISProverCommitmentAdvice, LigeSISProof,
    ExtSumCheckWithReductionProof,
};

/// Ligesis open with extension field SumCheck (128-bit security)
/// Uses extension field SumCheck and direct extension field opening in DeepFold
#[allow(non_snake_case)]
pub fn ligesis_open<F: PrimeField + HasQuadraticExtension>(
    prover_param: &LigeSISProverParam<F>,
    poly: &Arc<DenseMultilinearExtension<F>>,
    advice: &LigeSISProverCommitmentAdvice<F>,
    point: &[F],
    transcript: &mut IOPTranscript<F>,
) -> Result<LigeSISProof<F>, PCSError> {
    let LigeSISProverParam {
        eta,
        s_lambda,
        mu,
        log_m,
        log_n,
        c,
        ref rs,
        ref mat_a,
        ref mat_a_pad,
        ref com_mat_a_advice,
        ref deepfold_prover_param,
    } = *prover_param;
    let (m, n) = (1 << log_m, 1 << log_n);
    let rs_len = rs.get_k();
    let log_rs_len = rs_len.ilog2() as usize;

    assert_eq!(mu, log_m + log_n);
    assert!(poly.num_vars <= mu);

    // Pad polynomial and point if needed
    let poly_evals = if poly.num_vars < mu {
        resize_eval(&poly.evaluations, mu)
    } else {
        poly.evaluations.clone()
    };
    let point = resize_point(&point.to_vec(), mu);

    let mat_f = reshape(&poly_evals, m, n);

    let LigeSISProverCommitmentAdvice {
        mat_f_prime,
        mat_h: _,
        mat_h_chunks,
        com_mat_h_advices,
    } = advice;

    // Step 1
    let (z1, z2) = (point[log_n..].to_vec(), point[..log_n].to_vec());
    let eq_z1 = get_tensor(&z1);

    // Step 2: Commit to a
    let a: Vec<F> = (0..n)
        .map(|j| (0..m).map(|i| eq_z1[i] * mat_f[i][j]).sum())
        .collect();
    let a_pad = evals_to_arcpoly(&resize_eval(&a, deepfold_prover_param.max_mu));
    let (com_a, com_a_advice) = DeepFoldPCS::commit(deepfold_prover_param, &a_pad)?;

    // Step 3
    let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;

    // Step 4: Commit to bI
    let mat_f_prime_trans = transposition(mat_f_prime);
    let mat_bI = transposition(
        &I.iter()
            .map(|&i| decompose_vector(&mat_f_prime_trans[i]))
            .collect::<Vec<_>>(),
    );
    let bI_field = bool_vec_to_field_vec(&mat_bI.concat());
    let bI_poly = evals_to_arcpoly(&bI_field);
    let (com_bI_list, com_bI_advices, _bI_num_chunks_log2) =
        multi_commit(deepfold_prover_param, &bI_poly)?;
    let (bI_chunks, _) = split_polynomial(&bI_poly, deepfold_prover_param.max_mu);

    // Step 5: Get challenges
    let alpha1 = transcript
        .get_and_append_challenge_vectors(b"alpha1", (m * eta * s_lambda).ilog2() as usize)?;
    let alpha2 = transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
    let alpha3 = transcript.get_and_append_challenge_vectors(b"alpha3", log_rs_len)?;

    // Step 6: Extension field SumCheck for bI check (no reduction needed)
    let timer = start_timer!(|| "Ligesis.Open.ExtSumchecks");
    let bI_field_minus_one: Vec<F> = bI_field.iter().map(|&x| x - F::ONE).collect();
    let tensor_alpha1 = get_tensor(&alpha1);

    let bI_check_proof = run_ext_sumcheck::<F>(
        bI_field.len().ilog2() as usize,
        vec![&bI_field[..], &bI_field_minus_one[..], &tensor_alpha1[..]],
        F::ONE,
        transcript,
    )?;
    let r1_ext = bI_check_proof.ext_proof.point.clone();

    // Step 7: Check rs_a
    let rs_a = rs.encode(&a);
    let g = rs.get_generator();
    let rs_a_pad = evals_to_arcpoly(&resize_eval(&rs_a, deepfold_prover_param.max_mu));
    let (com_rs_a, com_rs_a_advice) = DeepFoldPCS::commit(deepfold_prover_param, &rs_a_pad)?;

    // Step 7.1: Extension field SumCheck for rs_a check
    let alpha3_mat_g = compute_alpha_mat_g(log_rs_len, log_n, &g, &alpha3);
    let alpha3_mat_g_n = alpha3_mat_g[log_rs_len][..n].to_vec();

    let rs_a_check_proof = run_ext_sumcheck::<F>(
        log_n,
        vec![&alpha3_mat_g_n[..], &a[..]],
        F::ONE,
        transcript,
    )?;
    let r6_ext = rs_a_check_proof.ext_proof.point.clone();
    // Get base field version for subsequent computations
    let r6: Vec<F> = r6_ext.iter().map(|x| F::ext_real(x)).collect();

    // Step 7.2: Extension field SumChecks for mat_g checks
    let mut cur_p = vec![r6.clone(), vec![F::ZERO; log_rs_len - log_n]].concat();
    let mut mat_g_check_proofs = Vec::new();
    for i in (2..=log_rs_len).rev() {
        let (x, b) = (cur_p[..i - 1].to_vec(), cur_p[i - 1]);
        let gi = g.pow([1u64 << (log_rs_len - i)]);
        let w: Vec<F> = (0..1 << (i - 1))
            .map(|z| {
                F::ONE - alpha3[log_rs_len - i]
                    + alpha3[log_rs_len - i]
                        * (gi.pow([z]) * (F::ONE - b) + gi.pow([z + (1 << (i - 1))]) * b)
            })
            .collect();
        let tensor_x = get_tensor(&x);

        let mat_g_check_proof = run_ext_sumcheck::<F>(
            i - 1,
            vec![&tensor_x[..], &alpha3_mat_g[i - 1][..], &w[..]],
            F::ONE,
            transcript,
        )?;

        // Use real part of extension field point for next iteration
        cur_p = mat_g_check_proof.ext_proof.point.iter().map(|x| F::ext_real(x)).collect();
        mat_g_check_proofs.push(mat_g_check_proof);
    }

    // Step 8
    let v = otimes(
        &get_tensor(&z1),
        &(0..eta)
            .map(|i| F::from(2u64).pow([i as u64]))
            .collect::<Vec<_>>(),
    );

    let r2 = transcript.get_and_append_challenge_vectors(b"r2", s_lambda.ilog2() as usize)?;
    let r3 = transcript.get_and_append_challenge_vectors(b"r3", (2 * n).ilog2() as usize)?;

    // Step 9: Extension field SumCheck for alpha2_a_bI_r2 check
    let alpha2_a = mat_mul(&vec![get_tensor(&alpha2)], mat_a)[0].clone();
    let bI_r2 =
        field_mat_mul_bool_mat(&vec![get_tensor(&r2)], &transposition(&mat_bI))[0].clone();

    let alpha2_a_bI_r2_check_proof = run_ext_sumcheck::<F>(
        mat_bI.len().ilog2() as usize,
        vec![&alpha2_a[..], &bI_r2[..]],
        F::ONE,
        transcript,
    )?;
    let r4_ext = alpha2_a_bI_r2_check_proof.ext_proof.point.clone();

    // Step 10: Extension field SumCheck for v_bI_r2 check
    let v_bI_r2_check_proof = run_ext_sumcheck::<F>(
        v.len().ilog2() as usize,
        vec![&v[..], &bI_r2[..]],
        F::ONE,
        transcript,
    )?;
    let r5_ext = v_bI_r2_check_proof.ext_proof.point.clone();
    end_timer!(timer);

    // Step 11: DeepFold batch open at extension field points
    // Convert base field points to extension field
    let z2_ext: Vec<F::Extension> = z2.iter().map(|&x| F::Extension::from_base(x)).collect();
    let r3_ext: Vec<F::Extension> = r3.iter().map(|&x| F::Extension::from_base(x)).collect();
    let alpha2_ext: Vec<F::Extension> = alpha2.iter().map(|&x| F::Extension::from_base(x)).collect();
    let alpha3_ext: Vec<F::Extension> = alpha3.iter().map(|&x| F::Extension::from_base(x)).collect();
    let r2_ext: Vec<F::Extension> = r2.iter().map(|&x| F::Extension::from_base(x)).collect();

    // Full point for mat_h is (r3, alpha2)
    let mat_h_full_point_ext: Vec<F::Extension> = vec![r3_ext.clone(), alpha2_ext.clone()].concat();

    // Expand mat_h chunks for batch_open
    let (mat_h_chunk_polys, mat_h_chunk_advices, mat_h_chunk_points) =
        expand_for_batch_open_ext::<F>(
            mat_h_chunks,
            com_mat_h_advices,
            &mat_h_full_point_ext,
            deepfold_prover_param.max_mu,
        );

    // Expand bI chunks for each opening point
    // bI is opened at 3 points: r1, (r2, r4), (r2, r5)
    let (bI_chunk_polys_r1, bI_chunk_advices_r1, bI_chunk_points_r1) =
        expand_for_batch_open_ext::<F>(
            &bI_chunks,
            &com_bI_advices,
            &r1_ext,
            deepfold_prover_param.max_mu,
        );
    let bI_r2_r4_point: Vec<F::Extension> = vec![r2_ext.clone(), r4_ext.clone()].concat();
    let (bI_chunk_polys_r2r4, bI_chunk_advices_r2r4, bI_chunk_points_r2r4) =
        expand_for_batch_open_ext::<F>(
            &bI_chunks,
            &com_bI_advices,
            &bI_r2_r4_point,
            deepfold_prover_param.max_mu,
        );
    let bI_r2_r5_point: Vec<F::Extension> = vec![r2_ext.clone(), r5_ext.clone()].concat();
    let (bI_chunk_polys_r2r5, bI_chunk_advices_r2r5, bI_chunk_points_r2r5) =
        expand_for_batch_open_ext::<F>(
            &bI_chunks,
            &com_bI_advices,
            &bI_r2_r5_point,
            deepfold_prover_param.max_mu,
        );

    // Build polys, advices, points vectors with expanded mat_h and bI chunks
    let mut polys: Vec<Arc<DenseMultilinearExtension<F>>> = vec![
        Arc::clone(&a_pad),
        Arc::clone(&a_pad),
        Arc::clone(&rs_a_pad),
        Arc::clone(&rs_a_pad),
    ];
    polys.extend(mat_h_chunk_polys);
    polys.push(Arc::clone(mat_a_pad));
    polys.extend(bI_chunk_polys_r1);
    polys.extend(bI_chunk_polys_r2r4);
    polys.extend(bI_chunk_polys_r2r5);

    let mut advices: Vec<&DeepFoldProverCommitmentAdvice<F>> = vec![
        &com_a_advice,
        &com_a_advice,
        &com_rs_a_advice,
        &com_rs_a_advice,
    ];
    for advice in &mat_h_chunk_advices {
        advices.push(advice);
    }
    advices.push(com_mat_a_advice);
    for advice in &bI_chunk_advices_r1 {
        advices.push(advice);
    }
    for advice in &bI_chunk_advices_r2r4 {
        advices.push(advice);
    }
    for advice in &bI_chunk_advices_r2r5 {
        advices.push(advice);
    }

    let timer = start_timer!(|| "Ligesis.Open.DeepFold");
    let mut points_ext: Vec<Vec<F::Extension>> = vec![
        resize_point_ext::<F>(&z2_ext, deepfold_prover_param.max_mu),
        resize_point_ext::<F>(&r6_ext, deepfold_prover_param.max_mu),
        resize_point_ext::<F>(&r3_ext, deepfold_prover_param.max_mu),
        resize_point_ext::<F>(&alpha3_ext, deepfold_prover_param.max_mu),
    ];
    for point in mat_h_chunk_points {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }
    points_ext.push(resize_point_ext::<F>(&vec![r4_ext.clone(), alpha2_ext.clone()].concat(), deepfold_prover_param.max_mu));
    for point in bI_chunk_points_r1 {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }
    for point in bI_chunk_points_r2r4 {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }
    for point in bI_chunk_points_r2r5 {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }

    let deepfold_batched_proof = DeepFoldPCS::batch_open_at_ext_point(
        deepfold_prover_param,
        polys,
        &advices,
        &points_ext,
        transcript,
    )?;
    end_timer!(timer);

    Ok(LigeSISProof {
        com_a,
        com_bI_list,
        com_rs_a,
        bI_check_proof,
        alpha2_a_bI_r2_check_proof,
        v_bI_r2_check_proof,
        rs_a_check_proof,
        mat_g_check_proofs,
        deepfold_batched_proof,
    })
}

/// Helper function to run extension field SumCheck (no reduction needed)
#[allow(non_snake_case)]
fn run_ext_sumcheck<F: PrimeField + HasQuadraticExtension>(
    num_vars: usize,
    evals_list: Vec<&[F]>,
    coeff: F,
    transcript: &mut IOPTranscript<F>,
) -> Result<ExtSumCheckWithReductionProof<F>, PCSError> {
    // Build and run extension field SumCheck
    let mut builder = ExtSumCheckBuilder::<F, F::Extension>::new(num_vars);
    let mles: Vec<Arc<DenseMultilinearExtension<F>>> = evals_list
        .iter()
        .map(|evals| evals_to_arcpoly(&evals.to_vec()))
        .collect();
    builder = builder.add_mle_list(mles, coeff)?;
    let ext_proof = builder.prove(transcript)?;

    Ok(ExtSumCheckWithReductionProof { ext_proof })
}

/// Helper to resize extension field point
fn resize_point_ext<F: PrimeField + HasQuadraticExtension>(point: &[F::Extension], target_len: usize) -> Vec<F::Extension> {
    let mut result = point.to_vec();
    while result.len() < target_len {
        result.push(F::Extension::from_base(F::ZERO));
    }
    result
}

/// Distributed Ligesis open with 128-bit security
/// Uses extension field SumCheck and direct extension field opening in DeepFold
#[allow(non_snake_case)]
pub fn ligesis_d_open<F: PrimeField + HasQuadraticExtension>(
    prover_param: &LigeSISProverParam<F>,
    poly: &Arc<DenseMultilinearExtension<F>>,
    advice: &LigeSISProverCommitmentAdvice<F>,
    point: &[F],
    transcript: &mut IOPTranscript<F>,
) -> Result<Option<LigeSISProof<F>>, PCSError> {
    let &LigeSISProverParam {
        eta,
        s_lambda,
        mu,
        log_m,
        log_n,
        c,
        ref rs,
        ref mat_a,
        ref mat_a_pad,
        com_mat_a_advice: _,
        ref deepfold_prover_param,
    } = prover_param;
    let num_party = Net::n_parties();
    let num_party_vars = Net::n_parties().log_2() as usize;
    let (m, n) = (1 << log_m, 1 << log_n);
    let rs_len = rs.get_k();
    let log_rs_len = rs_len.ilog2() as usize;
    let local_poly_size = (1 << deepfold_prover_param.max_mu) / num_party;

    assert_eq!(mu, log_m + log_n);
    assert!(poly.num_vars <= mu - num_party_vars);

    // Pad point if needed
    let point = resize_point(&point.to_vec(), mu);

    let mat_f = reshape(&poly.evaluations, m / num_party, n);

    let LigeSISProverCommitmentAdvice {
        mat_f_prime,
        mat_h: _,
        mat_h_chunks,
        com_mat_h_advices,
    } = advice;

    // Step 1
    let (z1, z2) = (point[log_n..].to_vec(), point[..log_n].to_vec());
    let (z1_0, z1_1) = (
        z1[log_m - num_party_vars..].to_vec(),
        z1[..log_m - num_party_vars].to_vec(),
    );
    let eq_z1_1 = get_tensor(&z1_1);
    let eq_z1_0 = get_tensor(&z1_0);

    // Step 2: Compute a and d_commit
    let timer = start_timer!(|| "DLigesis.Open.CommitA");
    let a_k = (0..n)
        .map(|j| {
            (0..m / num_party)
                .map(|i| eq_z1_1[i] * mat_f[i][j])
                .sum()
        })
        .collect::<Vec<F>>();
    let a_k_list = Net::send_to_master(&a_k);

    // Master computes full a and distributes for d_commit
    let (a_full, a_saved): (Vec<F>, Vec<F>) = if Net::am_master() {
        let a_k_list = a_k_list.ok_or(PCSError::UnexpectedNone("a_k_list".into()))?;
        let a: Vec<F> = (0..n)
            .map(|j| (0..num_party).map(|k| eq_z1_0[k] * a_k_list[k][j]).sum())
            .collect();
        (resize_eval(&a, deepfold_prover_param.max_mu), a)
    } else {
        (vec![], vec![])
    };

    // Distribute a for d_commit
    let local_a: Vec<F> = if Net::am_master() {
        let chunks: Vec<Vec<F>> = (0..num_party)
            .map(|k| a_full[k * local_poly_size..(k + 1) * local_poly_size].to_vec())
            .collect();
        Net::recv_from_master(Some(chunks))
    } else {
        Net::recv_from_master(None)
    };

    let a_pad = evals_to_arcpoly(&local_a);
    let (com_a_opt, com_a_advice) = DeepFoldPCS::d_commit(deepfold_prover_param, &a_pad)?;
    let com_a = if Net::am_master() { com_a_opt.unwrap() } else { DeepFoldCommitment::default() };
    end_timer!(timer);

    // Step 3: receive challenge indices
    let I = if Net::am_master() {
        let I = transcript.get_and_append_challenge_indices(b"I", s_lambda, 2 * n)?;
        Net::recv_from_master_uniform(Some(I.clone()));
        I
    } else {
        Net::recv_from_master_uniform(None)
    };

    // Step 4: Compute bI and d_multi_commit
    let timer = start_timer!(|| "DLigesis.Open.CommitBI");
    let mat_bI_k = {
        let mat_f_prime_trans = transposition(&mat_f_prime);
        transposition(
            &I.iter()
                .map(|&i| decompose_vector(&mat_f_prime_trans[i]))
                .collect::<Vec<_>>(),
        )
    };
    let mat_bI_k_list = Net::send_to_master(&mat_bI_k);

    // On master: compute full bI and split into chunks
    let (mat_bI, bI_chunks_master, bI_num_chunks) = if Net::am_master() {
        let mat_bI_k_list = mat_bI_k_list.ok_or(PCSError::UnexpectedNone("mat_bI_k_list".into()))?;
        let mat_bI = mat_bI_k_list.concat();
        let bI_field = bool_vec_to_field_vec(&mat_bI.concat());
        let bI_poly = evals_to_arcpoly(&bI_field);
        let (chunks, _) = split_polynomial(&bI_poly, deepfold_prover_param.max_mu);
        let num_chunks = chunks.len();
        (mat_bI, chunks, num_chunks)
    } else {
        (vec![], vec![], 0)
    };

    // Broadcast number of bI chunks to all parties
    let bI_num_chunks: usize = if Net::am_master() {
        Net::recv_from_master_uniform(Some(bI_num_chunks))
    } else {
        Net::recv_from_master_uniform(None)
    };

    // d_commit each bI chunk
    let mut com_bI_list = Vec::with_capacity(bI_num_chunks);
    let mut com_bI_advices = Vec::with_capacity(bI_num_chunks);
    let mut bI_chunks_local = Vec::with_capacity(bI_num_chunks);

    for chunk_idx in 0..bI_num_chunks {
        // Distribute this chunk's data
        let local_chunk: Vec<F> = if Net::am_master() {
            let chunk_full = &bI_chunks_master[chunk_idx].evaluations;
            let chunks_for_parties: Vec<Vec<F>> = (0..num_party)
                .map(|k| chunk_full[k * local_poly_size..(k + 1) * local_poly_size].to_vec())
                .collect();
            Net::recv_from_master(Some(chunks_for_parties))
        } else {
            Net::recv_from_master(None)
        };

        let local_chunk_poly = evals_to_arcpoly(&local_chunk);
        let (com_opt, advice) = DeepFoldPCS::d_commit(deepfold_prover_param, &local_chunk_poly)?;

        if Net::am_master() {
            com_bI_list.push(com_opt.unwrap());
        }
        com_bI_advices.push(advice);
        bI_chunks_local.push(local_chunk_poly);
    }
    end_timer!(timer);

    // Step 5: receive challenge vectors
    let (alpha1, alpha2, alpha3) = if Net::am_master() {
        let alpha1 = transcript.get_and_append_challenge_vectors(
            b"alpha1",
            (m * eta * s_lambda).ilog2() as usize,
        )?;
        let alpha2 =
            transcript.get_and_append_challenge_vectors(b"alpha2", c.ilog2() as usize)?;
        let alpha3 = transcript.get_and_append_challenge_vectors(b"alpha3", log_rs_len)?;
        Net::recv_from_master_uniform(Some((alpha1.clone(), alpha2.clone(), alpha3.clone())));
        (alpha1, alpha2, alpha3)
    } else {
        Net::recv_from_master_uniform(None)
    };

    // Step 6: Extension field SumCheck for bI check (on master)
    let timer = start_timer!(|| "DLigesis.Open.ExtSumchecks");
    let bI_field = if Net::am_master() {
        bool_vec_to_field_vec(&mat_bI.concat())
    } else {
        vec![]
    };
    let (bI_check_proof, r1_ext) = if Net::am_master() {
        let bI_field_minus_one: Vec<F> = bI_field.iter().map(|&x| x - F::ONE).collect();
        let tensor_alpha1 = get_tensor(&alpha1);
        let proof = run_ext_sumcheck::<F>(
            bI_field.len().ilog2() as usize,
            vec![&bI_field[..], &bI_field_minus_one[..], &tensor_alpha1[..]],
            F::ONE,
            transcript,
        )?;
        let r1 = proof.ext_proof.point.clone();
        (proof, r1)
    } else {
        (ExtSumCheckWithReductionProof { ext_proof: crate::ext_sumcheck::ExtSumCheckProof::default() }, vec![])
    };
    end_timer!(timer);

    // Step 7: Compute rs_a and d_commit
    let timer = start_timer!(|| "DLigesis.Open.CommitRSA");
    let (rs_a_full, rs_a_check_proof, r6_ext, mat_g_check_proofs) = if Net::am_master() {
        let rs_a = rs.encode(&a_saved);
        let g = rs.get_generator();

        // Step 7.1: Extension field SumCheck for rs_a check
        let alpha3_mat_g = compute_alpha_mat_g(log_rs_len as usize, log_n, &g, &alpha3);
        let alpha3_mat_g_n = alpha3_mat_g[log_rs_len][..n].to_vec();
        let rs_a_check_proof = run_ext_sumcheck::<F>(
            log_n,
            vec![&alpha3_mat_g_n[..], &a_saved[..]],
            F::ONE,
            transcript,
        )?;
        let r6_ext = rs_a_check_proof.ext_proof.point.clone();
        let r6: Vec<F> = r6_ext.iter().map(|x| F::ext_real(x)).collect();

        // Step 7.2: Extension field SumChecks for mat_g checks
        let mut cur_p = vec![r6.clone(), vec![F::ZERO; log_rs_len - log_n]].concat();
        let mut mat_g_check_proofs = Vec::new();
        for i in (2..=log_rs_len).rev() {
            let (x, b) = (cur_p[..i - 1].to_vec(), cur_p[i - 1]);
            let gi = g.pow([1u64 << (log_rs_len - i)]);
            let w: Vec<F> = (0..1 << (i - 1))
                .map(|z| {
                    F::ONE - alpha3[log_rs_len - i]
                        + alpha3[log_rs_len - i]
                            * (gi.pow([z]) * (F::ONE - b) + gi.pow([z + (1 << (i - 1))]) * b)
                })
                .collect();
            let tensor_x = get_tensor(&x);
            let mat_g_check_proof = run_ext_sumcheck::<F>(
                i - 1,
                vec![&tensor_x[..], &alpha3_mat_g[i - 1][..], &w[..]],
                F::ONE,
                transcript,
            )?;
            cur_p = mat_g_check_proof.ext_proof.point.iter().map(|x| F::ext_real(x)).collect();
            mat_g_check_proofs.push(mat_g_check_proof);
        }

        let rs_a_full = resize_eval(&rs_a, deepfold_prover_param.max_mu);
        (rs_a_full, rs_a_check_proof, r6_ext, mat_g_check_proofs)
    } else {
        (vec![], ExtSumCheckWithReductionProof { ext_proof: crate::ext_sumcheck::ExtSumCheckProof::default() }, vec![], vec![])
    };

    // Distribute rs_a for d_commit
    let local_rs_a: Vec<F> = if Net::am_master() {
        let chunks: Vec<Vec<F>> = (0..num_party)
            .map(|k| rs_a_full[k * local_poly_size..(k + 1) * local_poly_size].to_vec())
            .collect();
        Net::recv_from_master(Some(chunks))
    } else {
        Net::recv_from_master(None)
    };

    let rs_a_pad = evals_to_arcpoly(&local_rs_a);
    let (com_rs_a_opt, com_rs_a_advice) = DeepFoldPCS::d_commit(deepfold_prover_param, &rs_a_pad)?;
    let com_rs_a = if Net::am_master() { com_rs_a_opt.unwrap() } else { DeepFoldCommitment::default() };
    end_timer!(timer);

    // Step 8
    let (r2, r3, v, alpha2_a_bI_r2_check_proof, r4_ext, v_bI_r2_check_proof, r5_ext) = if Net::am_master() {
        let v = otimes(
            &get_tensor(&z1),
            &(0..eta)
                .map(|i| F::from(2u64).pow([i as u64]))
                .collect::<Vec<_>>(),
        );
        let r2 = transcript.get_and_append_challenge_vectors(b"r2", s_lambda.ilog2() as usize)?;
        let r3 = transcript.get_and_append_challenge_vectors(b"r3", (2 * n).ilog2() as usize)?;

        // Step 9: Extension field SumCheck for alpha2_a_bI_r2 check
        let alpha2_a = mat_mul(&vec![get_tensor(&alpha2)], &mat_a)[0].clone();
        let bI_r2 =
            field_mat_mul_bool_mat(&vec![get_tensor(&r2)], &transposition(&mat_bI))[0].clone();
        let alpha2_a_bI_r2_check_proof = run_ext_sumcheck::<F>(
            mat_bI.len().ilog2() as usize,
            vec![&alpha2_a[..], &bI_r2[..]],
            F::ONE,
            transcript,
        )?;
        let r4_ext = alpha2_a_bI_r2_check_proof.ext_proof.point.clone();

        // Step 10: Extension field SumCheck for v_bI_r2 check
        let v_bI_r2_check_proof = run_ext_sumcheck::<F>(
            v.len().ilog2() as usize,
            vec![&v[..], &bI_r2[..]],
            F::ONE,
            transcript,
        )?;
        let r5_ext = v_bI_r2_check_proof.ext_proof.point.clone();

        (r2, r3, v, alpha2_a_bI_r2_check_proof, r4_ext, v_bI_r2_check_proof, r5_ext)
    } else {
        (vec![], vec![], vec![], ExtSumCheckWithReductionProof { ext_proof: crate::ext_sumcheck::ExtSumCheckProof::default() }, vec![], ExtSumCheckWithReductionProof { ext_proof: crate::ext_sumcheck::ExtSumCheckProof::default() }, vec![])
    };

    // Distribute mat_a_pad for d_batch_open
    let local_mat_a: Vec<F> = if Net::am_master() {
        let mat_a_full = &mat_a_pad.evaluations;
        let chunks: Vec<Vec<F>> = (0..num_party)
            .map(|k| mat_a_full[k * local_poly_size..(k + 1) * local_poly_size].to_vec())
            .collect();
        Net::recv_from_master(Some(chunks))
    } else {
        Net::recv_from_master(None)
    };
    let local_mat_a_pad = evals_to_arcpoly(&local_mat_a);

    // d_commit mat_a for all parties to have proper advice
    let (_, com_mat_a_advice_dist) = DeepFoldPCS::d_commit(deepfold_prover_param, &local_mat_a_pad)?;

    // Compute and broadcast extension field points for d_batch_open
    // Including expanded mat_h chunk points
    let z2_ext: Vec<F::Extension> = z2.iter().map(|&x| F::Extension::from_base(x)).collect();
    let r3_ext: Vec<F::Extension> = r3.iter().map(|&x| F::Extension::from_base(x)).collect();
    let alpha2_ext: Vec<F::Extension> = alpha2.iter().map(|&x| F::Extension::from_base(x)).collect();
    let alpha3_ext: Vec<F::Extension> = alpha3.iter().map(|&x| F::Extension::from_base(x)).collect();
    let r2_ext: Vec<F::Extension> = r2.iter().map(|&x| F::Extension::from_base(x)).collect();

    // Full point for mat_h is (r3, alpha2)
    let mat_h_full_point_ext: Vec<F::Extension> = vec![r3_ext.clone(), alpha2_ext.clone()].concat();

    // Expand mat_h chunks for batch_open (works for all parties since they have their local chunks)
    let (mat_h_chunk_polys, mat_h_chunk_advices, mat_h_chunk_points) =
        expand_for_batch_open_ext::<F>(
            mat_h_chunks,
            com_mat_h_advices,
            &mat_h_full_point_ext,
            deepfold_prover_param.max_mu,
        );

    // Expand bI chunks for each opening point
    // bI is opened at 3 points: r1, (r2, r4), (r2, r5)
    let (bI_chunk_polys_r1, bI_chunk_advices_r1, bI_chunk_points_r1) =
        expand_for_batch_open_ext::<F>(
            &bI_chunks_local,
            &com_bI_advices,
            &r1_ext,
            deepfold_prover_param.max_mu,
        );
    let bI_r2_r4_point: Vec<F::Extension> = vec![r2_ext.clone(), r4_ext.clone()].concat();
    let (bI_chunk_polys_r2r4, bI_chunk_advices_r2r4, bI_chunk_points_r2r4) =
        expand_for_batch_open_ext::<F>(
            &bI_chunks_local,
            &com_bI_advices,
            &bI_r2_r4_point,
            deepfold_prover_param.max_mu,
        );
    let bI_r2_r5_point: Vec<F::Extension> = vec![r2_ext.clone(), r5_ext.clone()].concat();
    let (bI_chunk_polys_r2r5, bI_chunk_advices_r2r5, bI_chunk_points_r2r5) =
        expand_for_batch_open_ext::<F>(
            &bI_chunks_local,
            &com_bI_advices,
            &bI_r2_r5_point,
            deepfold_prover_param.max_mu,
        );

    // Build polys, advices vectors with expanded mat_h and bI chunks
    let mut polys: Vec<Arc<DenseMultilinearExtension<F>>> = vec![
        Arc::clone(&a_pad),
        Arc::clone(&a_pad),
        Arc::clone(&rs_a_pad),
        Arc::clone(&rs_a_pad),
    ];
    polys.extend(mat_h_chunk_polys);
    polys.push(Arc::clone(&local_mat_a_pad));
    polys.extend(bI_chunk_polys_r1);
    polys.extend(bI_chunk_polys_r2r4);
    polys.extend(bI_chunk_polys_r2r5);

    let mut advices: Vec<&DeepFoldProverCommitmentAdvice<F>> = vec![
        &com_a_advice,
        &com_a_advice,
        &com_rs_a_advice,
        &com_rs_a_advice,
    ];
    for advice in &mat_h_chunk_advices {
        advices.push(advice);
    }
    advices.push(&com_mat_a_advice_dist);
    for advice in &bI_chunk_advices_r1 {
        advices.push(advice);
    }
    for advice in &bI_chunk_advices_r2r4 {
        advices.push(advice);
    }
    for advice in &bI_chunk_advices_r2r5 {
        advices.push(advice);
    }

    // Build points vector
    let mut points_ext: Vec<Vec<F::Extension>> = vec![
        resize_point_ext::<F>(&z2_ext, deepfold_prover_param.max_mu),
        resize_point_ext::<F>(&r6_ext, deepfold_prover_param.max_mu),
        resize_point_ext::<F>(&r3_ext, deepfold_prover_param.max_mu),
        resize_point_ext::<F>(&alpha3_ext, deepfold_prover_param.max_mu),
    ];
    for point in mat_h_chunk_points {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }
    points_ext.push(resize_point_ext::<F>(&vec![r4_ext.clone(), alpha2_ext.clone()].concat(), deepfold_prover_param.max_mu));
    for point in bI_chunk_points_r1 {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }
    for point in bI_chunk_points_r2r4 {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }
    for point in bI_chunk_points_r2r5 {
        points_ext.push(resize_point_ext::<F>(&point, deepfold_prover_param.max_mu));
    }

    // Broadcast points to all parties
    let points_ext: Vec<Vec<F::Extension>> = if Net::am_master() {
        Net::recv_from_master_uniform(Some(points_ext.clone()));
        points_ext
    } else {
        Net::recv_from_master_uniform(None)
    };

    // Step 11: d_batch_open at extension field points
    let timer = start_timer!(|| "DLigesis.Open.DeepFold");

    let deepfold_batched_proof_opt = crate::deepfold::deepfold_d_batch_open_at_ext_point(
        deepfold_prover_param,
        polys,
        &advices,
        &points_ext,
        transcript,
    )?;
    end_timer!(timer);

    if Net::am_master() {
        Ok(Some(LigeSISProof {
            com_a,
            com_bI_list,
            com_rs_a,
            bI_check_proof,
            alpha2_a_bI_r2_check_proof,
            v_bI_r2_check_proof,
            rs_a_check_proof,
            mat_g_check_proofs,
            deepfold_batched_proof: deepfold_batched_proof_opt.unwrap(),
        }))
    } else {
        Ok(None)
    }
}
