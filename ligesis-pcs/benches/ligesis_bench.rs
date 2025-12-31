use std::time::Instant;
use ark_std::{test_rng, UniformRand};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use std::sync::Arc;

use ligesis_pcs::{
    ligesis::{LigeSISPCS, FGoldilocks},
    PolynomialCommitmentScheme,
};
use transcript::IOPTranscript;

type F = FGoldilocks;

// ========== 手动调整这里的参数 ==========
const MU: usize = 24;           // 多项式变量数
const ITERATIONS: usize = 1;    // 每个操作的迭代次数
// =======================================

fn main() {
    let mut rng = test_rng();

    println!("========================================");
    println!("LigeSIS PCS Benchmark");
    println!("mu = {}, iterations = {}", MU, ITERATIONS);
    println!("========================================\n");

    // Setup
    let start = Instant::now();
    let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, MU).unwrap();
    let (pp, vp) = LigeSISPCS::<F>::setup(&srs, 0.into(), 0.into()).unwrap();
    println!("Setup: {:?}", start.elapsed());

    let poly = Arc::new(DenseMultilinearExtension::<F>::rand(MU, &mut rng));
    let point: Vec<F> = (0..MU).map(|_| F::rand(&mut rng)).collect();

    // Commit
    let start = Instant::now();
    let mut com = None;
    let mut advice = None;
    for _ in 0..ITERATIONS {
        let (c, a) = LigeSISPCS::<F>::commit(&pp, &poly).unwrap();
        com = Some(c);
        advice = Some(a);
    }
    let commit_time = start.elapsed();
    println!("Commit (x{}): {:?} (avg: {:?})", ITERATIONS, commit_time, commit_time / ITERATIONS as u32);

    let com = com.unwrap();
    let advice = advice.unwrap();

    // Open
    let start = Instant::now();
    let mut proof = None;
    for _ in 0..ITERATIONS {
        let mut transcript = IOPTranscript::<F>::new(b"ligesis_pcs_bench");
        let p = LigeSISPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
        proof = Some(p);
    }
    let open_time = start.elapsed();
    println!("Open (x{}): {:?} (avg: {:?})", ITERATIONS, open_time, open_time / ITERATIONS as u32);

    let proof = proof.unwrap();
    let value = LigeSISPCS::<F>::compute_value_from_proof(MU - MU / 2, &point, &proof);

    // Verify
    let start = Instant::now();
    for _ in 0..ITERATIONS {
        let mut transcript = IOPTranscript::<F>::new(b"ligesis_pcs_bench");
        let res = LigeSISPCS::<F>::verify(&vp, &com, &point, &value, &proof, &mut transcript).unwrap();
        assert!(res);
    }
    let verify_time = start.elapsed();
    println!("Verify (x{}): {:?} (avg: {:?})", ITERATIONS, verify_time, verify_time / ITERATIONS as u32);

    println!("\n========================================");
    println!("Total (excluding setup): {:?}", commit_time + open_time + verify_time);
    println!("========================================");
}
