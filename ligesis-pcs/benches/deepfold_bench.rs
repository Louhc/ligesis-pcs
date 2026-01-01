use std::time::Instant;
use ark_bls12_381::Fr as F;
use ark_poly::DenseMultilinearExtension;
use ark_std::test_rng;
use std::sync::Arc;

use clap::Parser;
use ligesis_pcs::{
    deepfold::DeepFoldPCS,
    random_field_vector_from_rng,
    eval_mle_poly,
    PolynomialCommitmentScheme,
};
use transcript::IOPTranscript;

#[derive(Parser, Debug)]
#[command(name = "deepfold_bench")]
#[command(about = "DeepFold PCS Benchmark")]
struct Args {
    /// 多项式变量数
    #[arg(short, long, default_value_t = 20)]
    mu: usize,

    /// 每个操作的迭代次数
    #[arg(short, long, default_value_t = 1)]
    iterations: usize,

    /// 测试 batch open/verify（默认测试 single）
    #[arg(long = "test-batch")]
    test_batch: bool,

    /// batch open 的多项式数量
    #[arg(short, long, default_value_t = 3)]
    num_polys: usize,

    /// cargo bench 自动添加的参数，忽略
    #[arg(long, hide = true)]
    bench: bool,
}

fn bench_single(mu: usize, iterations: usize) {
    let mut rng = test_rng();

    println!("========================================");
    println!("DeepFold PCS Benchmark (Single)");
    println!("mu = {}, iterations = {}", mu, iterations);
    println!("========================================\n");

    // Setup
    let start = Instant::now();
    let srs = DeepFoldPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();
    let (pp, vp) = DeepFoldPCS::<F>::setup(srs).unwrap();
    println!("Setup: {:?}", start.elapsed());

    let evals = random_field_vector_from_rng::<F>(1 << mu, &mut rng);
    let poly = Arc::new(DenseMultilinearExtension::<F>::from_evaluations_vec(mu, evals));
    let point = random_field_vector_from_rng::<F>(mu, &mut rng);

    // Commit
    let start = Instant::now();
    let mut com = None;
    let mut advice = None;
    for _ in 0..iterations {
        let (c, a) = DeepFoldPCS::<F>::commit(&pp, &poly).unwrap();
        com = Some(c);
        advice = Some(a);
    }
    let commit_time = start.elapsed();
    println!(
        "Commit (x{}): {:?} (avg: {:?})",
        iterations,
        commit_time,
        commit_time / iterations as u32
    );

    let com = com.unwrap();
    let advice = advice.unwrap();

    // Open
    let start = Instant::now();
    let mut proof = None;
    for _ in 0..iterations {
        let mut transcript = IOPTranscript::<F>::new(b"deepfold_bench");
        let p = DeepFoldPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
        proof = Some(p);
    }
    let open_time = start.elapsed();
    println!(
        "Open (x{}): {:?} (avg: {:?})",
        iterations,
        open_time,
        open_time / iterations as u32
    );

    let proof = proof.unwrap();
    let value = DeepFoldPCS::<F>::compute_value_from_proof(&point, &proof);

    // Verify
    let start = Instant::now();
    for _ in 0..iterations {
        let mut transcript = IOPTranscript::<F>::new(b"deepfold_bench");
        let res =
            DeepFoldPCS::<F>::verify(&vp, &com, &point, &value, &proof, &mut transcript).unwrap();
        assert!(res);
    }
    let verify_time = start.elapsed();
    println!(
        "Verify (x{}): {:?} (avg: {:?})",
        iterations,
        verify_time,
        verify_time / iterations as u32
    );

    println!("\n========================================");
    println!("Total: {:?}", commit_time + open_time + verify_time);
    println!("========================================");
}

fn bench_batch(mu: usize, iterations: usize, num_polys: usize) {
    let mut rng = test_rng();

    println!("========================================");
    println!("DeepFold PCS Benchmark (Batch)");
    println!("mu = {}, iterations = {}, num_polys = {}", mu, iterations, num_polys);
    println!("========================================\n");

    // Setup
    let start = Instant::now();
    let srs = DeepFoldPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();
    let (pp, vp) = DeepFoldPCS::<F>::setup(srs).unwrap();
    println!("Setup: {:?}", start.elapsed());

    // Create multiple polynomials
    let polys: Vec<_> = (0..num_polys)
        .map(|_| {
            let evals = random_field_vector_from_rng::<F>(1 << mu, &mut rng);
            Arc::new(DenseMultilinearExtension::<F>::from_evaluations_vec(mu, evals))
        })
        .collect();

    let (coms, advices): (Vec<_>, Vec<_>) = polys
        .iter()
        .map(|poly| DeepFoldPCS::<F>::commit(&pp, poly).unwrap())
        .unzip();

    let points: Vec<Vec<F>> = (0..num_polys)
        .map(|_| random_field_vector_from_rng::<F>(mu, &mut rng))
        .collect();

    let evals: Vec<F> = polys
        .iter()
        .zip(points.iter())
        .map(|(poly, point)| eval_mle_poly(&poly.evaluations, point))
        .collect();

    // Batch Open
    let start = Instant::now();
    let mut batch_proof = None;
    for _ in 0..iterations {
        let mut transcript = IOPTranscript::<F>::new(b"deepfold_batch_bench");
        let advice_refs: Vec<_> = advices.iter().collect();
        let p = DeepFoldPCS::<F>::batch_open(
            &pp,
            polys.clone(),
            &advice_refs,
            &points,
            &evals,
            &mut transcript,
        )
        .unwrap();
        batch_proof = Some(p);
    }
    let batch_open_time = start.elapsed();
    println!(
        "Batch Open (x{}): {:?} (avg: {:?})",
        iterations,
        batch_open_time,
        batch_open_time / iterations as u32
    );

    let batch_proof = batch_proof.unwrap();

    // Batch Verify
    let start = Instant::now();
    for _ in 0..iterations {
        let mut transcript = IOPTranscript::<F>::new(b"deepfold_batch_bench");
        let res =
            DeepFoldPCS::<F>::batch_verify(&vp, &coms, &points, &batch_proof, &mut transcript)
                .unwrap();
        assert!(res);
    }
    let batch_verify_time = start.elapsed();
    println!(
        "Batch Verify (x{}): {:?} (avg: {:?})",
        iterations,
        batch_verify_time,
        batch_verify_time / iterations as u32
    );

    println!("\n========================================");
    println!("Total: {:?}", batch_open_time + batch_verify_time);
    println!("========================================");
}

fn main() {
    let args = Args::parse();

    if args.test_batch {
        bench_batch(args.mu, args.iterations, args.num_polys);
    } else {
        bench_single(args.mu, args.iterations);
    }
}
