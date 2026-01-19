use std::time::{Duration, Instant};
use ark_poly::DenseMultilinearExtension;
use ark_std::test_rng;
use std::sync::Arc;

use clap::Parser;
use ligesis_pcs::{
    ligero::LigeroPCS,
    random_field_vector_from_rng,
    PolynomialCommitmentScheme,
};
use transcript::IOPTranscript;

mod goldilocks {
    use ark_ff::fields::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "18446744069414584321"]
    #[generator = "7"]
    pub struct Config;
    pub type Fld = Fp64<MontBackend<Config, 1>>;
}
type F = goldilocks::Fld;

#[derive(Parser, Debug)]
#[command(name = "ligero_bench")]
#[command(about = "Ligero PCS Benchmark")]
struct Args {
    #[arg(short, long, default_value_t = 20)]
    mu: usize,

    #[arg(short, long, default_value_t = 1)]
    iterations: usize,

    #[arg(long, hide = true)]
    bench: bool,
}

fn fmt_duration(d: Duration) -> String {
    if d.as_secs() > 0 {
        format!("{:.3}s", d.as_secs_f64())
    } else if d.as_millis() > 0 {
        format!("{:.3}ms", d.as_secs_f64() * 1000.0)
    } else {
        format!("{:.3}us", d.as_secs_f64() * 1_000_000.0)
    }
}

fn main() {
    let args = Args::parse();
    let mu = args.mu;
    let iterations = args.iterations;

    let mut rng = test_rng();

    println!("========================================");
    println!("Ligero PCS Benchmark");
    println!("  mu = {}, iterations = {}", mu, iterations);
    println!("========================================");

    // Setup
    let start = Instant::now();
    let srs = LigeroPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();
    let (pp, vp) = LigeroPCS::<F>::setup(&srs).unwrap();
    println!("Setup:    {}", fmt_duration(start.elapsed()));

    // Prepare polynomial and point
    let evals = random_field_vector_from_rng::<F>(1 << mu, &mut rng);
    let poly = Arc::new(DenseMultilinearExtension::<F>::from_evaluations_vec(mu, evals));
    let point = random_field_vector_from_rng::<F>(mu, &mut rng);

    // Commit
    let start = Instant::now();
    let mut com = None;
    let mut advice = None;
    for _ in 0..iterations {
        let (c, a) = LigeroPCS::<F>::commit(&pp, &poly).unwrap();
        com = Some(c);
        advice = Some(a);
    }
    let commit_time = start.elapsed() / iterations as u32;
    println!("Commit:   {}", fmt_duration(commit_time));

    let com = com.unwrap();
    let advice = advice.unwrap();

    // Open
    let start = Instant::now();
    let mut proof = None;
    for _ in 0..iterations {
        let mut transcript = IOPTranscript::<F>::new(b"ligero_bench");
        let p = LigeroPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
        proof = Some(p);
    }
    let open_time = start.elapsed() / iterations as u32;
    println!("Open:     {}", fmt_duration(open_time));

    let proof = proof.unwrap();
    let log_m0 = mu / 2;
    let value = LigeroPCS::<F>::compute_value_from_proof(log_m0, &point, &proof);

    // Verify
    let start = Instant::now();
    for _ in 0..iterations {
        let mut transcript = IOPTranscript::<F>::new(b"ligero_bench");
        let res = LigeroPCS::<F>::verify(&vp, &com, &point, &value, &proof, &mut transcript).unwrap();
        assert!(res);
    }
    let verify_time = start.elapsed() / iterations as u32;
    println!("Verify:   {}", fmt_duration(verify_time));

    println!("----------------------------------------");
    println!("Total:    {}", fmt_duration(commit_time + open_time + verify_time));
    println!("========================================");
}
