use super::*;
use ark_std::{test_rng, UniformRand};
// use ark_bls12_381::Fr as F;
use ark_poly::MultilinearExtension;
use FGoldilocks as F;

extern crate test;
use test::Bencher;

fn inner_product<F: PrimeField>(a: &Vec<F>, b: &Vec<F>) -> F {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

fn eval_mle_poly<F: PrimeField>(f: &Vec<F>, point: &Vec<F>) -> F {
    inner_product(&f, &get_tensor(&point))
}

#[test]
fn test_ligesis_pcs() {
    let mut rng = test_rng();
    let mu = 18;

    let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();

    let mut transcript = IOPTranscript::<F>::new(b"ligesis_pcs_test");
    let mut transcript_clone = transcript.clone();

    let (pp, vp) = LigeSISPCS::<F>::setup(&srs, 0.into(), 0.into()).unwrap();
    let poly = Arc::new(DenseMultilinearExtension::<F>::rand(mu, &mut rng));

    let start = std::time::Instant::now();
    let (com, advice) = LigeSISPCS::<F>::commit(&pp, &poly).unwrap();
    println!("Commit: {} s", start.elapsed().as_secs_f64());

    let point = (0..mu).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();

    let start = std::time::Instant::now();
    let proof = LigeSISPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
    println!("Open: {} s", start.elapsed().as_secs_f64());

    let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);

    let start = std::time::Instant::now();
    let res =
        LigeSISPCS::<F>::verify(&vp, &com, &point, &value, &proof, &mut transcript_clone).unwrap();
    println!("Verify: {} s", start.elapsed().as_secs_f64());
    println!("Proof Size: {} B", proof.size_in_bytes());

    assert!(res);
    assert_eq!(eval_mle_poly(&poly.evaluations, &point), value);
    assert!(false);
}

#[bench]
fn bench_ligesis_pcs(b: &mut Bencher) {
    let mut rng = test_rng();
    let mu = 18;

    let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();

    let mut transcript = IOPTranscript::<F>::new(b"ligesis_pcs_test");
    let mut transcript_clone = transcript.clone();

    let (pp, vp) = LigeSISPCS::<F>::setup(&srs, 0.into(), 0.into()).unwrap();
    let poly = Arc::new(DenseMultilinearExtension::<F>::rand(mu, &mut rng));

    let start = std::time::Instant::now();
    let (com, advice) = LigeSISPCS::<F>::commit(&pp, &poly).unwrap();
    println!("Commit: {} s", start.elapsed().as_secs_f64());
    // b.iter(|| LigeSISPCS::<F>::commit(&pp, &poly));

    let point = (0..mu).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();

    let start = std::time::Instant::now();
    let proof = LigeSISPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
    println!("Open: {} s", start.elapsed().as_secs_f64());

    let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);

    let start = std::time::Instant::now();
    let res =
        LigeSISPCS::<F>::verify(&vp, &com, &point, &value, &proof, &mut transcript_clone).unwrap();
    println!("Verify: {} s", start.elapsed().as_secs_f64());
    println!("Proof Size: {} B", proof.size_in_bytes());

    assert!(res);
    assert_eq!(eval_mle_poly(&poly.evaluations, &point), value);
    assert!(res);
}
