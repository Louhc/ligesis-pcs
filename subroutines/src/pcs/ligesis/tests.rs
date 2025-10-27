use super::*;
use ark_std::{test_rng, UniformRand};
// use ark_bls12_381::Fr as F;
use FGoldilocks as F;
use ark_poly::MultilinearExtension;

fn inner_product<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> F {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

fn eval_mle_poly<F: PrimeField>( f: &Vec<F>, point: &Vec<F> ) -> F {
    inner_product(&f, &get_tensor(&point))
}

#[test]
fn test_ligesis_pcs() {
    let mut rng = test_rng();
    let mu = 24;

    let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();
    
    let mut transcript = IOPTranscript::<F>::new(b"ligesis_pcs_test");
    let mut transcript_clone = transcript.clone();

    let (pp, vp) = LigeSISPCS::<F>::setup(&srs, 0.into(), 0.into()).unwrap();
    let poly = Arc::new(DenseMultilinearExtension::<F>::rand(mu, &mut rng));

    let (com, advice) = LigeSISPCS::<F>::commit(&pp, &poly, &mut transcript).unwrap();
    
    let point = (0..mu).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let proof = LigeSISPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
    
    let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);
    
    let v_advice = LigeSISPCS::<F>::verifier_receive_commit(&vp, &com, &mut transcript_clone).unwrap();
    
    let res = LigeSISPCS::<F>::verify(&vp, &com, &point, &value, &v_advice, &proof, &mut transcript_clone).unwrap();
    
    assert!(res);
    assert_eq!(eval_mle_poly(&poly.evaluations, &point), value);
}
