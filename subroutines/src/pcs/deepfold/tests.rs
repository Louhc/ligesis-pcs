use crate::rand::random_field_vector_from_rng;

use super::*;
use ark_std::test_rng;
use ark_bls12_381::Fr as F;

#[test]
fn test_deepfold_pcs(){
    let mut rng = test_rng();
    let mu = 10;
    let srs = DeepFoldPCS::<F>::gen_srs_for_testing(&mut rng, mu).unwrap();
    let (pp, vp) = DeepFoldPCS::<F>::trim(srs, Some(mu), Some(mu)).unwrap();
    let poly = random_field_vector_from_rng::<F>(1 << mu, &mut rng);
    let poly = DenseMultilinearExtension::<F>::from_evaluations_vec(mu, poly);
    let poly_arc = Arc::new(poly);
    let mut transcript = IOPTranscript::new(b"test");
    let mut transcript_clone = transcript.clone();
    
    let (com, advice) = DeepFoldPCS::<F>::commit(&pp, &poly_arc, &mut transcript).unwrap();
    
    let point = random_field_vector_from_rng::<F>(mu, &mut rng);
    let proof = DeepFoldPCS::<F>::open(&pp, &poly_arc, &advice, &point, &mut transcript).unwrap();
    
    let value = DeepFoldPCS::compute_value_from_proof(&point, &proof);
    let result = DeepFoldPCS::<F>::verify(&vp, &com, &point, &value, &proof, &mut transcript_clone).unwrap();

    assert!(result);
    assert_eq!(eval_mle_poly(&poly_arc.evaluations, &point), value);
}