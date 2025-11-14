use super::*;
use ark_bls12_381::Fr as F;
use ark_poly::MultilinearExtension;
use ark_std::{test_rng, UniformRand};

#[test]
fn test_ligero_pcs() {
    let mut rng = test_rng();
    let srs = LigeroPCS::<F>::gen_srs_for_testing(&mut rng, 18).unwrap();
    let mut transcript = IOPTranscript::<F>::new(b"ligero_pcs_test");
    let mut transcript_clone = transcript.clone();
    let (pp, vp) = LigeroPCS::<F>::setup(&srs, 0.into(), 0.into()).unwrap();
    let poly = Arc::new(DenseMultilinearExtension::<F>::rand(18, &mut rng));

    let (com, advice) = LigeroPCS::<F>::commit(&pp, &poly, &mut transcript).unwrap();
    let point = (0..18).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
    let proof = LigeroPCS::<F>::open(&pp, &poly, &advice, &point, &mut transcript).unwrap();
    let value = LigeroPCS::<F>::compute_value_from_proof(pp.1, &point, &proof);

    let v_advice =
        LigeroPCS::<F>::verifier_receive_commit(&vp, &com, &mut transcript_clone).unwrap();
    let res = LigeroPCS::<F>::verify(
        &vp,
        &com,
        &point,
        &value,
        &v_advice,
        &proof,
        &mut transcript_clone,
    )
    .unwrap();
    assert!(res);
}
