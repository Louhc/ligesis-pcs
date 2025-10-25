use super::*;
use ark_std::{test_rng, UniformRand};
use ark_bls12_381::Fr as F;
use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    CryptographicSponge,
};
use ark_poly::MultilinearExtension;

pub fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    let full_rounds = 8;
    let partial_rounds = 31;
    let alpha = 17;

    let mds = vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ];

    let mut v = Vec::new();
    let mut ark_rng = test_rng();

    for _ in 0..(full_rounds + partial_rounds) {
        let mut res = Vec::new();

        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    let config = PoseidonConfig::new(full_rounds, partial_rounds, alpha, mds, v, 2, 1);
    PoseidonSponge::new(&config)
}

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
    
    let v_advice = LigeroPCS::<F>::verifier_receive_commit(&vp, &com, &mut transcript_clone).unwrap();
    let res = LigeroPCS::<F>::verify(&vp, &com, &point, &value, &v_advice, &proof, &mut transcript_clone).unwrap();
    assert!(res);
}