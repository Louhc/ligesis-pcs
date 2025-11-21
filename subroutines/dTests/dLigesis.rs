use arithmetic::math::Math;
use ark_bls12_381::Bls12_381;
use ark_ec::pairing::Pairing;
use ark_ff::{PrimeField, UniformRand};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use rand::Rng;
use std::{iter::zip, sync::Arc};
use subroutines::pcs::prelude::{LigeSISPCS, LigeSISSRS, PCSError, PolynomialCommitmentScheme};
use transcript::IOPTranscript;

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

mod common;
use common::{d_evaluate_mle, test_rng};
mod types;
use types::FGoldilocks as F;

fn test_multi<F: PrimeField>() -> Result<(), PCSError> {
    let mut rng = test_rng();
    let mu = 12;

    if Net::am_master() {
        println!(">   master: start");
        let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu)?;
        Net::recv_from_master_uniform(Some(srs.clone()));
        println!(">   master: srs distributed");

        let (pp, vp) = LigeSISPCS::<F>::setup(&srs, None, None)?;
        let poly = DenseMultilinearExtension::<F>::rand(mu, &mut rng);
        let num_party = Net::n_parties();
        let polys = (0..num_party)
            .map(|i| DenseMultilinearExtension::from_evaluations_slice(
                mu - (num_party.ilog2() as usize), 
                &poly.evaluations[i * (1 << mu) / num_party..(i + 1) * (1 << mu) / num_party]) 
            ).collect::<Vec<_>>();
        let poly_k = Arc::new(polys[0].clone());
        println!(">   master: poly distributing...");
        Net::recv_from_master(Some(polys.clone()));

        println!(">   master: poly distrbuted");

        let point = (0..mu).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
        Net::recv_from_master_uniform(Some(point.clone()));

        println!(">   master: point distributed");

        println!(">   master: start commit");
        let (com, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        println!(">   master: finish commit");

        let (com0, advice0) = LigeSISPCS::commit(&pp, &Arc::new(poly)).unwrap();

        let mut transcript = IOPTranscript::<F>::new(b"test");
        let mut transcript_clone = transcript.clone();
        println!(">   master: start open");
        let proof = LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap().unwrap();
        println!(">   master: finish open");
        let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);
        let result = LigeSISPCS::verify(&vp, &com.unwrap(), &point, &value, &proof, &mut transcript_clone)?;
        assert!(result);
    } else {
        println!(">   server({}): start", Net::party_id());
        // let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu)?;
        let srs = Net::recv_from_master_uniform::<LigeSISSRS<F>>(None);
        println!(">   server({}): srs received", Net::party_id());

        let (pp, vp) = LigeSISPCS::<F>::setup(&srs, None, None)?;
        let poly_k: DenseMultilinearExtension<F> = Net::recv_from_master(None);
        let poly_k: Arc<DenseMultilinearExtension<F>> = Arc::new(poly_k);
        println!(">   server({}): poly received", Net::party_id());
        
        let point: Vec<F> = Net::recv_from_master_uniform(None);
        println!(">   server({}): point received", Net::party_id());

        println!(">   server({}): start commit", Net::party_id());
        let (_, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        println!(">   server({}): finish commit", Net::party_id());

        let mut transcript = IOPTranscript::<F>::new(b"test");
        println!(">   server: start open");
        LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap();
        println!(">   server: finish open");
    };

    Ok(())
}

fn main() {
    common::network_run(|| {
        test_multi::<F>().unwrap();
    });
}