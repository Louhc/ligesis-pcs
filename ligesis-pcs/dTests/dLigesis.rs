use arithmetic::math::Math;
use ark_ff::{PrimeField, UniformRand};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use std::sync::Arc;
use ligesis_pcs::{LigeSISPCS, LigeSISSRS, PCSError, PolynomialCommitmentScheme};
use transcript::IOPTranscript;

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

mod common;
use common::test_rng;
mod types;
use types::FGoldilocks as F;

fn test_multi<F: PrimeField>() -> Result<(), PCSError> {
    let mut rng = test_rng();
    let start = std::time::Instant::now();
    
    if Net::am_master() {
        let mu = 28;
        println!(">   master: start ({} s)", start.elapsed().as_secs_f64());
        let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu)?;
        Net::recv_from_master_uniform(Some(srs.clone()));
        println!(">   master: srs distributed ({} s)", start.elapsed().as_secs_f64());
        let (pp, vp) = LigeSISPCS::<F>::setup(&srs, None, None)?;
        
        let num_party = Net::n_parties();
        let num_party_vars = num_party.ilog2() as usize;
        let poly_k = Arc::new(DenseMultilinearExtension::<F>::rand(mu - num_party_vars, &mut rng));
        

        println!(">   master: poly distrbuted ({} s)", start.elapsed().as_secs_f64());

        let point = (0..mu).map(|_| F::rand(&mut rng)).collect::<Vec<_>>();
        Net::recv_from_master_uniform(Some(point.clone()));

        println!(">   master: point distributed ({} s)", start.elapsed().as_secs_f64());

        println!(">   master: start commit ({} s)", start.elapsed().as_secs_f64());
        let (com, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        println!(">   master: finish commit ({} s)", start.elapsed().as_secs_f64());

        // println!(">   master: commit without distribution");
        // let start = std::time::Instant::now();
        // let (com0, advice0) = LigeSISPCS::commit(&pp, &Arc::new(poly)).unwrap();
        // println!(">   master: finish commit without distribution, {} s", start.elapsed().as_secs_f64());

        let mut transcript = IOPTranscript::<F>::new(b"test");
        let mut transcript_clone = transcript.clone();
        println!(">   master: start open ({} s)", start.elapsed().as_secs_f64());
        let proof = LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap().unwrap();
        println!(">   master: finish open ({} s)", start.elapsed().as_secs_f64());
        let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);
        let result = LigeSISPCS::verify(&vp, &com.unwrap(), &point, &value, &proof, &mut transcript_clone)?;
        println!(">   master: done ({} s)", start.elapsed().as_secs_f64());
        println!("result = {}", result);
        assert!(result);
    } else {
        println!(">   server({}): start ({} s)", Net::party_id(), start.elapsed().as_secs_f64());
        // let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu)?;
        let srs = Net::recv_from_master_uniform::<LigeSISSRS<F>>(None);
        println!(">   server({}): srs received ({} s)", Net::party_id(), start.elapsed().as_secs_f64());

        let (pp, vp) = LigeSISPCS::<F>::setup(&srs, None, None)?;
        let mu = srs.mu;
        let num_party = Net::n_parties();
        let num_party_vars = num_party.ilog2() as usize;

        let poly_k = Arc::new(DenseMultilinearExtension::<F>::rand(mu - num_party_vars, &mut rng));
        // let poly_k: DenseMultilinearExtension<F> = Net::recv_from_master(None);
        // let poly_k: Arc<DenseMultilinearExtension<F>> = Arc::new(poly_k);
        println!(">   server({}): poly received ({} s)", Net::party_id(), start.elapsed().as_secs_f64());
        
        let point: Vec<F> = Net::recv_from_master_uniform(None);
        println!(">   server({}): point received ({} s)", Net::party_id(), start.elapsed().as_secs_f64());

        println!(">   server({}): start commit ({} s)", Net::party_id(), start.elapsed().as_secs_f64());
        let (_, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        println!(">   server({}): finish commit ({} s)", Net::party_id(), start.elapsed().as_secs_f64());

        let mut transcript = IOPTranscript::<F>::new(b"test");
        println!(">   server: start open ({} s)", start.elapsed().as_secs_f64());
        LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap();
        println!(">   server: finish open ({} s)", start.elapsed().as_secs_f64());
    };

    Ok(())
}

fn main() {
    common::network_run(|| {
        test_multi::<F>().unwrap();
    });
}
