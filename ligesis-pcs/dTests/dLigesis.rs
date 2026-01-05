use arithmetic::math::Math;
use ark_ff::{PrimeField, UniformRand};
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use std::sync::Arc;
use std::time::Instant;
use ligesis_pcs::{LigeSISPCS, LigeSISSRS, PCSError, PolynomialCommitmentScheme, HasQuadraticExtension};
use transcript::IOPTranscript;

use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};

mod common;
use common::{test_rng, Opt};
// Use FGoldilocks from ligesis_pcs which has HasQuadraticExtension implemented
use ligesis_pcs::FGoldilocks as F;

fn test_multi<F: PrimeField + HasQuadraticExtension>(mu: usize) -> Result<(), PCSError> {
    let mut rng = test_rng();
    let num_party = Net::n_parties();
    let num_party_vars = num_party.ilog2() as usize;
    let party_id = Net::party_id();
    let should_print = party_id == 0 || party_id == 1;
    let global_start = Instant::now();

    macro_rules! log {
        ($($arg:tt)*) => {
            if should_print {
                print!("[P{}] ", party_id);
                println!($($arg)*);
            }
        };
    }

    macro_rules! log_step {
        ($step:expr, $elapsed:expr) => {
            if should_print {
                println!("[P{}] {:12} {:>10.3?}  (@ {:.3?})", party_id, $step, $elapsed, global_start.elapsed());
            }
        };
    }

    if Net::am_master() {
        log!("========================================");
        log!("LigeSIS Distributed Test");
        log!("  mu = {}, parties = {}", mu, num_party);
        log!("========================================");

        // Gen SRS
        let start = Instant::now();
        let srs = LigeSISPCS::<F>::gen_srs_for_testing(&mut rng, mu)?;
        log_step!("Gen SRS", start.elapsed());

        // Distribute SRS
        let start = Instant::now();
        Net::recv_from_master_uniform(Some(srs.clone()));
        log_step!("Dist SRS", start.elapsed());

        // Setup
        let start = Instant::now();
        let (pp, vp) = LigeSISPCS::<F>::setup(&srs)?;
        log_step!("Setup", start.elapsed());

        // Generate poly and point
        let poly_k = Arc::new(DenseMultilinearExtension::<F>::rand(mu - num_party_vars, &mut rng));
        let point: Vec<F> = (0..mu).map(|_| F::rand(&mut rng)).collect();
        Net::recv_from_master_uniform(Some(point.clone()));

        // Commit
        log!("--- Commit Phase ---");
        let start = Instant::now();
        let (com, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        log_step!("Commit", start.elapsed());

        // Open
        log!("--- Open Phase ---");
        let start = Instant::now();
        let mut transcript = IOPTranscript::<F>::new(b"test");
        let proof = LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap().unwrap();
        log_step!("Open", start.elapsed());

        // Verify
        log!("--- Verify Phase ---");
        let start = Instant::now();
        let mut transcript = IOPTranscript::<F>::new(b"test");
        let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);
        let result = LigeSISPCS::verify(&vp, &com.unwrap(), &point, &value, &proof, &mut transcript)?;
        log_step!("Verify", start.elapsed());

        log!("========================================");
        log!("Total: {:.3?}", global_start.elapsed());
        log!("Result: {}", if result { "PASS" } else { "FAIL" });
        log!("========================================");
        assert!(result);
    } else {
        // Non-master parties
        log!("--- Setup Phase ---");

        let start = Instant::now();
        let srs = Net::recv_from_master_uniform::<LigeSISSRS<F>>(None);
        log_step!("Recv SRS", start.elapsed());

        let start = Instant::now();
        let (pp, _vp) = LigeSISPCS::<F>::setup(&srs)?;
        log_step!("Setup", start.elapsed());

        let mu = srs.mu;
        let poly_k = Arc::new(DenseMultilinearExtension::<F>::rand(mu - num_party_vars, &mut rng));
        let point: Vec<F> = Net::recv_from_master_uniform(None);

        log!("--- Commit Phase ---");
        let start = Instant::now();
        let (_, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        log_step!("Commit", start.elapsed());

        log!("--- Open Phase ---");
        let start = Instant::now();
        let mut transcript = IOPTranscript::<F>::new(b"test");
        LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap();
        log_step!("Open", start.elapsed());

        log!("========================================");
        log!("Total: {:.3?}", global_start.elapsed());
        log!("========================================");
    };

    Ok(())
}

fn main() {
    common::network_run(|opt: Opt| {
        test_multi::<F>(opt.mu).unwrap();
    });
}
