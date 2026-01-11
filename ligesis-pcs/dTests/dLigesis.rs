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

        // Reset stats after SRS distribution (to measure only commit/open communication)
        Net::reset_stats();

        // Setup (distributed)
        let start = Instant::now();
        let (pp, vp_opt) = LigeSISPCS::<F>::d_setup(&srs)?;
        let vp = vp_opt.unwrap();
        log_step!("Setup", start.elapsed());
        let stats_after_setup = Net::stats();
        log!("  [DEBUG] After setup: sent={:.2} MB, recv={:.2} MB",
             stats_after_setup.bytes_sent as f64 / (1024.0 * 1024.0),
             stats_after_setup.bytes_recv as f64 / (1024.0 * 1024.0));

        // Generate poly and point
        let poly_k = Arc::new(DenseMultilinearExtension::<F>::rand(mu - num_party_vars, &mut rng));
        let point: Vec<F> = (0..mu).map(|_| F::rand(&mut rng)).collect();
        Net::recv_from_master_uniform(Some(point.clone()));

        // Commit
        log!("--- Commit Phase ---");
        Net::reset_stats();
        let start = Instant::now();
        let (com, advice) = LigeSISPCS::d_commit(&pp, &poly_k).unwrap();
        log_step!("Commit", start.elapsed());
        let stats_after_commit = Net::stats();
        log!("  [DEBUG] Commit phase: sent={:.2} MB, recv={:.2} MB",
             stats_after_commit.bytes_sent as f64 / (1024.0 * 1024.0),
             stats_after_commit.bytes_recv as f64 / (1024.0 * 1024.0));

        // Open
        log!("--- Open Phase ---");
        Net::reset_stats();
        let start = Instant::now();
        let mut transcript = IOPTranscript::<F>::new(b"test");
        let proof = LigeSISPCS::d_open(&pp, &poly_k, &advice, &point, &mut transcript).unwrap().unwrap();
        log_step!("Open", start.elapsed());
        let stats_after_open = Net::stats();
        log!("  [DEBUG] Open phase: sent={:.2} MB, recv={:.2} MB",
             stats_after_open.bytes_sent as f64 / (1024.0 * 1024.0),
             stats_after_open.bytes_recv as f64 / (1024.0 * 1024.0));

        // Verify
        log!("--- Verify Phase ---");
        let start = Instant::now();
        let mut transcript = IOPTranscript::<F>::new(b"test");
        let value = LigeSISPCS::<F>::compute_value_from_proof(mu - mu / 2, &point, &proof);
        let result = LigeSISPCS::verify(&vp, &com.unwrap(), &point, &value, &proof, &mut transcript)?;
        log_step!("Verify", start.elapsed());

        // Output communication statistics (sum of all phases)
        let total_sent = stats_after_setup.bytes_sent + stats_after_commit.bytes_sent + stats_after_open.bytes_sent;
        let total_recv = stats_after_setup.bytes_recv + stats_after_commit.bytes_recv + stats_after_open.bytes_recv;
        let total_comm_bytes = total_sent + total_recv;
        let total_comm_mb = total_comm_bytes as f64 / (1024.0 * 1024.0);
        log!("========================================");
        log!("Communication Stats (master, total):");
        log!("  bytes_sent: {} ({:.2} MB)", total_sent, total_sent as f64 / (1024.0 * 1024.0));
        log!("  bytes_recv: {} ({:.2} MB)", total_recv, total_recv as f64 / (1024.0 * 1024.0));
        log!("  total:      {} ({:.2} MB)", total_comm_bytes, total_comm_mb);
        println!("COMM_TOTAL_BYTES: {}", total_comm_bytes);
        println!("COMM_TOTAL_MB: {:.2}", total_comm_mb);
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

        // Reset stats after SRS distribution (to measure only commit/open communication)
        Net::reset_stats();

        let start = Instant::now();
        let (pp, _vp) = LigeSISPCS::<F>::d_setup(&srs)?;
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

        // Output communication statistics for worker
        let stats = Net::stats();
        let total_comm_bytes = stats.bytes_sent + stats.bytes_recv;
        log!("========================================");
        log!("Communication Stats (worker {}):", party_id);
        log!("  bytes_sent: {} ({:.2} MB)", stats.bytes_sent, stats.bytes_sent as f64 / (1024.0 * 1024.0));
        log!("  bytes_recv: {} ({:.2} MB)", stats.bytes_recv, stats.bytes_recv as f64 / (1024.0 * 1024.0));
        log!("  total:      {} ({:.2} MB)", total_comm_bytes, total_comm_bytes as f64 / (1024.0 * 1024.0));
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
