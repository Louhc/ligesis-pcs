use crate::{
    errors::PCSError, hash::*, rscode::*, utils::*,
    IOPProof, PolyIOP, PolynomialCommitmentScheme,
    sumcheck::SumCheck,
};
use arithmetic::{VPAuxInfo, VirtualPolynomial};
use ark_ff::PrimeField;
use ark_poly::{DenseMultilinearExtension, EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{
    borrow::Borrow, end_timer, marker::PhantomData, rand::Rng, start_timer, sync::Arc, vec,
    vec::Vec,
};
use deNetwork::{DeMultiNet as Net, DeNet, DeSerNet};
use transcript::IOPTranscript;

#[cfg(test)]
mod tests;
mod utils;
use utils::*;

/// DeepFold Polynomial Commitment Scheme
pub struct DeepFoldPCS<F: PrimeField> {
    #[doc(hidden)]
    phantom: PhantomData<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, Copy)]
pub struct DeepFoldSRS<F: PrimeField> {
    pub max_mu: usize,
    pub l0: GeneralEvaluationDomain<F>,
    pub s: usize,
}

impl<F: PrimeField> Default for DeepFoldSRS<F> {
    fn default() -> Self {
        DeepFoldSRS {
            max_mu: 0,
            l0: GeneralEvaluationDomain::<F>::new(1).unwrap(),
            s: 0,
        }
    }
}

#[derive(Clone)]
pub struct DeepFoldProverParam<F: PrimeField> {
    pub max_mu: usize,
    pub l0: GeneralEvaluationDomain<F>,
    pub s: usize,
}

#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct DeepFoldVerifierParam<F: PrimeField> {
    pub max_mu: usize,
    pub len_l0: usize,
    pub g: F,
    pub s: usize,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
/// proof of opening
pub struct DeepFoldProof<F: PrimeField> {
    pub linear_polys: Vec<Vec<(F, F)>>,
    pub mt_roots: Vec<Byte32>,
    pub f_mu: F,
    pub mt_proofs: Vec<Vec<(usize, (F, F), Vec<F>, Vec<Byte32>)>>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
pub struct DeepFoldBatchedProof<F: PrimeField> {
    pub deepfold_proof: DeepFoldProof<F>,
    pub sum_check_proof: IOPProof<F>,
    pub mt_proofs_for_mt0: Vec<Vec<(Vec<F>, Vec<Byte32>)>>,
    pub evals: Vec<F>,
    pub sum_check_evals: Vec<F>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct DeepFoldProverCommitmentAdvice<F: PrimeField> {
    pub f0: Vec<F>,
    pub mt0: MerkleTree,
    pub v0: Vec<F>,
    /// Full polynomial evaluations for distributed setting (only master has this)
    pub f_tilde: Vec<F>,
    /// Upper tree for distributed setting (only master has this)
    pub upper_tree: Option<MerkleTree>,
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq, Default)]
pub struct DeepFoldCommitment {
    pub mu: usize,
    pub rt0: Byte32,
}

impl<F: PrimeField> DeepFoldPCS<F> {
    pub fn compute_value_from_proof(point: &Vec<F>, proof: &DeepFoldProof<F>) -> F {
        eval_linear_poly(&proof.linear_polys[0][0], &point[0])
    }

    /// Compute the claimed value from a distributed proof
    /// In the distributed setting, the proof structure is the same, so we just use the first linear polynomial
    pub fn compute_value_from_proof_distributed(point: &Vec<F>, proof: &DeepFoldProof<F>, _num_party: usize) -> F {
        eval_linear_poly(&proof.linear_polys[0][0], &point[0])
    }
}

impl<F: PrimeField> PolynomialCommitmentScheme<F> for DeepFoldPCS<F> {
    // Parameters
    type ProverParam = DeepFoldProverParam<F>;
    type VerifierParam = DeepFoldVerifierParam<F>;
    type SRS = DeepFoldSRS<F>;
    // Polynomial and its associated types
    type Polynomial = Arc<DenseMultilinearExtension<F>>;
    type ProverCommitmentAdvice = DeepFoldProverCommitmentAdvice<F>;
    type Point = Vec<F>;
    type Evaluation = F;
    // Commitments and proofs
    type Commitment = DeepFoldCommitment; // merkle tree root
    type Proof = DeepFoldProof<F>; // merkle tree paths, columes of `E`
    type BatchProof = DeepFoldBatchedProof<F>; //

    fn gen_srs_for_testing<R: Rng>(_rng: &mut R, log_size: usize) -> Result<Self::SRS, PCSError> {
        let max_mu = log_size;
        let len_l0 = (1 << max_mu) * 8;
        let l0 = GeneralEvaluationDomain::<F>::new(len_l0).unwrap();
        let s = 33;
        Ok(DeepFoldSRS { max_mu, l0, s })
    }

    fn setup(
        srs: impl Borrow<Self::SRS>,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), PCSError> {
        let srs = srs.borrow();
        Ok((
            DeepFoldProverParam {
                max_mu: srs.max_mu,
                l0: srs.l0,
                s: srs.s,
            },
            DeepFoldVerifierParam {
                max_mu: srs.max_mu,
                len_l0: srs.l0.size(),
                g: srs.l0.element(1),
                s: srs.s,
            },
        ))
    }

    fn commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<(Self::Commitment, Self::ProverCommitmentAdvice), PCSError> {
        let &Self::ProverParam { max_mu, l0, s } = prover_param.borrow();
        let mu = poly.num_vars;
        assert!(mu <= max_mu);

        let f0 = evals_to_coeffs(mu, &poly.evaluations);
        let v0 = l0.fft(&f0);

        // let mt0 = MerkleTree::new(&v0.iter().map(|&x|
        // compute_sha256_row(&[x])).collect());
        let mt0 = build_merkle_tree(&v0);

        let rt0 = mt0.root();
        Ok((
            DeepFoldCommitment { mu, rt0 },
            DeepFoldProverCommitmentAdvice { f0, mt0, v0, f_tilde: poly.evaluations.clone(), upper_tree: None },
        ))
    }

    /// Distributed commit: each party has local polynomial evaluations
    /// Each party builds local subtree, master builds upper tree from collected roots
    /// Returns (Option<Commitment>, Advice) - commitment is Some only for master
    fn d_commit(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
    ) -> Result<(Option<Self::Commitment>, Self::ProverCommitmentAdvice), PCSError> {
        let &Self::ProverParam { max_mu, l0, s } = prover_param.borrow();
        let num_party = Net::n_parties();
        let num_party_vars = num_party.ilog2() as usize;

        // Each party has local evaluations of size 2^local_mu
        let local_mu = poly.num_vars;
        let mu = local_mu + num_party_vars;
        assert!(mu <= max_mu);

        // Step 1: Gather all evaluations to master
        let timer = start_timer!(|| "DCommit.GatherEvals");
        let all_evals_opt = Net::send_to_master(&poly.evaluations);
        end_timer!(timer);

        // Step 2: Master computes full f0, v0, and distributes leaf hashes
        let timer = start_timer!(|| "DCommit.FFT");
        let (f0, v0, f_tilde, local_leaves, leaf_size): (Vec<F>, Vec<F>, Vec<F>, Vec<Byte32>, usize) = if Net::am_master() {
            let all_evals: Vec<Vec<F>> = all_evals_opt.unwrap();
            let full_evals: Vec<F> = all_evals.into_iter().flatten().collect();

            // Compute full coefficients and FFT
            let f0 = evals_to_coeffs(mu, &full_evals);
            let v0 = l0.fft(&f0);

            // Compute leaf hashes from v0
            let (all_leaves, leaf_size) = compute_leaf_hashes(&v0);

            // Split leaf hashes into chunks for each party
            let chunk_size = all_leaves.len() / num_party;
            let leaf_chunks: Vec<Vec<Byte32>> = (0..num_party)
                .map(|i| all_leaves[i * chunk_size..(i + 1) * chunk_size].to_vec())
                .collect();

            // Distribute leaf hash chunks to each party
            let local_leaves = Net::recv_from_master(Some(leaf_chunks));

            // Also broadcast leaf_size
            Net::recv_from_master_uniform(Some(leaf_size));

            (f0, v0, full_evals, local_leaves, leaf_size)
        } else {
            // Workers receive their portion of leaf hashes
            let local_leaves: Vec<Byte32> = Net::recv_from_master(None);
            let leaf_size: usize = Net::recv_from_master_uniform(None);
            (vec![], vec![], vec![], local_leaves, leaf_size)
        };
        end_timer!(timer);

        // Step 3: Each party builds local subtree
        let timer = start_timer!(|| "DCommit.DMerkle");
        let local_mt0 = MerkleTree::with_leaf_size(&local_leaves, leaf_size);
        let local_root = local_mt0.root();

        // Gather all local roots to master to build upper tree
        let all_roots_opt = Net::send_to_master(&local_root);
        end_timer!(timer);

        if Net::am_master() {
            let all_roots: Vec<Byte32> = all_roots_opt.unwrap();

            // Build upper tree from all party roots
            let upper_tree = MerkleTree::with_leaf_size(&all_roots, leaf_size);
            let rt0 = upper_tree.root();

            Ok((
                Some(DeepFoldCommitment { mu, rt0 }),
                DeepFoldProverCommitmentAdvice { f0, mt0: local_mt0, v0, f_tilde, upper_tree: Some(upper_tree) },
            ))
        } else {
            Ok((
                None,
                DeepFoldProverCommitmentAdvice {
                    f0: vec![],
                    mt0: local_mt0,
                    v0: vec![],
                    f_tilde: vec![],
                    upper_tree: None,
                },
            ))
        }
    }

    fn open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::Proof, PCSError> {
        let &Self::ProverParam { max_mu, l0, s } = prover_param.borrow();

        let mu = poly.num_vars;

        assert!(mu <= max_mu);

        let Self::ProverCommitmentAdvice { f0, mt0, v0, f_tilde: _, upper_tree: _ } = advice.clone();
        let mut a = vec![Vec::new()];
        let mut f_tilde = vec![poly.evaluations.clone()];
        let mut f = vec![f0];
        let mut alpha = vec![F::ZERO];
        let mut linear_polys = Vec::new();
        let mut l = vec![l0];
        l.append(
            &mut (1..mu + 1)
                .map(|i| GeneralEvaluationDomain::<F>::new(l0.size() >> i).unwrap())
                .collect::<Vec<_>>(),
        );
        let mut v = vec![v0];
        let mut mt_roots = vec![mt0.root().clone()];
        let mut mt = vec![mt0];
        let mut mt_proofs = Vec::new();
        let mut f_mu = F::ZERO;
        let mut r = vec![F::ZERO];

        // Step 1
        a[0].push(point.clone());

        // Step 2
        for i in 1..mu + 1 {
            // Step 2.a
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);
            a[i - 1].push(get_alpha_powers::<F>(alpha[i], mu - i + 1));
            let (f0, f1) = split_even_odd(&f_tilde[i - 1]);
            let (fe, fo) = split_even_odd(&f[i - 1]);
            // Step 2.b
            if i == mu {
                linear_polys.push(vec![(f_tilde[i - 1][0], f_tilde[i - 1][1])]);
            } else {
                linear_polys.push(
                    a[i - 1]
                        .iter()
                        .map(|w| {
                            assert!(!w.is_empty());
                            let w_tensor = get_tensor(&w[1..].to_vec());
                            (inner_product(&w_tensor, &f0), inner_product(&w_tensor, &f1))
                        })
                        .collect::<Vec<_>>(),
                );

                a.push(a[i - 1].iter().map(|w| w[1..].to_vec()).collect::<Vec<_>>());
            }
            // Step 2.c
            let ri = transcript.get_and_append_challenge(b"r")?;
            r.push(ri);
            // Step 2.d
            f.push(vector_add(&fe, &scalar_vector_product(ri, &fo)));
            f_tilde.push(vector_add(
                &scalar_vector_product(F::ONE - ri, &f0),
                &scalar_vector_product(ri, &f1),
            ));
            // Step 2.e
            v.push(l[i].fft(&f[i]));
            if i == mu {
                f_mu = v[i][0];
            } else {
                let mti = build_merkle_tree(&v[i]);
                mt_roots.push(mti.root().clone());
                mt.push(mti);
            }
        }
        // Step 4
        for t in 0..s {
            // Step 4.a
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, l[0].size())?[0];
            // Step 4.b
            mt_proofs.push(Vec::new());
            for i in 0..mu {
                mt_proofs[t].push(open_merkle_tree_at_conjugate_points(&mt[i], &v[i], beta));
                if beta >= l[i + 1].size() {
                    beta -= l[i + 1].size();
                }
            }
        }
        Ok(DeepFoldProof {
            linear_polys,
            mt_roots,
            f_mu,
            mt_proofs,
        })
    }

    fn d_open(
        prover_param: impl Borrow<Self::ProverParam>,
        poly: &Self::Polynomial,
        advice: &Self::ProverCommitmentAdvice,
        point: &Self::Point,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Option<Self::Proof>, PCSError> {
        let &Self::ProverParam { max_mu, l0, s } = prover_param.borrow();
        let num_party = Net::n_parties();
        let num_party_vars = num_party.ilog2() as usize;

        // Each party has local evaluations of size 2^local_mu
        let local_mu = poly.num_vars;
        let mu = local_mu + num_party_vars;
        assert!(mu <= max_mu);

        // Initialize structures - use f0, v0, mt0, f_tilde, upper_tree from advice (computed in d_commit)
        let Self::ProverCommitmentAdvice { f0, mt0, v0, f_tilde: advice_f_tilde, upper_tree } = advice.clone();
        let mut l = vec![l0];
        l.append(
            &mut (1..mu + 1)
                .map(|i| GeneralEvaluationDomain::<F>::new(l0.size() >> i).unwrap())
                .collect::<Vec<_>>(),
        );

        let timer = start_timer!(|| "DOpen.ComputeProof");

        // Each party stores local subtrees, master also stores upper trees
        let mut local_mts: Vec<MerkleTree> = vec![mt0.clone()];
        let mut upper_mts: Vec<Option<MerkleTree>> = vec![upper_tree.clone()];
        let mut is_distributed: Vec<bool> = vec![true]; // Track which rounds use distributed Merkle
        let mut mt_roots: Vec<Byte32> = vec![];

        // Get mt_roots[0] from upper_tree (master) or placeholder (workers)
        if Net::am_master() {
            mt_roots.push(upper_tree.as_ref().unwrap().root());
        } else {
            mt_roots.push(Byte32::default());
        }

        // Master-only data for polynomial computation
        let mut a = vec![Vec::new()];
        let mut f_tilde: Vec<Vec<F>> = if Net::am_master() {
            vec![advice_f_tilde]  // Use f_tilde from advice (computed in d_commit)
        } else {
            Vec::new()
        };
        let mut f = vec![f0.clone()];
        let mut alpha = vec![F::ZERO];
        let mut linear_polys = Vec::new();
        let mut v = vec![v0.clone()];
        let mut f_mu = F::ZERO;
        let mut r = vec![F::ZERO];

        // Step 1
        a[0].push(point.clone());

        // Step 2: Main loop - all parties participate in building distributed Merkle trees
        for i in 1..mu + 1 {
            // Step 2.a: Get alpha challenge - master generates and broadcasts to workers
            let alpha_i = if Net::am_master() {
                let a = transcript.get_and_append_challenge(b"alpha")?;
                Net::recv_from_master_uniform(Some(a));
                a
            } else {
                Net::recv_from_master_uniform(None)
            };
            alpha.push(alpha_i);

            if Net::am_master() {
                a[i - 1].push(get_alpha_powers::<F>(alpha[i], mu - i + 1));
                let (f0_split, f1) = split_even_odd(&f_tilde[i - 1]);
                let (fe, fo) = split_even_odd(&f[i - 1]);

                // Step 2.b: Compute linear_polys (master only)
                if i == mu {
                    linear_polys.push(vec![(f_tilde[i - 1][0], f_tilde[i - 1][1])]);
                } else {
                    linear_polys.push(
                        a[i - 1]
                            .iter()
                            .map(|w| {
                                assert!(!w.is_empty());
                                let w_tensor = get_tensor(&w[1..].to_vec());
                                (inner_product(&w_tensor, &f0_split), inner_product(&w_tensor, &f1))
                            })
                            .collect::<Vec<_>>(),
                    );
                    a.push(a[i - 1].iter().map(|w| w[1..].to_vec()).collect::<Vec<_>>());
                }

                // Step 2.c: Get r challenge - master generates and broadcasts to workers
                let ri = transcript.get_and_append_challenge(b"r")?;
                Net::recv_from_master_uniform(Some(ri));

                // Step 2.d: Compute f[i] and f_tilde[i]
                r.push(ri);
                f.push(vector_add(&fe, &scalar_vector_product(ri, &fo)));
                f_tilde.push(vector_add(
                    &scalar_vector_product(F::ONE - ri, &f0_split),
                    &scalar_vector_product(ri, &f1),
                ));

                // Step 2.e: Compute v[i] = FFT(f[i])
                let vi = l[i].fft(&f[i]);
                v.push(vi.clone());

                if i == mu {
                    f_mu = v[i][0];
                } else {
                    // Check if we can use distributed Merkle tree
                    let (all_leaves, leaf_size) = compute_leaf_hashes(&vi);
                    let can_distribute = all_leaves.len() >= num_party;

                    // First broadcast can_distribute flag so workers know what to expect
                    Net::recv_from_master_uniform(Some(can_distribute));

                    if can_distribute {
                        // Distribute leaf hashes for distributed Merkle tree
                        let chunk_size = all_leaves.len() / num_party;
                        let leaf_chunks: Vec<Vec<Byte32>> = (0..num_party)
                            .map(|j| all_leaves[j * chunk_size..(j + 1) * chunk_size].to_vec())
                            .collect();

                        let local_leaves: Vec<Byte32> = Net::recv_from_master(Some(leaf_chunks));
                        Net::recv_from_master_uniform(Some(leaf_size));

                        // Build local Merkle tree
                        let local_mt = MerkleTree::with_leaf_size(&local_leaves, leaf_size);
                        let local_root = local_mt.root();
                        local_mts.push(local_mt);

                        // Gather local roots to build upper tree
                        let all_roots: Vec<Byte32> = Net::send_to_master(&local_root).unwrap();
                        let upper_tree = MerkleTree::with_leaf_size(&all_roots, leaf_size);
                        mt_roots.push(upper_tree.root());
                        upper_mts.push(Some(upper_tree));
                        is_distributed.push(true);
                    } else {
                        // Too few leaves for distribution - master builds full tree alone
                        let full_mt = MerkleTree::with_leaf_size(&all_leaves, leaf_size);
                        mt_roots.push(full_mt.root());
                        local_mts.push(full_mt);
                        upper_mts.push(None);
                        is_distributed.push(false);
                    }
                }
            } else {
                // Workers: receive r challenge from master
                let ri: F = Net::recv_from_master_uniform(None);
                r.push(ri);

                // Workers: participate in distributed Merkle tree construction
                if i != mu {
                    // First receive can_distribute flag
                    let can_distribute: bool = Net::recv_from_master_uniform(None);

                    if can_distribute {
                        // Receive leaf hashes for distributed Merkle tree
                        let local_leaves: Vec<Byte32> = Net::recv_from_master(None);
                        let leaf_size: usize = Net::recv_from_master_uniform(None);

                        // Build local Merkle tree
                        let local_mt = MerkleTree::with_leaf_size(&local_leaves, leaf_size);
                        let local_root = local_mt.root();
                        local_mts.push(local_mt);

                        // Send local root to master
                        Net::send_to_master(&local_root);

                        // Workers don't have upper trees
                        upper_mts.push(None);
                        is_distributed.push(true);
                    } else {
                        // Non-distributed mode: workers just push placeholders
                        local_mts.push(MerkleTree::default());
                        upper_mts.push(None);
                        is_distributed.push(false);
                    }
                }
            }
        }
        end_timer!(timer);

        // Step 4: Generate merkle proofs
        let timer = start_timer!(|| "DOpen.GenProofs");
        let mut mt_proofs = Vec::new();
        for t in 0..s {
            // All parties need to participate in transcript to get beta
            let mut beta = if Net::am_master() {
                let b = transcript.get_and_append_challenge_indices(b"beta", 1, l[0].size())?[0];
                Net::recv_from_master_uniform(Some(b));
                b
            } else {
                Net::recv_from_master_uniform(None)
            };

            let mut proofs_for_t = Vec::new();
            for i in 0..mu {
                let vi_len = l[i].size();
                let leaf_size = local_mts[i].leaf_size();
                let step = vi_len / leaf_size;
                let local_beta = beta % step;

                if is_distributed[i] {
                    // Use d_prove for distributed proof generation
                    let proof_opt = MerkleTree::d_prove(local_beta, &local_mts[i], upper_mts[i].as_ref());

                    if Net::am_master() {
                        let merkle_proof = proof_opt.unwrap();
                        let beta_prime = if beta >= vi_len / 2 {
                            beta - vi_len / 2
                        } else {
                            beta + vi_len / 2
                        };
                        let leaf_elements = get_leaf_elements(&v[i], local_beta, step, leaf_size);
                        proofs_for_t.push((beta, (v[i][beta], v[i][beta_prime]), leaf_elements, merkle_proof));
                    }
                } else if Net::am_master() {
                    // Non-distributed: master uses regular prove
                    let merkle_proof = local_mts[i].prove(local_beta);
                    let beta_prime = if beta >= vi_len / 2 {
                        beta - vi_len / 2
                    } else {
                        beta + vi_len / 2
                    };
                    let leaf_elements = get_leaf_elements(&v[i], local_beta, step, leaf_size);
                    proofs_for_t.push((beta, (v[i][beta], v[i][beta_prime]), leaf_elements, merkle_proof));
                }

                if beta >= l[i + 1].size() {
                    beta -= l[i + 1].size();
                }
            }
            if Net::am_master() {
                mt_proofs.push(proofs_for_t);
            }
        }
        end_timer!(timer);

        if Net::am_master() {
            Ok(Some(DeepFoldProof {
                linear_polys,
                mt_roots,
                f_mu,
                mt_proofs,
            }))
        } else {
            Ok(None)
        }
    }

    fn batch_open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomials: Vec<Self::Polynomial>,
        advices: &[&Self::ProverCommitmentAdvice],
        points: &[Self::Point],
        _evals: &[Self::Evaluation],
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Self::BatchProof, PCSError> {
        let &Self::ProverParam { max_mu, l0, s } = prover_param.borrow();
        let num_poly = polynomials.len();
        let mu = max_mu;
        assert!(polynomials.iter().all(|poly| poly.num_vars == mu));
        assert!(points.iter().all(|point| point.len() == mu));
        assert!(points.len() == num_poly && advices.len() == num_poly);
        let mt0_list = advices.iter().map(|advice| &advice.mt0).collect::<Vec<_>>();

        // SumCheck Phase
        let timer = start_timer!(|| "DeepFold.Sumcheck");
        let r = transcript.get_and_append_challenge(b"batched_sumcheck")?;
        let mut sum_check = VirtualPolynomial::new(max_mu);
        for i in 0..num_poly {
            sum_check
                .add_mle_list(
                    [
                        evals_to_arcpoly(&polynomials[i].evaluations),
                        evals_to_arcpoly(&get_tensor(&points[i])),
                    ],
                    r.pow([i as u64]),
                )
                .map_err(|e| PCSError::VirtualPolynomialError(format!("{:?}", e)))?;
        }
        let sum_check_proof = <PolyIOP<F> as SumCheck<F>>::prove(sum_check, transcript).map_err(|e| PCSError::SumCheckError(format!("{:?}", e)))?;
        let point = sum_check_proof.point.clone();
        let sum_check_evals = polynomials
            .iter()
            .map(|poly| eval_mle_poly(&poly.evaluations, &point))
            .collect::<Vec<_>>();
        end_timer!(timer);

        // Batched Open Phase - Compute combined polynomial WITHOUT building mt0
        let timer = start_timer!(|| "DeepFold.BatchedOpen");
        let gamma = transcript.get_and_append_challenge_vectors(b"gamma", num_poly)?;
        let poly_evals: Vec<F> = (0..1 << max_mu)
            .map(|i| {
                (0..num_poly)
                    .map(|j| gamma[j] * polynomials[j].evaluations[i])
                    .sum::<F>()
            })
            .collect();

        // Compute f0 and v0 for combined polynomial (needed for subsequent rounds)
        let f0 = evals_to_coeffs(mu, &poly_evals);
        let v0 = l0.fft(&f0);

        // Initialize structures - NO mt0 needed
        let mut a = vec![Vec::new()];
        let mut f_tilde = vec![poly_evals];
        let mut f = vec![f0];
        let mut alpha = vec![F::ZERO];
        let mut linear_polys = Vec::new();
        let mut l = vec![l0];
        l.append(
            &mut (1..mu + 1)
                .map(|i| GeneralEvaluationDomain::<F>::new(l0.size() >> i).unwrap())
                .collect::<Vec<_>>(),
        );
        let mut v = vec![v0];
        let mut mt_roots = Vec::new(); // Will be filled starting from round 1
        let mut mt = Vec::new();
        let mut f_mu = F::ZERO;
        let mut r_vals = vec![F::ZERO];

        // Step 1
        a[0].push(point.clone());

        // Step 2
        for i in 1..mu + 1 {
            // Step 2.a
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);
            a[i - 1].push(get_alpha_powers::<F>(alpha[i], mu - i + 1));
            let (f0_split, f1) = split_even_odd(&f_tilde[i - 1]);
            let (fe, fo) = split_even_odd(&f[i - 1]);
            // Step 2.b
            if i == mu {
                linear_polys.push(vec![(f_tilde[i - 1][0], f_tilde[i - 1][1])]);
            } else {
                linear_polys.push(
                    a[i - 1]
                        .iter()
                        .map(|w| {
                            assert!(!w.is_empty());
                            let w_tensor = get_tensor(&w[1..].to_vec());
                            (inner_product(&w_tensor, &f0_split), inner_product(&w_tensor, &f1))
                        })
                        .collect::<Vec<_>>(),
                );

                a.push(a[i - 1].iter().map(|w| w[1..].to_vec()).collect::<Vec<_>>());
            }
            // Step 2.c
            let ri = transcript.get_and_append_challenge(b"r")?;
            r_vals.push(ri);
            // Step 2.d
            f.push(vector_add(&fe, &scalar_vector_product(ri, &fo)));
            f_tilde.push(vector_add(
                &scalar_vector_product(F::ONE - ri, &f0_split),
                &scalar_vector_product(ri, &f1),
            ));
            // Step 2.e
            v.push(l[i].fft(&f[i]));
            if i == mu {
                f_mu = v[i][0];
            } else {
                // Build merkle trees starting from i=1 (skip mt0)
                let mti = build_merkle_tree(&v[i]);
                mt_roots.push(mti.root().clone());
                mt.push(mti);
            }
        }

        // Step 4: Generate merkle proofs (starting from index 1)
        let mut mt_proofs = Vec::new();
        for t in 0..s {
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, l[0].size())?[0];
            mt_proofs.push(Vec::new());

            // For i=0, just store the values without merkle proof (will be verified via linear combination)
            let leaf_size = LEAF_SIZE.min(v[0].len());
            let step = v[0].len() / leaf_size;
            let local_beta = beta % step;
            let beta_prime = if beta >= v[0].len() / 2 {
                beta - v[0].len() / 2
            } else {
                beta + v[0].len() / 2
            };
            mt_proofs[t].push((
                beta, // Store original position
                (v[0][beta], v[0][beta_prime]),
                get_leaf_elements(&v[0], local_beta, step, leaf_size),
                vec![], // No merkle path needed
            ));
            if beta >= l[1].size() {
                beta -= l[1].size();
            }

            // For i=1..mu-1, generate full merkle proofs
            for i in 1..mu {
                mt_proofs[t].push(open_merkle_tree_at_conjugate_points(&mt[i - 1], &v[i], beta));
                if beta >= l[i + 1].size() {
                    beta -= l[i + 1].size();
                }
            }
        }
        end_timer!(timer);

        // Additional proofs for individual mt0s
        let timer = start_timer!(|| "DeepFold.Mt0Proofs");
        let mut mt_proofs_for_mt0 = Vec::new();
        let idx = (0..num_poly).filter(
            |&i| (0..i).all(|j| polynomials[i] != polynomials[j])
        ).collect::<Vec<_>>();
        for t in 0..s {
            mt_proofs_for_mt0.push(Vec::new());
            let x0 = mt_proofs[t][0].0;

            for (ki, &k) in idx.iter().enumerate() {
                let leaf_size = mt0_list[k].leaf_size();
                let step = l0.size() / leaf_size;
                let local_x0 = x0 % step;
                mt_proofs_for_mt0[t].push((
                    get_leaf_elements(&advices[k].v0, local_x0, step, leaf_size),
                    mt0_list[k].prove(local_x0),
                ));
            }
        }

        let evals = polynomials
            .iter()
            .zip(points.iter())
            .map(|(poly, point)| eval_mle_poly(&poly.evaluations, point))
            .collect::<Vec<_>>();
        end_timer!(timer);

        Ok(Self::BatchProof {
            deepfold_proof: DeepFoldProof {
                linear_polys,
                mt_roots,
                f_mu,
                mt_proofs,
            },
            sum_check_proof,
            mt_proofs_for_mt0,
            evals,
            sum_check_evals,
        })
    }

    fn d_batch_open(
        prover_param: impl Borrow<Self::ProverParam>,
        polynomials: Vec<Self::Polynomial>,
        advices: &[&Self::ProverCommitmentAdvice],
        points: &[Self::Point],
        _evals: &[Self::Evaluation],
        transcript: &mut IOPTranscript<F>,
    ) -> Result<Option<Self::BatchProof>, PCSError> {
        let &Self::ProverParam { max_mu, l0, s } = prover_param.borrow();
        let num_party = Net::n_parties();
        let num_party_vars = num_party.ilog2() as usize;
        let num_poly = polynomials.len();

        // Each party has local evaluations of size 2^local_mu
        let local_mu = polynomials[0].num_vars;
        let mu = local_mu + num_party_vars;
        assert!(mu <= max_mu);
        assert!(polynomials.iter().all(|poly| poly.num_vars == local_mu));
        assert!(points.iter().all(|point| point.len() == mu));
        assert!(points.len() == num_poly && advices.len() == num_poly);

        // Step 1: Gather all polynomial evaluations to master for SumCheck
        let timer = start_timer!(|| "DBatchOpen.GatherEvals");
        let all_poly_evals: Vec<Option<Vec<Vec<F>>>> = polynomials
            .iter()
            .map(|poly| Net::send_to_master(&poly.evaluations))
            .collect();
        end_timer!(timer);

        // Initialize structures for all parties
        let mt0_list = advices.iter().map(|advice| &advice.mt0).collect::<Vec<_>>();
        let mut l = vec![l0];
        l.append(
            &mut (1..mu + 1)
                .map(|i| GeneralEvaluationDomain::<F>::new(l0.size() >> i).unwrap())
                .collect::<Vec<_>>(),
        );

        // All parties need these for distributed Merkle tree construction
        let mut mt: Vec<MerkleTree> = Vec::new();
        let mut mt_roots: Vec<Byte32> = Vec::new();

        // Master-only data
        let mut a = vec![Vec::new()];
        let mut f_tilde: Vec<Vec<F>> = Vec::new();
        let mut f: Vec<Vec<F>> = Vec::new();
        let mut alpha = vec![F::ZERO];
        let mut linear_polys = Vec::new();
        let mut v: Vec<Vec<F>> = Vec::new();
        let mut f_mu = F::ZERO;
        let mut sum_check_proof: Option<IOPProof<F>> = None;
        let mut sum_check_evals: Vec<F> = Vec::new();
        let mut point: Vec<F> = Vec::new();
        let mut gamma: Vec<F> = Vec::new();
        let mut full_poly_evals: Vec<Vec<F>> = Vec::new(); // Save for evals computation

        let timer = start_timer!(|| "DBatchOpen.Compute");
        if Net::am_master() {
            // Reconstruct full polynomials on master
            full_poly_evals = all_poly_evals
                .into_iter()
                .map(|evals_opt| {
                    let all_evals: Vec<Vec<F>> = evals_opt.unwrap();
                    all_evals.into_iter().flatten().collect()
                })
                .collect();

            // SumCheck Phase
            let r = transcript.get_and_append_challenge(b"batched_sumcheck")?;
            let mut sum_check = VirtualPolynomial::new(mu);
            for i in 0..num_poly {
                sum_check
                    .add_mle_list(
                        [
                            evals_to_arcpoly(&full_poly_evals[i]),
                            evals_to_arcpoly(&get_tensor(&points[i])),
                        ],
                        r.pow([i as u64]),
                    )
                    .map_err(|e| PCSError::VirtualPolynomialError(format!("{:?}", e)))?;
            }
            let sc_proof = <PolyIOP<F> as SumCheck<F>>::prove(sum_check, transcript)
                .map_err(|e| PCSError::SumCheckError(format!("{:?}", e)))?;
            point = sc_proof.point.clone();
            sum_check_evals = full_poly_evals
                .iter()
                .map(|evals| eval_mle_poly(evals, &point))
                .collect();
            sum_check_proof = Some(sc_proof);

            // Batched Open Phase
            gamma = transcript.get_and_append_challenge_vectors(b"gamma", num_poly)?;
            let poly_evals: Vec<F> = (0..1 << mu)
                .map(|i| {
                    (0..num_poly)
                        .map(|j| gamma[j] * full_poly_evals[j][i])
                        .sum::<F>()
                })
                .collect();

            // Compute f0 and v0 for combined polynomial
            let f0 = evals_to_coeffs(mu, &poly_evals);
            let v0 = l0.fft(&f0);

            f_tilde.push(poly_evals);
            f.push(f0);
            v.push(v0);
            a[0].push(point.clone());
        }

        // Step 2: Main loop - all parties participate in dMerkle
        for i in 1..mu + 1 {
            // Get challenges (all parties must sync on transcript)
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);

            if Net::am_master() {
                a[i - 1].push(get_alpha_powers::<F>(alpha[i], mu - i + 1));
                let (f0_split, f1) = split_even_odd(&f_tilde[i - 1]);
                let (fe, fo) = split_even_odd(&f[i - 1]);

                // Compute linear_polys
                if i == mu {
                    linear_polys.push(vec![(f_tilde[i - 1][0], f_tilde[i - 1][1])]);
                } else {
                    linear_polys.push(
                        a[i - 1]
                            .iter()
                            .map(|w| {
                                assert!(!w.is_empty());
                                let w_tensor = get_tensor(&w[1..].to_vec());
                                (inner_product(&w_tensor, &f0_split), inner_product(&w_tensor, &f1))
                            })
                            .collect::<Vec<_>>(),
                    );
                    a.push(a[i - 1].iter().map(|w| w[1..].to_vec()).collect::<Vec<_>>());
                }

                // Get r challenge
                let ri = transcript.get_and_append_challenge(b"r")?;

                // Compute f[i] and f_tilde[i]
                f.push(vector_add(&fe, &scalar_vector_product(ri, &fo)));
                f_tilde.push(vector_add(
                    &scalar_vector_product(F::ONE - ri, &f0_split),
                    &scalar_vector_product(ri, &f1),
                ));

                // Compute v[i] = FFT(f[i])
                let vi = l[i].fft(&f[i]);
                v.push(vi.clone());

                if i == mu {
                    f_mu = v[i][0];
                } else {
                    // Build full Merkle tree on master for proof generation
                    let mti = build_merkle_tree(&vi);
                    mt_roots.push(mti.root().clone());
                    mt.push(mti);

                    // Workers still need to sync on leaf_size broadcast
                    let leaf_size = LEAF_SIZE.min(vi.len());
                    Net::recv_from_master_uniform(Some(leaf_size));
                }
            } else {
                // Workers: sync on transcript challenges
                let _ri = transcript.get_and_append_challenge(b"r")?;

                if i != mu {
                    // Receive leaf_size to stay in sync
                    let _leaf_size: usize = Net::recv_from_master_uniform(None);

                    mt_roots.push(Byte32::default());
                    mt.push(MerkleTree::default());
                }
            }
        }
        end_timer!(timer);

        // Only master generates proofs
        if !Net::am_master() {
            return Ok(None);
        }

        // Step 4: Generate merkle proofs
        let timer = start_timer!(|| "DBatchOpen.GenProofs");
        let mut mt_proofs = Vec::new();
        for t in 0..s {
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, l[0].size())?[0];
            mt_proofs.push(Vec::new());

            // For i=0, store values without merkle proof
            let leaf_size = LEAF_SIZE.min(v[0].len());
            let step = v[0].len() / leaf_size;
            let local_beta = beta % step;
            let beta_prime = if beta >= v[0].len() / 2 {
                beta - v[0].len() / 2
            } else {
                beta + v[0].len() / 2
            };
            mt_proofs[t].push((
                beta,
                (v[0][beta], v[0][beta_prime]),
                get_leaf_elements(&v[0], local_beta, step, leaf_size),
                vec![],
            ));
            if beta >= l[1].size() {
                beta -= l[1].size();
            }

            for i in 1..mu {
                mt_proofs[t].push(open_merkle_tree_at_conjugate_points(&mt[i - 1], &v[i], beta));
                if beta >= l[i + 1].size() {
                    beta -= l[i + 1].size();
                }
            }
        }
        end_timer!(timer);

        // Additional proofs for individual mt0s
        let timer = start_timer!(|| "DBatchOpen.Mt0Proofs");
        let mut mt_proofs_for_mt0 = Vec::new();
        // Deduplicate based on mt0 roots (same as verifier uses commitments)
        let idx = (0..num_poly)
            .filter(|&i| (0..i).all(|j| advices[i].mt0.root() != advices[j].mt0.root()))
            .collect::<Vec<_>>();
        for t in 0..s {
            mt_proofs_for_mt0.push(Vec::new());
            let x0 = mt_proofs[t][0].0;

            for (_, &k) in idx.iter().enumerate() {
                let leaf_size = mt0_list[k].leaf_size();
                let step = l0.size() / leaf_size;
                let local_x0 = x0 % step;
                mt_proofs_for_mt0[t].push((
                    get_leaf_elements(&advices[k].v0, local_x0, step, leaf_size),
                    mt0_list[k].prove(local_x0),
                ));
            }
        }

        // Compute evals: evaluation of each polynomial at its corresponding point
        let evals: Vec<F> = full_poly_evals
            .iter()
            .zip(points.iter())
            .map(|(poly_evals, pt)| eval_mle_poly(poly_evals, pt))
            .collect();

        end_timer!(timer);

        Ok(Some(Self::BatchProof {
            deepfold_proof: DeepFoldProof {
                linear_polys,
                mt_roots,
                f_mu,
                mt_proofs,
            },
            sum_check_proof: sum_check_proof.unwrap(),
            mt_proofs_for_mt0,
            evals,
            sum_check_evals,
        }))
    }

    fn verify(
        verifier_param: &Self::VerifierParam,
        com: &Self::Commitment,
        point: &Self::Point,
        value: &F,
        proof: &Self::Proof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        let Self::VerifierParam {
            max_mu,
            len_l0,
            g,
            s,
        } = verifier_param.clone();
        let Self::Commitment { mu, rt0 } = com.clone();
        // let mu = max_mu;
        // let point = resize_point(&point, mu);
        assert!(mu <= max_mu);
        let Self::Proof {
            linear_polys,
            mt_roots,
            f_mu,
            mt_proofs,
        } = proof.clone();

        if rt0 != mt_roots[0] {
            eprintln!("VERIFY FAIL: rt0 != mt_roots[0]");
            return Ok(false);
        }

        let mut alpha = vec![F::ZERO];
        let mut r = vec![F::ZERO];

        for _ in 1..mu + 1 {
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);
            r.push(transcript.get_and_append_challenge(b"r")?);
        }

        if eval_linear_poly(&linear_polys[0][0], &point[0]) != *value
            || eval_linear_poly(&linear_polys[mu - 1][0], &r[mu]) != f_mu
        {
            eprintln!("VERIFY FAIL: linear poly check");
            eprintln!("  eval_linear_poly(&linear_polys[0][0], &point[0])={:?}", eval_linear_poly(&linear_polys[0][0], &point[0]));
            eprintln!("  value={:?}", value);
            eprintln!("  eval_linear_poly(&linear_polys[mu - 1][0], &r[mu])={:?}", eval_linear_poly(&linear_polys[mu - 1][0], &r[mu]));
            eprintln!("  f_mu={:?}", f_mu);
            return Ok(false);
        }

        for i in 1..mu {
            for j in 0..linear_polys[i - 1].len() {
                let k = if i < mu - 1 { j } else { 0 };
                let w1 = if j == 0 {
                    point[i]
                } else {
                    alpha[j].pow([1 << (i + 1 - j) as u64])
                };
                if eval_linear_poly(&linear_polys[i - 1][j], &r[i])
                    != eval_linear_poly(&linear_polys[i][k], &w1)
                {
                    eprintln!("VERIFY FAIL: linear poly consistency check at i={}, j={}", i, j);
                    return Ok(false);
                }
            }
        }

        for t in 0..s {
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, len_l0)?[0];
            let mut beta_point = g.pow([beta as u64]);
            for i in 0..mu {
                let offset = len_l0 >> (i + 1);
                if !verify_merkle_tree_at_conjugate_points(
                    len_l0 >> i,
                    &mt_roots[i],
                    beta,
                    &mt_proofs[t][i].1,
                    &mt_proofs[t][i].2,
                    &mt_proofs[t][i].3,
                ) {
                    eprintln!("VERIFY FAIL: merkle proof at t={}, i={}, beta={}", t, i, beta);
                    return Ok(false);
                }

                let next_beta = if beta >= offset { beta - offset } else { beta };
                let val = if i < mu - 1 {
                    mt_proofs[t][i + 1].1.0
                } else {
                    f_mu
                };

                if !is_collinear(
                    (beta_point, mt_proofs[t][i].1 .0),
                    (-beta_point, mt_proofs[t][i].1 .1),
                    (r[i + 1], val),
                ) {
                    eprintln!("VERIFY FAIL: collinear check at t={}, i={}", t, i);
                    return Ok(false);
                }

                beta = next_beta;
                beta_point *= beta_point;
            }
        }

        Ok(true)
    }

    fn batch_verify(
        verifier_param: &Self::VerifierParam,
        commitments: &[Self::Commitment],
        points: &[Self::Point],
        batch_proof: &Self::BatchProof,
        transcript: &mut IOPTranscript<F>,
    ) -> Result<bool, PCSError> {
        let Self::VerifierParam {
            max_mu,
            len_l0,
            g,
            s,
        } = verifier_param.clone();
        let mu = max_mu;
        assert!(commitments.iter().all(|com| com.mu == mu));
        let num_poly = commitments.len();
        assert!(points.iter().all(|point| point.len() == mu));
        assert!(points.len() == num_poly);
        let Self::BatchProof {
            deepfold_proof,
            sum_check_proof,
            mt_proofs_for_mt0,
            evals,
            sum_check_evals,
        } = batch_proof.clone();

        let DeepFoldProof {
            ref linear_polys,
            ref mt_roots,
            f_mu,
            ref mt_proofs,
        } = deepfold_proof;

        // Sumcheck Phase
        let r_batch = transcript.get_and_append_challenge(b"batched_sumcheck")?;
        let sum_check_sum = <PolyIOP<F> as SumCheck<F>>::extract_sum(&sum_check_proof);
        if sum_check_sum
            != (0..num_poly)
                .map(|k| r_batch.pow([k as u64]) * evals[k])
                .sum::<F>()
        {
            return Ok(false);
        }
        let sum_check_claim = <PolyIOP<F> as SumCheck<F>>::verify(
            sum_check_sum,
            &sum_check_proof,
            &VPAuxInfo {
                max_degree: 2,
                num_variables: mu,
                phantom: PhantomData::<F>::default(),
            },
            transcript,
        )
        .map_err(|e| PCSError::SumCheckError(format!("{:?}", e)))?;
        let point = sum_check_proof.point.clone();
        if sum_check_claim.expected_evaluation
            != (0..num_poly)
                .map(|k| r_batch.pow([k as u64]) * eval_mle_eq(&point, &points[k]) * sum_check_evals[k])
                .sum::<F>()
        {
            return Ok(false);
        }

        // Batched Open Phase
        let gamma = transcript.get_and_append_challenge_vectors(b"gamma", num_poly)?;
        let value: F = (0..num_poly)
            .map(|k| gamma[k] * sum_check_evals[k])
            .sum();

        // Get challenges (same as verify())
        let mut alpha = vec![F::ZERO];
        let mut r = vec![F::ZERO];
        for _ in 1..mu + 1 {
            alpha.push(transcript.get_and_append_challenge(b"alpha")?);
            r.push(transcript.get_and_append_challenge(b"r")?);
        }

        // Verify linear polynomial relationships
        if eval_linear_poly(&linear_polys[0][0], &point[0]) != value
            || eval_linear_poly(&linear_polys[mu - 1][0], &r[mu]) != f_mu
        {
            return Ok(false);
        }

        for i in 1..mu {
            for j in 0..linear_polys[i - 1].len() {
                let k = if i < mu - 1 { j } else { 0 };
                let w1 = if j == 0 {
                    point[i]
                } else {
                    alpha[j].pow([1 << (i + 1 - j) as u64])
                };
                if eval_linear_poly(&linear_polys[i - 1][j], &r[i])
                    != eval_linear_poly(&linear_polys[i][k], &w1)
                {
                    return Ok(false);
                }
            }
        }

        // Verify merkle proofs and collinearity (skip mt0 merkle verification)
        for t in 0..s {
            let mut beta = transcript.get_and_append_challenge_indices(b"beta", 1, len_l0)?[0];
            let mut beta_point = g.pow([beta as u64]);

            for i in 0..mu {
                let offset = len_l0 >> (i + 1);

                // For i=0: skip merkle verification (will be checked via linear combination)
                // For i>=1: verify merkle proof using mt_roots[i-1]
                if i > 0 {
                    if !verify_merkle_tree_at_conjugate_points(
                        len_l0 >> i,
                        &mt_roots[i - 1],  // mt_roots shifted by 1
                        beta,
                        &mt_proofs[t][i].1,
                        &mt_proofs[t][i].2,
                        &mt_proofs[t][i].3,
                    ) {
                        return Ok(false);
                    }
                }

                let next_beta = if beta >= offset { beta - offset } else { beta };
                let val = if i < mu - 1 {
                    mt_proofs[t][i + 1].1.0
                } else {
                    f_mu
                };

                if !is_collinear(
                    (beta_point, mt_proofs[t][i].1.0),
                    (-beta_point, mt_proofs[t][i].1.1),
                    (r[i + 1], val),
                ) {
                    return Ok(false);
                }

                beta = next_beta;
                beta_point *= beta_point;
            }
        }

        // Additional checks for individual mt0s via linear combination
        let idx = (0..num_poly).filter(
            |&i| (0..i).all(|j| commitments[i] != commitments[j])
        ).collect::<Vec<_>>();
        let mut flag = vec![]; let mut cnt = 0;
        for i in 0..num_poly {
            for j in 0..=i {
                if i == j {
                    flag.push(cnt); cnt += 1;
                } else if commitments[i] == commitments[j] {
                    flag.push(flag[j]); break;
                }
            }
        }
        for t in 0..s {
            let mut sum = F::ZERO;
            let x = mt_proofs[t][0].0;
            for (ki, &k) in idx.iter().enumerate() {
                let leaf_size = mt_proofs_for_mt0[t][ki].0.len();
                let step = len_l0 / leaf_size;
                if !MerkleTree::verify(
                    &commitments[k].rt0,
                    x % step,
                    &compute_sha256_row(&mt_proofs_for_mt0[t][ki].0),
                    &mt_proofs_for_mt0[t][ki].1,
                ) {
                    return Ok(false);
                }
            }
            // Use the first proof's leaf_size for the sum computation
            let leaf_size = mt_proofs_for_mt0[t][0].0.len();
            let step = len_l0 / leaf_size;
            for (k, &ki) in flag.iter().enumerate() {
                sum += gamma[k] * mt_proofs_for_mt0[t][ki].0[x / step];
            }
            if sum != mt_proofs[t][0].1.0 {
                return Ok(false);
            }
        }

        Ok(true)
    }
}
