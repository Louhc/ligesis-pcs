use crate::pcs::prelude::*;
use ark_ff::PrimeField;

pub fn build_merkle_tree<F: PrimeField>(v: &Vec<F>) -> MerkleTree {
    assert!(v.len() >= 4);
    let step = v.len() / 4;
    MerkleTree::new(
        &(0..v.len() / 4)
            .map(|i| compute_sha256_row(&[v[i], v[i + step], v[i + step * 2], v[i + step * 3]]))
            .collect(),
    )
}

pub fn open_merkle_tree_at_conjugate_points<F: PrimeField>(
    mt: &MerkleTree,
    v: &Vec<F>,
    x: usize,
) -> (usize, (F, F), [F; 4], Vec<Byte32>) {
    let step = v.len() / 4;
    let x0 = x % step;
    let x_prime = if x >= v.len() / 2 {
        x - v.len() / 2
    } else {
        x + v.len() / 2
    };
    (
        x.clone(),
        (v[x], v[x_prime]),
        [v[x0], v[x0 + step], v[x0 + step * 2], v[x0 + step * 3]],
        mt.prove(x0),
    )
}

pub fn verify_merkle_tree_at_conjugate_points<F: PrimeField>(
    n: usize,
    root: &Byte32,
    x: usize,
    v: &(F, F),
    w: &[F; 4],
    proof: &Vec<Byte32>,
) -> bool {
    let step = n / 4;
    let x0 = x % step;
    let x_prime = if x >= n / 2 { x - n / 2 } else { x + n / 2 };
    if v.0 != w[x / step] || v.1 != w[x_prime / step] {
        assert!(false);
        return false;
    }
    MerkleTree::verify(root, x0, &compute_sha256_row(w), proof)
}
