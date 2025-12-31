use crate::hash::*;
use ark_ff::PrimeField;

pub fn build_merkle_tree<F: PrimeField>(v: &Vec<F>) -> MerkleTree {
    assert!(v.len() >= 8);
    // let mut v = v.clone();
    // if v.len() < 8 { v.resize(8, F::ZERO); }
    let step = v.len() / 8;
    MerkleTree::new(
        &(0..v.len() / 8)
            .map(|i| compute_sha256_row(&[v[i], v[i + step], v[i + step * 2], v[i + step * 3], v[i + step * 4], v[i + step * 5], v[i + step * 6], v[i + step * 7]]))
            .collect(),
    )
}

pub fn open_merkle_tree_at_conjugate_points<F: PrimeField>(
    mt: &MerkleTree,
    v: &Vec<F>,
    x: usize,
) -> (usize, (F, F), [F; 8], Vec<Byte32>) {
    assert!(v.len() >= 8);
    // let mut v = v.clone();
    // if v.len() < 8 { v.resize(8, F::ZERO); }

    let step = v.len() / 8;
    let x0 = x % step;
    let x_prime = if x >= v.len() / 2 {
        x - v.len() / 2
    } else {
        x + v.len() / 2
    };
    (
        x.clone(),
        (v[x], v[x_prime]),
        [v[x0], v[x0 + step], v[x0 + step * 2], v[x0 + step * 3], v[x0 + step * 4], v[x0 + step * 5], v[x0 + step * 6], v[x0 + step * 7]],
        mt.prove(x0),
    )
}

pub fn verify_merkle_tree_at_conjugate_points<F: PrimeField>(
    n: usize,
    root: &Byte32,
    x: usize,
    v: &(F, F),
    w: &[F; 8],
    proof: &Vec<Byte32>,
) -> bool {
    assert!(n >= 8);
    // let n = if n < 8 { 8 } else { n };
    let step = n / 8;
    let x0 = x % step;
    let x_prime = if x >= n / 2 { x - n / 2 } else { x + n / 2 };
    if v.0 != w[x / step] || v.1 != w[x_prime / step] {
        assert!(false);
        return false;
    }
    MerkleTree::verify(root, x0, &compute_sha256_row(w), proof)
}
