use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use sha2::{Digest, Sha256};

pub type Byte32 = [u8; 32];

pub fn serialize<F: Field>(data: &[F]) -> Vec<u8> {
    let mut serialized = Vec::new();
    for element in data {
        element
            .serialize_with_mode(&mut serialized, ark_serialize::Compress::Yes)
            .expect("Serialization fails");
    }
    serialized
}

pub fn compute_sha256(data: &[u8]) -> Byte32 {
    let mut hasher = Sha256::new();
    hasher.update(data);
    hasher.finalize().into()
}

pub fn compute_sha256_row<F: Field>(data: &[F]) -> Byte32 {
    compute_sha256(&serialize(data))
}

fn hash_pair(left: &Byte32, right: &Byte32) -> Byte32 {
    let mut buf = [0u8; 64];
    buf[..32].copy_from_slice(left);
    buf[32..].copy_from_slice(right);
    compute_sha256(&buf)
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
pub struct MerkleTree {
    n: usize,
    digest_layers: Vec<Vec<Byte32>>,
}

impl Default for MerkleTree {
    fn default() -> Self {
        MerkleTree {
            n: 0,
            digest_layers: Vec::new(),
        }
    }
}

pub type MerkleTreeProof = Vec<Byte32>;

impl MerkleTree {
    pub fn new(leaves: &Vec<Byte32>) -> Self {
        let mut digest_layers: Vec<Vec<Byte32>> = Vec::new();
        let n = leaves.len().next_power_of_two().trailing_zeros() as usize + 1;

        digest_layers.push({
            let mut digest_layer = leaves.clone();
            digest_layer.resize(1 << (n - 1), [0u8; 32]);
            digest_layer
        });

        for i in 1..n {
            let prev = &digest_layers[i - 1];
            let digest_layer: Vec<Byte32> = (0..(1usize << (n - i - 1)))
                .map(|j| hash_pair(&prev[j << 1], &prev[j << 1 | 1]))
                .collect();
            digest_layers.push(digest_layer);
        }

        Self { n, digest_layers }
    }

    pub fn root(&self) -> Byte32 {
        *self.digest_layers.last().unwrap().first().unwrap()
    }

    pub fn prove(&self, pos: usize) -> Vec<Byte32> {
        let mut proof = Vec::with_capacity(self.n - 1);
        let mut j = pos;
        for i in 0..(self.n - 1) {
            proof.push(self.digest_layers[i][j ^ 1]);
            j >>= 1;
        }
        proof
    }

    pub fn verify(root: &Byte32, pos: usize, val: &Byte32, proof: &MerkleTreeProof) -> bool {
        let mut now = *val;
        let mut j = pos;
        for sibling in proof {
            now = if j & 1 == 0 {
                hash_pair(&now, sibling)
            } else {
                hash_pair(sibling, &now)
            };
            j >>= 1;
        }
        root == &now
    }
}

#[cfg(test)]
mod tests;
