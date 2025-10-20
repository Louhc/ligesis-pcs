use ark_ff::{Field, PrimeField};
use ark_std::{
    marker::PhantomData, cfg_into_iter
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use sha2::{digest, Digest, Sha256};

pub type Byte32 = [u8; 32];

pub mod sis;
pub use sis::*;

pub fn decompose<F: Field>(data: &[F]) -> Vec<u8> {
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
    let result = hasher.finalize();
    result.into()
}

pub fn compute_sha256_row<F: Field>(data: &[F]) -> Byte32 {
    compute_sha256(&decompose(data))
}

pub trait Hasher<I, O, C> {
    fn new( config: &C ) -> Self;
    fn compute( &self, data : &I ) -> O;
}

pub struct MTSHA256Hasher;
impl Hasher<(Byte32, Byte32), Byte32, ()> for MTSHA256Hasher{
    fn new( _: &() ) -> Self {
        MTSHA256Hasher{}
    }
    fn compute( &self, data: &(Byte32, Byte32) ) -> Byte32 {
        let (data0, data1) = data;
        let mut d = data0.to_vec();
        d.extend_from_slice(data1);
        compute_sha256(&d)
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone, Debug, PartialEq, Eq)]
pub struct MerkleTree {
    n: usize,
    digest_layers: Vec<Vec<Byte32>>,
}

impl Default for MerkleTree {
    fn default() -> Self {
        MerkleTree { n: 0, digest_layers: Vec::<Vec<Byte32>>::new() }
    }
}

pub type MerkleTreeProof = Vec<Byte32>;

impl MerkleTree {
    pub fn new( leaves: &Vec<Byte32> ) -> Self {
        let mut digest_layers: Vec<Vec<Byte32>> = Vec::new();
        let n = if leaves.len() == 1 { 1usize } else { ((leaves.len() - 1).ilog2() + 2) as usize };
        let hasher = MTSHA256Hasher::new(&());
        
        digest_layers.push({
            let mut digest_layer= leaves.clone();
            digest_layer.resize(1 << (n - 1), [0u8; 32]);
            digest_layer
        });

        for i in 1..n {
            digest_layers.push({
                let mut digest_layer: Vec<Byte32> = Vec::new();
                for j in 0..(1usize << (n - i - 1)) {
                    digest_layer.push(hasher.compute(&(
                        digest_layers[i - 1][j << 1], 
                        digest_layers[i - 1][j << 1 | 1]
                    )));
                }
                digest_layer
            });
        }
        
        Self { n, digest_layers }
    }

    pub fn root( &self ) -> Byte32 {
        self.digest_layers.last().unwrap().first().unwrap().clone()
    }

    pub fn prove( &self, pos: usize ) -> Vec<Byte32> {
        let mut proof: Vec<Byte32> = Vec::new();
        let mut j = pos;
        for i in 0..(self.n - 1) {
            proof.push(self.digest_layers[i][j ^ 1]);
            j >>= 1;
        }
        proof
    }

    pub fn verify( root: &Byte32, pos: usize, val: &Byte32, proof: &MerkleTreeProof ) -> bool {
        let hasher = MTSHA256Hasher::new(&());
        let mut now = val.clone();
        let mut j = pos;
        for i in 0..proof.len() {
            now = if j & 1 == 0 { hasher.compute(&(now, proof[i])) } 
                else { hasher.compute(&(proof[i], now)) };
            j >>= 1;
        }
        root == &now
    }
}

#[cfg(test)]
mod tests;