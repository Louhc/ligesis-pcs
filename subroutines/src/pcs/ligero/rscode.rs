use ark_ff::FftField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

pub struct ReedSolomon<F: FftField> {
    n: usize,
    domain: GeneralEvaluationDomain<F>,
}

impl<F: FftField> ReedSolomon<F> {
    pub fn new(n: usize, ecc: usize) -> Self {
        let domain = GeneralEvaluationDomain::<F>::new(n + ecc).unwrap();
        Self { n, domain }
    }

    pub fn encode(&self, v: &Vec<F>) -> Vec<F> {
        self.domain.fft(v)
    }

    pub fn get_n(&self) -> usize {
        self.n
    }

    pub fn get_k(&self) -> usize {
        self.domain.size()
    }
}