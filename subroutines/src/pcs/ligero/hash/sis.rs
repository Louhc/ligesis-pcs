use super::*;
use ark_ff::PrimeField;

pub struct SISHasher<F: PrimeField>{
    n: usize,
    m: usize,
    b: F,
    mat_a: Vec<Vec<F>>,
}

impl<F: PrimeField> Hasher<Vec<F>, Vec<F>, (Vec<Vec<F>>, F)> for SISHasher<F> {
    fn new( config: &(Vec<Vec<F>>, F) ) -> Self {
        let (mat_a, b) = (config.0.clone(), config.1.clone());
        let (n, m) = (mat_a.len(), mat_a[0].len());
        assert!(mat_a.len() == n 
            && mat_a.iter().all(|row| row.len() == m));
        SISHasher{n, m, b, mat_a}
    }

    fn compute( &self, v: &Vec<F> ) -> Vec<F> {
        assert_eq!(v.len(), self.m);
        assert!(v.iter().all(|&x| x < self.b));
        (0..self.n).map(|i| (0..self.m).map(
            |j| self.mat_a[i][j] * v[j]
        ).sum()).collect()
    }
}

pub struct SISBitHasher<F: PrimeField>{
    n: usize,
    m: usize,
    mat_a: Vec<Vec<F>>,
}
impl<F: PrimeField> Hasher<Vec<u8>, Vec<F>, Vec<Vec<F>>> for SISBitHasher<F> {
    fn new( config: &Vec<Vec<F>> ) -> Self {
        let mat_a = config.clone();
        let (n, m) = (mat_a.len(), mat_a[0].len());
        assert!(mat_a.len() == n 
            && mat_a.iter().all(|row| row.len() == m));
        SISBitHasher{n, m, mat_a}
    }

    fn compute( &self, data: &Vec<u8> ) -> Vec<F> {
        let mut res = vec![F::ZERO; self.n];
        for i in 0..self.m {
            if data[i] == 1 {
                for j in 0..self.n {
                    res[j] += self.mat_a[j][i];
                }
            }
        }
        res
    }
}