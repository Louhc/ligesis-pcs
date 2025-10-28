use arithmetic::VirtualPolynomial;
use ark_ff::{PrimeField, BigInteger, BigInt};
use ark_poly::DenseMultilinearExtension;
use std::sync::Arc;

use crate::IOPProof;

pub fn get_alpha_powers<F: PrimeField>( alpha: F, mu: usize ) -> Vec<F> {
    let mut res = Vec::new();
    let mut t = alpha;
    for _ in 0..mu {
        res.push(t);
        t = t * t;
    }
    res
}

pub fn get_tensor<F: PrimeField>( r: &Vec<F> ) -> Vec<F> {
    let mut res = vec![F::ONE];
    for i in 0..r.len() {
        let mut new_res = Vec::new();
        for &x in res.iter() {
            new_res.push(x * (F::ONE - r[i]));
        }
        for &x in res.iter() {
            new_res.push(x * r[i]);
        }
        res = new_res;
    }
    res
}

pub fn split_even_odd<F: PrimeField>( v: &Vec<F> ) -> (Vec<F>, Vec<F>) {
    let mut even = Vec::new();
    let mut odd = Vec::new();
    for i in 0..v.len() {
        if i % 2 == 0 {
            even.push(v[i]);
        } else {
            odd.push(v[i]);
        }
    }
    (even, odd)
}

pub fn hadamard_product<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> Vec<F> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect()
}

pub fn inner_product<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> F {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

pub fn scalar_vector_product<F: PrimeField>( scalar: F, v: &Vec<F> ) -> Vec<F> {
    v.iter().map(|&x| scalar * x).collect()
}

pub fn vector_add<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> Vec<F> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}

pub fn evals_to_coeffs<F: PrimeField>( mu: usize, v: &Vec<F> ) -> Vec<F> {
    let mut u = v.clone();
    for j in 0..mu {
        for i in 0..(1 << mu) {
            if i & (1 << j) != 0 {
                u[i] = u[i] - u[i ^ (1 << j)];
            }
        }
    }
    u
}

pub fn eval_linear_poly<F: PrimeField>( f: &(F, F), point: &F ) -> F {
    f.0 * (F::ONE - *point) + f.1 * *point
}

pub fn is_collinear<F: PrimeField>( p0: (F, F), p1: (F, F), p2: (F, F) ) -> bool {
    let (x0, y0) = p0;
    let (x1, y1) = p1;
    let (x2, y2) = p2;
    return (y1 - y0) * (x2 - x1) == (y2 - y1) * (x1 - x0);
}

pub fn eval_univar_poly<F: PrimeField>( f: &Vec<F>, alpha: &F ) -> F {
    (0..f.len()).map(
        |i| f[i] * alpha.pow([i as u64])
    ).sum()
}

pub fn eval_mle_poly<F: PrimeField>( f: &Vec<F>, point: &Vec<F> ) -> F {
    inner_product(&f, &get_tensor(&point))
}

pub fn reshape<F: PrimeField>( a: &Vec<F>, n: usize, m: usize ) -> Vec<Vec<F>> {
    (0..n).map(
        |i| (0..m).map(
            |j| if i * m + j < a.len() { a[i * m + j] } else { F::ZERO }
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn transposition<F: Copy>( mat: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    (0..mat[0].len()).map(
        |i| (0..mat.len()).map(
            |j| mat[j][i]
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn decompose<F: PrimeField>( x: &F ) -> Vec<bool> {
    x.into_bigint().to_bits_be()
}

pub fn decompose_vector<F: PrimeField>( v: &Vec<F> ) -> Vec<bool> {
    v.iter().map(|x| decompose(x)).collect::<Vec<_>>().concat()
}

pub fn decompose_mat_by_col<F: PrimeField>( mat: &Vec<Vec<F>> ) -> Vec<Vec<bool>> {
    transposition(&transposition(&mat)
        .iter()
        .map(|col| decompose_vector(col))
        .collect::<Vec<_>>())
}

pub fn mat_mul<F: PrimeField>( a: &Vec<Vec<F>>, b: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert!(m == b.len());
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| a[i][k] * b[k][j]).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn field_mat_mul_bool_mat<F: PrimeField>( a: &Vec<Vec<F>>, b: &Vec<Vec<bool>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert_eq!(m, b.len());
    
    let c: Vec<Vec<u128>> = (0..n).map(
        |i| (0..m).map(
            |j| u64::from_le_bytes(a[i][j].into_bigint().to_bytes_le().try_into().unwrap()) as u128
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();

    let mut res = (0..n).map(|_| vec![0u128; p]).collect::<Vec<_>>();

    for j in 0..m {
        for k in 0..p {
            if b[j][k] {
                for i in 0..n {
                    res[i][k] = res[i][k] + c[i][j];
                }
            }
        }
    }

    (0..n).map(
        |i| (0..p).map(
            |k| F::from(res[i][k])
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn field_mat_mul_trans_bool_mat<F: PrimeField>( a: &Vec<Vec<F>>, b: &Vec<Vec<bool>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b.len();
    assert_eq!(m, b[0].len());
    
    let c: Vec<Vec<u128>> = (0..n).map(
        |i| (0..m).map(
            |j| u64::from_le_bytes(a[i][j].into_bigint().to_bytes_le().try_into().unwrap()) as u128
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();

    let mut res = (0..n).map(|_| vec![0u128; p]).collect::<Vec<_>>();

    for k in 0..p {
        for j in 0..m {
            if b[k][j] {
                for i in 0..n {
                    res[i][k] = res[i][k] + c[i][j];
                }
            }
        }
    }

    (0..n).map(
        |i| (0..p).map(
            |k| F::from(res[i][k])
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn bool_mat_mul_field_mat<F: PrimeField>( a: &Vec<Vec<bool>>, b: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert_eq!(m, b.len());

    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert_eq!(m, b.len());
    
    let c: Vec<Vec<u128>> = (0..n).map(
        |i| (0..m).map(
            |j| u64::from_le_bytes(b[i][j].into_bigint().to_bytes_le().try_into().unwrap()) as u128
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();


    let mut res = (0..n).map(|_| vec![0u128; p]).collect::<Vec<_>>();

    for i in 0..n {
        for j in 0..m {
            if a[i][j] {
                for k in 0..p {
                    res[i][k] = res[i][k] + c[j][k];
                }
            }
        }
    }

    (0..n).map(
        |i| (0..p).map(
            |k| F::from(res[i][k])
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()

    // (0..n).map(
    //     |i| (0..p).map(
    //         |j| (0..m).map(|k| if a[i][k] { b[k][j] } else { F::ZERO }).sum::<F>()
    //     ).collect::<Vec<_>>()
    // ).collect::<Vec<_>>()
}

pub fn evals_to_arcpoly<F: PrimeField>( a: &Vec<F> ) -> Arc<DenseMultilinearExtension<F>> {
    Arc::new(DenseMultilinearExtension::<F>::from_evaluations_vec(a.len().ilog2() as usize, a.clone()))
}

pub fn otimes<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> Vec<F> {
    a.iter().map(
        |x| b.iter().map(
            |y| (*x) * (*y)
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>().concat()
}

pub fn bool_vec_to_field_vec<F: PrimeField>( a: &Vec<bool> ) -> Vec<F> {
    (0..a.len()).map(
        |i| if a[i] {F::ONE} else {F::ZERO}
    ).collect::<Vec<_>>()
}

pub fn eval_mle_eq<F: PrimeField>( a: &Vec<F>, b: &Vec<F> ) -> F {
    a.iter().zip(b.iter()).map(
        |(&x, &y)| x * y + (F::ONE - x) * (F::ONE - y)
    ).product()
}

pub fn eval_mat_g_mle<F: PrimeField>( log_m: usize, log_n: usize, g: F, a: &Vec<F>, b: &Vec<F> ) -> F {
    let (m, n) = (1 << log_m, 1 << log_n);
    let ta = get_tensor(&a);
    let tb = get_tensor(&b);
    let mut res = F::ZERO;
    for i in 0..m {
        for j in 0..n {
            res += ta[i] * g.pow([(i * j) as u64]) * tb[j];
        }
    }
    res
}

pub fn compute_alpha_mat_g<F: PrimeField>( log_m: usize, log_n: usize, g: &F, alpha: &Vec<F> ) -> Vec<Vec<F>> {
    let mut alpha_mat_g = vec![vec![F::ONE]];
    for i in 0..log_m {
        let gi = g.pow([1u64 << (log_m - i - 1)]);
        let mut x = F::ONE;
        alpha_mat_g.push(Vec::new());
        for j in 0..(1 << (i + 1)) {
            let v = alpha_mat_g[i][j % (1 << i)] * (F::ONE - alpha[log_m - i - 1] + alpha[log_m - i - 1] * x);
            alpha_mat_g[i + 1].push(v);
            x *= gi;
        }
    }
    alpha_mat_g
}

#[cfg(test)]
mod tests {
    use crate::{random_field_vector_from_rng, ReedSolomon};

    use super::*;
    use ark_bls12_381::Fr as F;
    use ark_ff::Field;
    use ark_std::test_rng;

    #[test]
    fn test_compute_alpha_mat_g(){
        let mut rng = test_rng();
        let (log_m, log_n) = (6, 5);
        let (m, n) = (1 << log_m, 1 << log_n);
        let rs = ReedSolomon::<F>::new(n, m);
        let g = rs.get_generator();
        assert_eq!(g.pow([1u64 << log_m]), F::ONE);
        let mut alpha = random_field_vector_from_rng(log_m, &mut rng);
        let alpha_mat_g = compute_alpha_mat_g(log_m, log_n, &g, &alpha);

        let mat_g = (0..m).map(
            |i| (0..n).map(
                |j| g.pow([(i*j) as u64])
            ).collect::<Vec<_>>()
        ).collect::<Vec<_>>();

        let a = mat_mul(
            &vec![get_tensor(&alpha)], 
            &mat_g)[0].clone();

        assert_eq!(alpha_mat_g[log_m][..n].to_vec(), a);
    }
}