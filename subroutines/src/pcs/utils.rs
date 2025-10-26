use ark_ff::{PrimeField, BigInteger};
use ark_poly::DenseMultilinearExtension;
use std::sync::Arc;

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
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| if b[k][j] { a[i][k] } else { F::ZERO }).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
}

pub fn bool_mat_mul_field_mat<F: PrimeField>( a: &Vec<Vec<bool>>, b: &Vec<Vec<F>> ) -> Vec<Vec<F>> {
    let n = a.len();
    let m = a[0].len();
    let p = b[0].len();
    assert_eq!(m, b.len());
    (0..n).map(
        |i| (0..p).map(
            |j| (0..m).map(|k| if a[i][k] { b[k][j] } else { F::ZERO }).sum::<F>()
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>()
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