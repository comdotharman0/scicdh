use scicdh::probability::*;
#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    fn uniform_rv() -> RandomVariable {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let px = vec![0.2, 0.2, 0.2, 0.2, 0.2];
        RandomVariable::new(x, px)
    }

    // ===== Probability =====

    #[test]
    fn test_factorial_zero() {
        assert!(approx_eq(Probability::factorial(0).unwrap(), 1.0));
    }

    #[test]
    fn test_factorial_five() {
        // 1*2*3*4*5 = 120
        assert!(approx_eq(Probability::factorial(5).unwrap(), 120.0));
    }

    #[test]
    fn test_factorial_over_limit_errors() {
        assert!(Probability::factorial(102).is_err());
    }

    #[test]
    fn test_n_p_r() {
        // P(5,2) = 5*4 = 20
        assert!(approx_eq(Probability::n_p_r(5, 2).unwrap(), 20.0));
    }

    #[test]
    fn test_n_p_r_zero_n_errors() {
        assert!(Probability::n_p_r(0, 2).is_err());
    }

    #[test]
    fn test_n_c_r() {
        // C(5,2) = 5!/(2!*3!) = 120/12 = 10
        assert!(approx_eq(Probability::n_c_r(5, 2).unwrap(), 10.0));
    }

    #[test]
    fn test_get_probability() {
        // 3 favourable out of 6 = 0.5
        assert!(approx_eq(
            Probability::get_probability(6.0, 3.0).unwrap(),
            0.5
        ));
    }

    #[test]
    fn test_distribution_function() {
        // probs = [0.2, 0.3, 0.5]
        // CDF   = [0.2, 0.5, 1.0]
        let probs = vec![0.2, 0.3, 0.5];
        let cdf = Probability::distribution_function(&probs).unwrap();
        assert!(approx_eq(cdf[0], 0.2));
        assert!(approx_eq(cdf[1], 0.5));
        assert!(approx_eq(cdf[2], 1.0));
    }

    #[test]
    fn test_random_variable_func() {
        // f(x) = x^2 applied to [1,2,3] = [1,4,9]
        let set = vec![1.0, 2.0, 3.0];
        let result = Probability::random_variable(
            &set, |x| x.powf(2.0)
        ).unwrap();
        assert!(approx_eq(result[0], 1.0));
        assert!(approx_eq(result[1], 4.0));
        assert!(approx_eq(result[2], 9.0));
    }

    // ===== RandomVariable =====
    // Uniform on {0,1,2,3,4}, P(x)=0.2 each
    // E[X] = 0.2*(0+1+2+3+4) = 2.0
    // E[X^2] = 0.2*(0+1+4+9+16) = 6.0
    // Var = 6.0 - 4.0 = 2.0

    #[test]
    fn test_rv_mean() {
        assert!(approx_eq(uniform_rv().mean().unwrap(), 2.0));
    }

    #[test]
    fn test_rv_variance() {
        assert!(approx_eq(uniform_rv().variance().unwrap(), 2.0));
    }

    #[test]
    fn test_rv_std_deviation() {
        // sqrt(2.0)
        assert!(approx_eq(
            uniform_rv().std_deviation().unwrap(),
            2.0_f64.sqrt()
        ));
    }

    #[test]
    fn test_rv_mgf_at_zero() {
        // MGF(0) = E[e^0] = E[1] = 1.0 always
        assert!(approx_eq(
            uniform_rv().moment_generating_func(0.0).unwrap(),
            1.0
        ));
    }

    #[test]
    fn test_rv_expectation_from_func() {
        // E[X^2] = 0.2*(0+1+4+9+16) = 6.0
        let result = uniform_rv()
            .expectation_from_func(&|x:f64| x.powf(2.0))
            .unwrap();
        assert!(approx_eq(result, 6.0));
    }

    #[test]
    fn test_rv_expectation_from_set() {
        // h(X) = [0,1,4,9,16], P = [0.2,0.2,0.2,0.2,0.2]
        // E = 0.2*(0+1+4+9+16) = 6.0
        let hx = vec![0.0, 1.0, 4.0, 9.0, 16.0];
        let result = uniform_rv()
            .expectation_from_set(&hx)
            .unwrap();
        assert!(approx_eq(result, 6.0));
    }

    #[test]
    fn test_characteristic_func_errors() {
        // Must return Err — complex numbers not supported
        #[allow(deprecated)]
        let result = uniform_rv().characteristic_func(1.0);
        assert!(result.is_err());
    }

    // ===== Binomial =====
    // B(n=10, p=0.5)
    // Mean = 10*0.5 = 5.0
    // Variance = 10*0.5*0.5 = 2.5
    // P(X=0) = C(10,0)*0.5^0*0.5^10 = 1/1024

    #[test]
    fn test_binomial_mean() {
        let b = Binomial::new(10, vec![0], 0.5);
        assert!(approx_eq(b.mean().unwrap(), 5.0));
    }

    #[test]
    fn test_binomial_variance() {
        let b = Binomial::new(10, vec![0], 0.5);
        assert!(approx_eq(b.variance().unwrap(), 2.5));
    }

    #[test]
    fn test_binomial_pmf_at_zero() {
        // P(X=0) = (0.5)^10 = 1/1024
        let b = Binomial::new(10, vec![0], 0.5);
        let probs = b.get_probability_set().unwrap();
        assert!(approx_eq(probs[0], 1.0 / 1024.0));
    }

    #[test]
    fn test_binomial_pmf_sums_to_one() {
        // Full distribution B(4, 0.5): x=0,1,2,3,4
        // Must sum to 1.0 — law of total probability
        let b = Binomial::new(4, vec![0, 1, 2, 3, 4], 0.5);
        let sum: f64 = b.get_probability_set()
            .unwrap()
            .iter()
            .sum();
        assert!(approx_eq(sum, 1.0));
    }

    // ===== Geometric =====
    // Geometric(p=0.5)
    // Mean = 1/0.5 = 2.0
    // Variance = 0.5/0.25 = 2.0
    // P(X=1) = 0.5^1 * 0.5^0 = 0.5

    #[test]
    fn test_geometric_mean() {
        let g = Geometric::new(vec![1], 0.5);
        assert!(approx_eq(g.mean().unwrap(), 2.0));
    }

    #[test]
    fn test_geometric_variance() {
        let g = Geometric::new(vec![1], 0.5);
        assert!(approx_eq(g.variance().unwrap(), 2.0));
    }

    #[test]
    fn test_geometric_pmf_at_one() {
        // P(X=1) = p*(1-p)^0 = 0.5
        let g = Geometric::new(vec![1], 0.5);
        let probs = g.get_probability_set().unwrap();
        assert!(approx_eq(probs[0], 0.5));
    }

    // ===== Pascal =====
    // Pascal(r=2, p=0.5)
    // Mean = r/p = 2/0.5 = 4.0  (YOU calculated this)
    // Variance = r/p^2 = 2/0.25 = 8.0  (YOU calculated this)
    // P(X=2) = C(1,1)*0.5^2*0.5^0 = 0.25

    #[test]
    fn test_pascal_mean() {
        let p = Pascal::new(2, vec![2, 3], 0.5);
        assert!(approx_eq(p.mean().unwrap(), 4.0));
    }

    #[test]
    fn test_pascal_variance() {
        let p = Pascal::new(2, vec![2, 3], 0.5);
        assert!(approx_eq(p.variance().unwrap(), 8.0));
    }

    #[test]
    fn test_pascal_pmf_at_r() {
        // P(X=r) = C(r-1,r-1)*p^r*(1-p)^0 = p^r
        // For r=2, p=0.5: P(X=2) = 0.5^2 = 0.25
        let p = Pascal::new(2, vec![2], 0.5);
        let probs = p.get_probability_set().unwrap();
        assert!(approx_eq(probs[0], 0.25));
    }

    // ===== HyperGeometric =====
    // N=10, T=4, n=3, x=1
    // Mean = n*T/N = 3*4/10 = 1.2  (YOU calculated this)
    // Variance = mean*(1-T/N)*((N-n)/(N-1))
    //          = 1.2*(1-0.4)*(7/9)
    //          = 1.2*0.6*0.7778
    //          = 0.56

    #[test]
    fn test_hypergeometric_mean() {
        let h = HyperGeometric::new(10, 4, 3, 1);
        assert!(approx_eq(h.mean().unwrap(), 1.2));
    }

    #[test]
    fn test_hypergeometric_variance() {
        // 1.2 * 0.6 * (7.0/9.0) = 0.5600000...
        let h = HyperGeometric::new(10, 4, 3, 1);
        let expected = 1.2 * (1.0 - 4.0/10.0) * (7.0/9.0);
        assert!(approx_eq(h.variance().unwrap(), expected));
    }

    #[test]
    fn test_hypergeometric_pmf() {
        // P(X=1) = C(4,1)*C(6,2)/C(10,3)
        // = 4*15/120 = 60/120 = 0.5
        let h = HyperGeometric::new(10, 4, 3, 1);
        assert!(approx_eq(h.get_probability().unwrap(), 0.5));
    }

    // ===== RandomVector =====

    #[test]
    fn test_random_vector_rejects_single_rv() {
        let rv = uniform_rv();
        let jp = JointProbability::new(vec![vec![1.0]]);
        assert!(RandomVector::new(vec![rv], jp).is_err());
    }

    #[test]
    fn test_random_vector_length_mismatch_rejected() {
        // Tests the bug fix: !i.variable.len()==len was always false
        // After fixing to i.variable.len() != len, this must return Err
        let rv1 = RandomVariable::new(
            vec![0.0, 1.0],
            vec![0.5, 0.5]
        );
        let rv2 = RandomVariable::new(
            vec![0.0, 1.0, 2.0],
            vec![0.33, 0.33, 0.34]
        );
        let jp = JointProbability::new(vec![vec![1.0]]);
        assert!(RandomVector::new(vec![rv1, rv2], jp).is_err());
    }

    #[test]
    fn test_joint_probability_validates_sum() {
        // 4 elements each 0.25 = sum 1.0
        let jp = JointProbability::new(vec![
            vec![0.25, 0.25],
            vec![0.25, 0.25],
        ]);
        let sum = jp.matrix.iter()
            .fold(0.0_f64, |acc, row| {
                acc + row.iter().sum::<f64>()
            });
        assert!(approx_eq(sum, 1.0));
    }
}
