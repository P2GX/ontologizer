use std::cmp::{max, min};

pub struct Hypergeometric {
    log_factorials: Vec<f64>,
}

impl Hypergeometric {
    pub fn new() -> Self {
        let mut log_facts = vec![];
        // NOTE: log(0!) = 0 and log(1!) = 0
        log_facts.push(0.0); // log(0!)
        log_facts.push(0.0); // log(1!)
        Hypergeometric {
            log_factorials: log_facts,
        }
    }

    pub fn log_factorial(&mut self, i: usize) -> Result<f64, String> {
        // Ensure log_factorial[i] is computed; if not, extend the table.
        if i >= self.log_factorials.len() {
            for j in self.log_factorials.len()..=i {
                // current last entry equals ln((j-1)!)
                let last = *self.log_factorials.last().unwrap();
                // push ln(j!) = ln((j-1)!) + ln(j)
                self.log_factorials.push(last + (j as f64).ln());
            }
        }
        self.log_factorials
            .get(i)
            .copied()
            .ok_or_else(|| format!("Could not calculate log factorial for i={}", i))
    }

    /// Compute ln(n choose k) with simple validation.
    pub fn log_n_choose_k(&mut self, n: usize, k: usize) -> Result<f64, String> {
        let result = self.log_factorial(n)? - self.log_factorial(k)? - self.log_factorial(n - k)?;
        Ok(result)
    }

    #[allow(non_snake_case)]
    fn support(&mut self, n: usize, K: usize, N: usize) -> (usize, usize) {
        // Compute the support of the distributio, i.e. the range of values for which P(k) != 0.
        let k_min = max(0isize, (n as isize) + (K as isize) - (N as isize)) as usize;
        // Example: n=100, N=200, K=120 -> cannot draw less than 20 successes
        let k_max = min(n, K);
        // Example: n=100, K=80, -> cannot draw more than 80 successes
        (k_min, k_max)
    }

    /// dhyper: Probability density function of the hypergeometric distribution.
    ///
    /// Parameters
    /// - `k`: number of successes draws
    /// - `n`: number of draws
    /// - `K`: total number of success items
    /// - `N`: total number of items
    ///
    /// Returns `P(X = k)`, the probability of observing exactly `k` successes in `n` draws.
    #[allow(non_snake_case)]
    pub fn dhyper(&mut self, k: usize, n: usize, K: usize, N: usize) -> Result<f64, String> {
        // Domain checks
        if n > N || K > N {
            return Ok(0.0);
        }
        let (k_min, k_max) = self.support(n, K, N);
        if k < k_min || k > k_max {
            return Ok(0.0);
        }
        let log_prob = self.log_n_choose_k(K, k)? + self.log_n_choose_k(N - K, n - k)?
            - self.log_n_choose_k(N, n)?;
        Ok(log_prob.exp())
    }

    /// Forward term ratio r_fwd(k) = t(k+1)/t(k)
    /// where t(k) = P(X=k) for Hypergeometric(n, K, N)
    #[allow(non_snake_case)]
    fn ratio_forward(&mut self, k: usize, n: usize, K: usize, N: usize) -> f64 {
        // (K - x)/(x + 1) * (n - x)/(N - K - n + x + 1)
        let a = (K - k) as f64 / (k + 1) as f64;
        let b = (n - k) as f64 / (N - K - n + k + 1) as f64;
        a * b
    }

    /// Backward term ratio r_back(k) = t(k-1)/t(k)
    /// where t(k) = P(X=k) for Hypergeometric(n, K, N)
    #[allow(non_snake_case)]
    fn ratio_backward(&mut self, k: usize, n: usize, K: usize, N: usize) -> f64 {
        // x/(K - x + 1) * (N - K - n + x)/(n - x + 1)
        let a = k as f64 / (K - k + 1) as f64;
        let b = (N - K - n + k) as f64 / (n - k + 1) as f64;
        a * b
    }

    /// phyper - Cumulative density function of the hypergeometric distribution.
    ///
    /// Parameters
    /// - `k`: number of successes draws
    /// - `n`: number of draws
    /// - `K`: total number of success items
    /// - `N`: total number of items
    /// - `lower_tail` if true, then P(X <= x) is calculated, otherwise P(X > x) is calculated.
    ///
    /// Returns `P(X <= k)` or `P(X > k)`, the probability of observing up to (or less than) `k` successes in `n` draws.
    #[allow(non_snake_case)]
    pub fn phyper(
        &mut self,
        k: usize,
        n: usize,
        K: usize,
        N: usize,
        lower_tail: bool,
    ) -> Result<f64, String> {
        // Domain checks
        if n > N || K > N {
            return Ok(0.0);
        }
        let (k_min, k_max) = self.support(n, K, N);
        if k < k_min {
            return Ok(0.0);
        }
        if k > k_max {
            return Ok(1.0);
        }

        let mut lower_sum = 0f64;
        let mut upper_sum = 0f64;
        let mut c = 0f64;
        let mut y;
        let mut t;

        // let lower_len = k - lo + 1;
        // let upper_len = up - k;

        if lower_tail {
            for x in k_min..=k { // includes k
                let tx = self.dhyper(x, n, K, N)?;
                y = tx - c;
                t = lower_sum + y;
                c = (t - lower_sum) - y;
                lower_sum = t;
            }
            Ok(lower_sum)
        } else {
            for x in k + 1..=k_max { // excludes k
                let tx = self.dhyper(x, n, K, N)?;
                y = tx - c;
                t = upper_sum + y;
                c = (t - upper_sum) - y;
                upper_sum = t;
            }
            Ok(upper_sum)
        }

        /*
        let lower_len = k - lo + 1;
        let upper_len = up - k;
        // start at t(k) = P(X=k)
        let mut tk = self.dhyper(k, n, K, N)?;
        let mut k = k;
        if lower_len <= upper_len {
            // Sum lower tail: t(k) + t(k-1) + ... + t(lo)
            let mut lower_sum = tk;
            let mut c = 0;
            let mut y = tk;
            while k > lo {
                // t(k-1) = t(k) * r_back(k)
                tk = self.ratio_backward(k, n, K, N) * tk;
                lower_sum += tk;
                k -= 1;
            }
            if lower_tail {
                Ok(lower_sum)
            } else {
                Ok(1.0 - lower_sum)
            }
        } else {
            // sum upper tail: t(k+1) + ... + t(up)
            let mut upper_sum = 0f64;
            let mut k = k;
            while k < up {
                tk = self.ratio_forward(k, n, K, N) * tk;
                upper_sum += tk;
                k += 1;
            }
            if lower_tail {
                Ok(1.0 - upper_sum)
            } else {
                Ok(upper_sum)
            }
        }*/
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use float_eq::float_eq;
    use rstest::rstest;

    #[test]
    fn test_lfactorial() {
        let tests: Vec<(usize, i64)> = vec![
            (1, 1),
            (2, 2),
            (3, 6),
            (4, 24),
            (5, 120),
            (6, 720),
            (7, 5040),
            (8, 40_320),
            (9, 362_880),
            (10, 3_628_800),
            (11, 39_916_800),
            (12, 479_001_600),
            (13, 6_227_020_800),
            (14, 87_178_291_200),
            (15, 1_307_674_368_000),
        ];
        let mut hgeom = Hypergeometric::new();
        for test in tests {
            let lf = hgeom.log_factorial(test.0);
            assert!(lf.is_ok());
            let fact = lf.unwrap().exp();
            assert!(float_eq!(test.1 as f64, fact, rmax <= 1e-6));
        }
    }

    #[test]
    fn test_l_n_choose_k() {
        // Test log N-choose-K for K=4, N=20
        let mut hgeom = Hypergeometric::new();
        let lf4 = hgeom.log_factorial(4).unwrap();
        let lf16 = hgeom.log_factorial(16).unwrap();
        let lf20 = hgeom.log_factorial(20).unwrap();
        // lfactorial(20) =>  42.33562
        assert!(float_eq!(42.33562, lf20, rmax <= 1e-6));
        // lfactorial(16) =>  30.67186
        assert!(float_eq!(30.67186, lf16, rmax <= 1e-6));
        // lfactorial(4) =>  3.178054
        assert!(float_eq!(3.178054, lf4, rmax <= 1e-6));
        let nck = lf20 - lf4 - lf16; // calculate by hand
        // lchoose(20, 4) =>  8.485703 in R
        let expected_r = 8.485703;
        assert!(float_eq!(expected_r, nck, rmax <= 1e-6));
        let mylck = hgeom.log_n_choose_k(20, 4).unwrap();
        assert!(float_eq!(expected_r, mylck, rmax <= 1e-6));
    }

    #[test]
    fn test_dhyper() {
        let mut hgeom = Hypergeometric::new();
        // Let's first valid log-NchooseK
        // in R,  dhyper(4,20,45,10) yields 0.2204457
        // This means we choose 4 white balls from an urn with 20 white and 45 black balls when we take a total of 10 balls
        let a = hgeom.log_n_choose_k(20, 4);
        let b = hgeom.log_n_choose_k(45, 6);
        let c = hgeom.log_n_choose_k(65, 10);
        let lck = a.unwrap() + b.unwrap() - c.unwrap();
        // Exponentiate the result to get the result
        let n_choose_k = lck.exp();
        let expected = 0.2204457;
        assert!(float_eq!(expected, n_choose_k, rmax <= 1e-6));
        // Now test our version
        let result = hgeom.dhyper(4, 10, 20, 65);
        assert!(result.is_ok());
        let our_n_choose_k = result.unwrap();
        assert!(float_eq!(expected, our_n_choose_k, rmax <= 1e-6));
    }

    #[test]
    fn test_dhyper_edge() {
        // We have ten balls, all white, and we draw 10. All ten must be white
        let mut hgeom = Hypergeometric::new();
        let result = hgeom.dhyper(10, 10, 10, 10);
        assert!(result.is_ok());
        let our_n_choose_k = result.unwrap();
        assert!(float_eq!(1f64, our_n_choose_k, rmax <= 1e-6));
    }
    /// Tests Hypergeometric::dhyper against reference values from R's dhyper()
    #[rstest]
    #[case(3, 5, 4, 18, 0.04248366)] // valid case
    #[case(4, 10, 20, 65, 0.2204457)]
    #[case(5, 7, 6, 46, 8.74363e-05)]
    #[case(10, 10, 10, 10, 1.0)] // edge case: only white balls in urn, all drawn → P = 1
    #[case(0, 10, 0, 10, 1.0)] // edge case: only black balls in urn, all drawn → P = 1
    #[case(4, 5, 3, 17, 0.0)] // invalid case: x > m (drawing more white balls than available)
    #[case(6, 5, 8, 22, 0.0)] // invalid case: x > k (drawing more white balls than total draws)
    #[case(0, 5, 4, 6, 0.0)] // invalid case: k - x > n (drawing more black balls than available)
    #[allow(non_snake_case)]
    fn test_dhyper_cases(
        #[case] k: usize,
        #[case] n: usize,
        #[case] K: usize,
        #[case] N: usize,
        #[case] expected: f64,
    ) {
        let mut hgeom = Hypergeometric::new();
        let result = hgeom.dhyper(k, n, K, N).unwrap();
        assert!(
            float_eq!(result, expected, rmax <= 1e-6),
            "Expected {}, got {:.15}",
            expected,
            result
        );
    }

    #[test]
    fn test_phyper() {
        let mut hgeom = Hypergeometric::new();
        let result = hgeom.phyper(3 - 1, 190, 4, 1526, false);
        // assertTrue(result > 0.0069 && result < 0.0070); -- from ontologizer code.
        println!("{}", result.unwrap());

        let result = hgeom.phyper(66 - 1, 262, 599, 8376, false);
        // assertTrue(result > 0.0069 && result < 0.0070); -- from ontologizer code.
        println!("{}", result.unwrap());

        let result = hgeom.phyper(72 - 1, 262, 899, 8376, false);
        // assertTrue(result > 0.0069 && result < 0.0070); -- from ontologizer code.
        println!("{}", result.unwrap());
    }
}
