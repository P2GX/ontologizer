use std::cmp::{max, min};

pub trait DiscreteDistribution {
    /// Probability Mass Function: P(X = k)
    fn pmf(&self, k: usize) -> f64;

    /// Cumulative Distribution Function: P(X <= k)
    #[allow(unused)]
    fn cdf(&self, k: usize) -> f64;

    /// Survival Function: P(X > k)
    fn sf(&self, k: usize) -> f64;
}

pub struct LogFactorialCache {
    log_factorials: Vec<f64>,
}

impl LogFactorialCache {
    pub fn new(max_n: usize) -> Self {
        let mut log_facts = Vec::with_capacity(max_n + 1);
        log_facts.push(0.0); // log(0!)

        let mut current = 0.0;
        for i in 1..=max_n {
            current += (i as f64).ln();
            log_facts.push(current);
        }

        Self {
            log_factorials: log_facts,
        }
    }

    pub fn log_factorial(&self, n: usize) -> f64 {
        // Ensure log_factorial[i] is computed; if not, extend the table.
        if n < self.log_factorials.len() {
            self.log_factorials[n]
        } else {
            let mut res = *self.log_factorials.last().unwrap();
            for i in self.log_factorials.len()..=n {
                res += (i as f64).ln();
            }
            res
        }
    }

    pub fn log_n_choose_k(&self, n: usize, k: usize) -> f64 {
        self.log_factorial(n) - self.log_factorial(k) - self.log_factorial(n - k)
    }
}

#[allow(non_snake_case)]
pub struct Hypergeometric<'a> {
    n: usize, // draws (study size)
    K: usize, // total successes (annotation count)
    N: usize, // total items (population size)
    cache: &'a LogFactorialCache,
}

#[allow(non_snake_case, unused)]
impl<'a> Hypergeometric<'a> {
    pub fn new(n: usize, K: usize, N: usize, cache: &'a LogFactorialCache) -> Self {
        Self { n, K, N, cache }
    }

    pub fn n(&self) -> usize {
        self.n
    }
    pub fn K(&self) -> usize {
        self.K
    }
    pub fn N(&self) -> usize {
        self.N
    }
    fn support(&self) -> (usize, usize) {
        // Compute the support of the distributio, i.e. the range of values for which P(k) != 0.
        let k_min = max(
            0isize,
            (self.n as isize) + (self.K as isize) - (self.N as isize),
        ) as usize;
        // Example: n=100, N=200, K=120 -> cannot draw less than 20 successes
        let k_max = min(self.n, self.K);
        // Example: n=100, K=80, -> cannot draw more than 80 successes
        (k_min, k_max)
    }
}
impl<'a> DiscreteDistribution for Hypergeometric<'a> {
    /// pdf: Probability mass function of the hypergeometric distribution.
    ///
    /// Parameters
    /// - `k`: number of successes draws
    /// Returns `P(X = k)`, the probability of observing exactly `k` successes in `n` draws.
    #[allow(non_snake_case)]
    fn pmf(&self, k: usize) -> f64 {
        if self.n > self.N || self.K > self.N {
            return 0.0;
        }

        let (k_min, k_max) = self.support();
        if k < k_min || k > k_max {
            return 0.0;
        }

        let log_prob = self.cache.log_n_choose_k(self.K, k)
            + self.cache.log_n_choose_k(self.N - self.K, self.n - k)
            - self.cache.log_n_choose_k(self.N, self.n);
        log_prob.exp()
    }

    /// cdf - Cumulative distribution function of the hypergeometric distribution.
    ///
    /// Parameters
    /// - `k`: number of successes draws
    /// - `lower_tail` if true, then P(X <= x) is calculated, otherwise P(X > x) is calculated.
    ///
    /// Returns `P(X <= k)` or `P(X > k)`, the probability of observing up to (or less than) `k` successes in `n` draws.
    #[allow(non_snake_case)]
    fn cdf(&self, k: usize) -> f64 {
        // Domain checks
        if self.n > self.N || self.K > self.N {
            return 0.0;
        }

        let (k_min, k_max) = self.support();
        if k < k_min {
            return 0.0;
        }
        if k >= k_max {
            return 1.0;
        }

        let mut sum = 0f64;
        let mut c = 0f64;
        let mut y;
        let mut t;

        for x in k_min..=k {
            // includes k
            let tx = self.pmf(x);
            y = tx - c;
            t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        sum
    }

    fn sf(&self, k: usize) -> f64 {
        // Domain checks
        if self.n > self.N || self.K > self.N {
            return 0.0;
        }

        let (k_min, k_max) = self.support();
        if k < k_min {
            return 1.0;
        }
        if k >= k_max {
            return 0.0;
        }

        let mut sum = 0f64;
        let mut c = 0f64;
        let mut y;
        let mut t;

        for x in k + 1..=k_max {
            // excludes k
            let tx = self.pmf(x);
            y = tx - c;
            t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        sum
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use float_eq::float_eq;
    use rstest::rstest;

    #[test]
    fn test_log_factorials() {
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
        let cache = LogFactorialCache::new(15);
        for test in tests {
            let lf = cache.log_factorial(test.0);
            let fact = lf.exp();
            assert!(float_eq!(test.1 as f64, fact, rmax <= 1e-6));
        }
    }

    #[test]
    fn test_log_n_choose_k() {
        // Test log N-choose-K for K=4, N=20
        let cache = LogFactorialCache::new(25);
        let lf4 = cache.log_factorial(4);
        let lf16 = cache.log_factorial(16);
        let lf20 = cache.log_factorial(20);
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
        let mylck = cache.log_n_choose_k(20, 4);
        assert!(float_eq!(expected_r, mylck, rmax <= 1e-6));
    }

    #[test]
    fn test_pmf() {
        let cache = LogFactorialCache::new(100);
        // Let's first valid log-NchooseK
        // in R,  dhyper(4,20,45,10) yields 0.2204457
        // This means we choose 4 white balls from an urn with 20 white and 45 black balls when we take a total of 10 balls
        let a = cache.log_n_choose_k(20, 4);
        let b = cache.log_n_choose_k(45, 6);
        let c = cache.log_n_choose_k(65, 10);
        let lck = a + b - c;
        // Exponentiate the result to get the result
        let n_choose_k = lck.exp();
        let expected = 0.2204457;
        assert!(float_eq!(expected, n_choose_k, rmax <= 1e-6));
        // Now test our version
        let hypergeo = Hypergeometric::new(10, 20, 65, &cache);
        // let result = hgeom.pmf(4, 10, 20, 65);
        let n_choose_k = hypergeo.pmf(4);
        assert!(float_eq!(expected, n_choose_k, rmax <= 1e-6));
    }

    #[test]
    fn test_pmf_edge() {
        // We have ten balls, all white, and we draw 10. All ten must be white
        let cache = LogFactorialCache::new(10);
        let hypergeometric = Hypergeometric::new(10, 10, 10, &cache);
        let n_choose_k = hypergeometric.pmf(10);
        assert!(float_eq!(1f64, n_choose_k, rmax <= 1e-6));
    }

    // Tests Hypergeometric::dhyper against reference values from R's dhyper()
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
    fn test_pmf_cases(
        #[case] k: usize,
        #[case] n: usize,
        #[case] K: usize,
        #[case] N: usize,
        #[case] expected: f64,
    ) {
        let cache = LogFactorialCache::new(N);
        let hypergeometric = Hypergeometric::new(n, K, N, &cache);
        let result = hypergeometric.pmf(k);
        assert!(
            float_eq!(result, expected, rmax <= 1e-6),
            "Expected {}, got {:.15}",
            expected,
            result
        );
    }

    #[test]
    fn test_sf() {
        let cache = LogFactorialCache::new(10_000);

        let n = 190;
        let K = 4;
        let N = 1526;
        let hypergeometric = Hypergeometric::new(n, K, N, &cache);
        let result = hypergeometric.sf(3 - 1);
        // assertTrue(result > 0.0069 && result < 0.0070); -- from ontologizer code.
        println!("{}", result);
    }

    #[test]
    fn test_sf_large_sums() {
        let cache = LogFactorialCache::new(10_000);
        let n = 262;
        let K = 599;
        let N = 8376;
        let hypergeometric = Hypergeometric::new(n, K, N, &cache);
        let result = hypergeometric.sf(66 - 1);
        // assertTrue(result > 0.0)
        println!("{}", result);

        let n = 262;
        let K = 899;
        let N = 8376;
        let hypergeometric = Hypergeometric::new(n, K, N, &cache);
        let result = hypergeometric.sf(72 - 1);
        // assertTrue(result > 0.0069 && result < 0.0070); -- from ontologizer code.
        println!("{}", result);
    }
}
