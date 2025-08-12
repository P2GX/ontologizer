use std::cmp::min;

pub struct Hypergeometric {
    // log factorial
    lfactorial: Vec<f64>,
}

impl Hypergeometric {
    pub fn new() -> Self {
        let mut lfac = vec![];
        // by convention, we let log(0) = 0, because 0! = 1
        lfac.push(0 as f64); // log fac(0)
        lfac.push(0 as f64); // log fac(1)
        Hypergeometric { lfactorial: lfac }
    }

    pub fn logfact(&mut self, i: usize) -> Result<f64, String> {
        /*
         * Make sure value is already in lfactorial. If not, calculate all
         * values up to that for i
         */
        if i >= self.lfactorial.len() {
            for j in self.lfactorial.len()..=i {
                let last_logfact = self.lfactorial.last().unwrap();
                self.lfactorial.push(last_logfact + (j as f64).ln());
            }
        }
        self.lfactorial
            .get(i)
            .copied()
            .ok_or_else(|| format!("Could not calculate log factorial for j={}", i))
    }

    /**
     * <P>
     * For the hypergeometric distribution note the following.
     * </P>
     * <UL>
     * <LI>We set up the problem as sampling a set of genes (study genes, for
     * instance, the set of upregulated genes in some microarray experiment)
     * from a larger set of genes (say, the set of all genes of a species). The
     * sampling is done without replacement. </LI>
     * <LI>For each GO term, we can conceive of the population as being divided
     * into genes annotated to this term and those genes that are not annotated
     * to the term.</LI>
     * <LI>The probability of having a certain number of terms annotated to the
     * term in the study set can then be calculated by the hypergeometric
     * distribution. See GeneMerge by Castillo-Davis et al (Bioinformatics).
     * </LI>
     * <LI>To do this, we need to divide the population into two groups, genes
     * with and without annotation. The arguments to the function supply us with
     * <B>n</B>, the total number of genes in the population group, and <B>p
     * </B>, the proportion of genes with annotation to the term in question.
     * We can then calculate the number of genes in the population annotated to
     * the term by <B>round(n*p)</B>, and the number of genes not annotated to
     * the term by <B>round(n*(1-p))</B>.</LI>
     * </UL>
     *
     *
     * @param n
     *            Number of population genes
     * @param p
     *            Proportion of population genes
     * @param k
     *            Number of study genes
     * @param r
     *            Number of study genes in group
     */
    pub fn phypergeometric(&mut self, n: usize, p: f64, k: usize, r: usize) -> Result<f64, String> {
        /*
         * Study group cannot be larger than population. If this happens there
         * is probably something wrong with the input data, but returning 1.0
         * prevents confusing and wrong output.
         */
        if k >= n {
            return Ok(1 as f64);
        }

        if r < 1 {
            return Ok(1 as f64); // Not valid for r < 2, less than 2 study genes.
        }

        let q = 1.0 - p;
        let np = (n as f64 * p).round() as usize; // Round to nearest int
        let nq = (n as f64 * q).round() as usize;

        let log_n_choose_k = self.lNchooseK(n, k)?;
        let mut top = k;
        if np < k {
            top = np;
        }

        let mut lfoo = self.lNchooseK(np, top)? + self.lNchooseK(nq, k - top)?;

        let mut sum: f64 = 0.0;

        for i in (r..=top).rev() {
            sum += (lfoo - log_n_choose_k).exp();
            if i > r {
                lfoo = lfoo
                    + (i as f64) / (np - i + 1) as f64
                    + (nq - k + i) as f64 / (k - i + 1) as f64
                //+ java.lang.Math.log((double) i / (double) (np - i + 1))
                //+ java.lang.Math.log((double) (nq - k + i)
                //		/ (double) (k - i + 1));
            }
        }
        return Ok(sum);
    }

    /// n-choose-k in log space
    pub fn lNchooseK(&mut self, n: usize, k: usize) -> Result<f64, String> {
        let result = self.logfact(n)? - self.logfact(k)? - self.logfact(n - k)?;
        Ok(result)
    }

    /// dhyper: Density for the hypergeometric function
    ///
    /// # Parameters
    /// - `x`: number of white balls drawn without replacement
    /// - `m``: the number of white balls in the urn.
    /// - `n``: the number of black balls in the urn.
    /// - `k``: the number of balls drawn from the urn, hence must be in 0,1,..., m+n.
    ///
    /// # Returns
    /// - Returns the density function as a Result.
    pub fn dhyper(&mut self, x: usize, m: usize, n: usize, k: usize) -> Result<f64, String> {
        // It is not possible to draw more white balls from an urn containing M white balls. Hence the probability is 0.
        if x > m {
            return Ok(0 as f64);
        }
        // it is also not possible to draw more white balls than the number of draws. Probability is 0
        if x > k {
            return Ok(0 as f64);
        }
        // it is also not possible to draw more black balls (k-x) than there are within the urn.
        if k - x > n {
            return Ok(0 as f64);
        }
        // #way to choose x white balls + #ways to choose k-x black balls less #ways to choose k balls from total (m+n)
        let result =
            self.lNchooseK(m, x)? + self.lNchooseK(n, k - x)? - self.lNchooseK(m + n, k)?;
        Ok(result.exp())
    }

    /// phyper - Calculates P(X &gt; x) where X is the hypergeometric distribution with indices N,M,n.
    ///
    /// If lower_tail is set to true, then P(X &lt;= x) is calculated.
    ///
    /// # Parameters
    /// - `x`:  x number of white balls drawn without replacement
    /// - `N number of balls in the urn
    /// - `M number of white balls in the urn
    /// - `n number of balls drawn from the urn
    /// - `lower_tail`` if true, then P(X &lt;= x) is calculated, otherwise P(X &gt; x) is calculated.
    ///
    /// # Returns
    /// - Returns the density function as a Result.
    pub fn phyper(
        &mut self,
        x: usize,
        N: usize,
        M: usize,
        n: usize,
        lower_tail: bool,
    ) -> Result<f64, String> {
        let mut up;
        let mut p = 0 as f64;

        up = min(n, M);

        if x < up / 2 {
            for i in (0..=x).rev() {
                p += self.dhyper(i, M, N - M, n)?;
            }

            if lower_tail {
                return Ok(p);
            } else {
                return Ok(1 as f64 - p);
            }
        } else {
            for i in x + 1..=up {
                p += self.dhyper(i, M, N - M, n)?;
            }
            if lower_tail {
                return Ok(1 as f64 - p);
            } else {
                return Ok(p);
            }
        }
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
            let lf = hgeom.logfact(test.0);
            assert!(lf.is_ok());
            let fact = lf.unwrap().exp();
            assert!(float_eq!(test.1 as f64, fact, rmax <= 1e-6));
        }
    }

    #[test]
    fn test_l_n_choose_k() {
        // Test log N-choose-K for K=4, N=20
        let mut hgeom = Hypergeometric::new();
        let lf4 = hgeom.logfact(4).unwrap();
        let lf16 = hgeom.logfact(16).unwrap();
        let lf20 = hgeom.logfact(20).unwrap();
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
        let mylck = hgeom.lNchooseK(20, 4).unwrap();
        assert!(float_eq!(expected_r, mylck, rmax <= 1e-6));
    }

    #[test]
    fn test_dhyper() {
        let mut hgeom = Hypergeometric::new();
        // Let's first valid log-NchooseK
        // in R,  dhyper(4,20,45,10) yields 0.2204457
        // This means we choose 4 white balls from an urn with 20 white and 45 black balls when we take a total of 10 balls
        let mut a = hgeom.lNchooseK(20, 4);
        let b = hgeom.lNchooseK(45, 6);
        let c = hgeom.lNchooseK(65, 10);
        let lck = a.unwrap() + b.unwrap() - c.unwrap();
        // Exponentiate the result to get the result
        let n_choose_k = lck.exp();
        let expected = 0.2204457;
        assert!(float_eq!(expected, n_choose_k, rmax <= 1e-6));
        // Now test our version
        let result = hgeom.dhyper(4, 20, 45, 10);
        assert!(result.is_ok());
        let our_n_choose_k = result.unwrap();
        assert!(float_eq!(expected, our_n_choose_k, rmax <= 1e-6));
    }

    #[test]
    fn test_dhyper_edge() {
        // We have ten balls, all white, and we draw 10. All ten must be white
        let mut hgeom = Hypergeometric::new();
        let result = hgeom.dhyper(10, 10, 0, 10);
        assert!(result.is_ok());
        let our_n_choose_k = result.unwrap();
        assert!(float_eq!(1 as f64, our_n_choose_k, rmax <= 1e-6));
    }
    /// Tests Hypergeometric::dhyper against reference values from R's dhyper()
    #[rstest]
    #[case(3, 4, 14, 5, 0.04248366)] // valid case
    #[case(4, 20, 45, 10, 0.2204457)]
    #[case(5, 6, 40, 7, 8.74363e-05)]
    #[case(10, 10, 0, 10, 1.0)] // edge case: only white balls in urn, all drawn → P = 1
    #[case(0, 0, 10, 10, 1.0)] // edge case: only black balls in urn, all drawn → P = 1
    #[case(4, 3, 14, 5, 0.0)] // invalid case: x > m (drawing more white balls than available)
    #[case(6, 8, 14, 5, 0.0)] // invalid case: x > k (drawing more white balls than total draws)
    #[case(0, 4, 2, 5, 0.0)] // invalid case: k - x > n (drawing more black balls than available)
    fn test_dhyper_cases(
        #[case] x: usize,
        #[case] m: usize,
        #[case] n: usize,
        #[case] k: usize,
        #[case] expected: f64,
    ) {
        let mut hgeom = Hypergeometric::new();
        let result = hgeom.dhyper(x, m, n, k).unwrap();
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
        let result = hgeom.phyper(2, 1526, 4, 190, false);
        // assertTrue(result > 0.0069 && result < 0.0070); -- from ontologizer code.
        println!("{}", result.unwrap());
    }
}
