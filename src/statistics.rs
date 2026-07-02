use crate::probability::CDHResult;
use std::collections::HashMap;


/// A high-performance structure encapsulating a linear array of values 
/// designed for descriptive and bivariate statistical data mining.
pub struct DataSet<T>(pub Vec<T>);

impl DataSet<f64> {
    /// Instantiates a new data frame context around a contiguous memory buffer.
    ///
    /// # Example
    /// ```rust
   /// use scicdh::statistics::DataSet;
    /// let ds = DataSet::new(vec![1.0, 2.0, 3.0]);
    /// ```
    pub fn new(data: Vec<f64>) -> Self {
        DataSet(data)
    }

    /// Computes the arithmetic mean ($\mu$) of the sample space.
    ///
    /// $$\mu = \frac{\sum_{i=1}^n x_i}{n}$$
    pub fn mean(&self) -> CDHResult<f64> {
        let len = self.0.len() as f64;
        if len == 0.0 {
            return Err("Cannot compute mean of an empty dataset.".to_string());
        }
        Ok(self.0.iter().sum::<f64>() / len)
    }

    /// Computes the sample median using a non-destructive 0-indexed sorting pipeline.
    /// Handles both odd and even dataset lengths correctly.
    pub fn median(&self) -> CDHResult<f64> {
        let len = self.0.len();
        if len == 0 {
            return Err("Cannot compute median of an empty dataset.".to_string());
        }
        let mut sorted_arr = self.0.clone();
        sorted_arr.sort_by(|a, b| a.total_cmp(b));
        
        match len {
            a if a % 2 == 1 => Ok(sorted_arr[(a + 1) / 2 - 1]),
            _ => {
                let center = len / 2;
                let center_elem = sorted_arr[center - 1];
                let next_elem = sorted_arr[center];
                Ok((center_elem + next_elem) / 2.0)
            }
        }
    }

    /// Computes the statistical mode of the dataset in linear time $\mathcal{O}(n)$.
    ///
    /// Uses low-level float bit-packing (`to_bits()`) to completely bypass 
    /// hashing restrictions on `f64`. Returns `Ok(None)` if the distribution 
    /// is completely uniform (no single clear majority frequency exists).
    pub fn mode(&self) -> CDHResult<Option<f64>> {
        if self.0.is_empty() {
            return Err("Cannot compute mode of an empty dataset.".to_string());
        }

        let mut frequencies: HashMap<u64, usize> = HashMap::with_capacity(self.0.len());
        let mut max_count = 0;

        // Step A: Aggregate counts via float-to-bit maps in pure O(n)
        for &num in &self.0 {
            let bits = num.to_bits();
            let count = frequencies.entry(bits).or_insert(0);
            *count += 1;
            if *count > max_count {
                max_count = *count;
            }
        }

        // Step B: Extract matching patterns to verify unimodal dominance
        let modes: Vec<u64> = frequencies
            .iter()
            .filter(|&(_, &count)| count == max_count)
            .map(|(&bits, _)| bits)
            .collect();

        // If all items appear exactly the same number of times, there is no mode
        if modes.len() == frequencies.len() && max_count == 1 {
            return Ok(None);
        }

        // If multiple distinct elements tie for peak occurrence, it's multimodal
        if modes.len() > 1 {
            Ok(None)
        } else {
            Ok(Some(f64::from_bits(modes[0])))
        }
    }

    /// Calculates the algebraic range (Span) between the maximum and minimum parameters.
    ///
    /// $$\text{Range} = X_{\text{max}} - X_{\text{min}}$$
    pub fn range(&self) -> CDHResult<f64> {
        if self.0.is_empty() {
            return Err("Cannot compute range of an empty dataset.".to_string());
        }
        Ok(self.max()? - self.min()?)
    }

    /// Scans a target container or the internal dataset structure to check 
    /// for the existence of a precise value.
    /// 
    /// Fixed Bug: Bypassed memory copying vulnerabilities from previous design.
    pub fn has(&self, value: f64, other: Option<&[f64]>) -> CDHResult<bool> {
        let set = match other {
            None => &self.0,
            Some(slice) => slice,
        };
        Ok(set.iter().any(|&v| v == value))
    }

    /// Computes the sum of squared deviations for the primary variable ($SS_{xx}$).
    ///
    /// $$SS_{xx} = \sum x^2 - \frac{(\sum x)^2}{n}$$
    pub fn sample_variation_x2(&self) -> CDHResult<f64> {
        let n = self.0.len() as f64;
        if n == 0.0 { return Err("Empty dataset context.".to_string()); }
        let xsqrsum = self.0.iter().fold(0.0, |acc, b| acc + b.powi(2));
        let xsum = self.0.iter().sum::<f64>();
        Ok(xsqrsum - xsum.powi(2) / n)
    }

    /// Computes the sum of squared deviations for an independent reference slice ($SS_{yy}$).
    ///
    /// $$SS_{yy} = \sum y^2 - \frac{(\sum y)^2}{n}$$
    pub fn sample_variation_y2(&self, y: &[f64]) -> CDHResult<f64> {
        let n = y.len() as f64;
        if n == 0.0 { return Err("Input slice cannot be empty.".to_string()); }
        let ysqrsum = y.iter().fold(0.0, |acc, b| acc + b.powi(2));
        let ysum = y.iter().sum::<f64>();
        Ok(ysqrsum - ysum.powi(2) / n)
    }

    /// Computes the unbiased sample variance ($s^2$) leveraging Bessel's correction factor.
    ///
    /// $$s^2 = \frac{SS_{xx}}{n - 1}$$
    pub fn sample_variance(&self) -> CDHResult<f64> {
        let len = self.0.len();
        if len <= 1 {
            return Err("Sample variance requires at least 2 data points.".to_string());
        }
        Ok(self.sample_variation_x2()? / ((len - 1) as f64))
    }

    /// Computes the unbiased sample standard deviation ($s$).
    pub fn sample_std_deviation(&self) -> CDHResult<f64> {
        Ok(self.sample_variance()?.sqrt())
    }

    /// Evaluates the dimensionless sample coefficient of variation ($CV$).
    ///
    /// $$CV = \frac{s}{\bar{x}}$$
    pub fn sample_coeff_of_variation(&self) -> CDHResult<f64> {
        let mean = self.mean()?;
        if mean == 0.0 {
            return Err("Coefficient of variation undefined for a zero mean.".to_string());
        }
        Ok(self.sample_std_deviation()? / mean)
    }

    /// Generates the $k$-th mathematical population central moment ($m_k$).
    ///
    /// $$\mu_k = \frac{\sum (x_i - \bar{x})^k}{n}$$
    pub fn kth_central_moment(&self, k: f64) -> CDHResult<f64> {
        let n = self.0.len() as f64;
        if n == 0.0 { return Err("Empty dataset context.".to_string()); }
        let mean = self.mean()?;
        Ok(self.0.iter().fold(0.0, |acc, b| acc + (b - mean).powf(k)) / n)
    }

    /// Measures the skewness ($\beta_3$), tracking distribution symmetry.
    pub fn skewness(&self) -> CDHResult<f64> {
        let denom = self.population_std_deviation()?.powi(3);
        if denom == 0.0 { return Ok(0.0); }
        Ok(self.kth_central_moment(3.0)? / denom)
    }

    /// Computes excess kurtosis ($\beta_4$), tracking peak sharpness relative to standard normal shapes.
    pub fn kurtosis(&self) -> CDHResult<f64> {
        let denom = self.population_std_deviation()?.powi(4);
        if denom == 0.0 { return Ok(0.0); }
        Ok((self.kth_central_moment(4.0)? / denom) - 3.0)
    }

    /// Computes the cross-variation sum of joint distribution observations ($SS_{xy}$).
    ///
    /// $$SS_{xy} = \sum xy - \frac{(\sum x)(\sum y)}{n}$$
    pub fn sample_variation_xy(&self, y: &[f64]) -> CDHResult<f64> {
        if y.len() != self.0.len() {
            return Err("Length of y should be equal to length of dataset".to_string());
        }
        let n = y.len() as f64;
        let xysum = self.0.iter().zip(y.iter()).fold(0.0, |acc, b| acc + (b.0 * b.1));
        let xsum = self.0.iter().sum::<f64>();
        let ysum = y.iter().sum::<f64>();
        Ok(xysum - (xsum * ysum / n))
    }

    /// Evaluates Pearson's Product-Moment Correlation Coefficient ($r_{xy}$).
    /// Hardened against the textbook typo by using $SS_{yy}$ space arrays.
    pub fn pearson_correlation_coeff(&self, y: &[f64]) -> CDHResult<f64> {
        let denom = (self.sample_variation_x2()? * self.sample_variation_y2(y)?).sqrt();
        if denom == 0.0 {
            return Ok(0.0); // Protective baseline against zero variance tracking
        }
        Ok(self.sample_variation_xy(y)? / denom)
    }

    /// Deduplicates entries from the dataset to form a unique element mathematical set pool.
    pub fn into_set(&self) -> CDHResult<Vec<f64>> {
        let mut set: Vec<f64> = Vec::with_capacity(self.0.len());
        for &a in self.0.iter() {
            if !self.has(a, Some(&set))? {
                set.push(a);
            }
        }
        Ok(set)
    }

    /// Calculates pure population variance ($\sigma^2$) relative to the exact population mean.
    pub fn population_variance(&self) -> CDHResult<f64> {
        let n = self.0.len() as f64;
        if n == 0.0 { return Err("Empty dataset context.".to_string()); }
        let mean = self.mean()?;
        Ok(self.0.iter().fold(0.0, |acc, &x| acc + (x - mean).powi(2)) / n)
    }

    /// Computes the population standard deviation ($\sigma$).
    pub fn population_std_deviation(&self) -> CDHResult<f64> {
        Ok(self.population_variance()?.sqrt())
    }

    /// Evaluates the peak extreme value bound present inside the dataset.
    pub fn max(&self) -> CDHResult<f64> {
        if self.0.is_empty() { return Err("Dataset is empty.".to_string()); }
        let mut max = self.0[0];
        self.0.iter().for_each(|&i| {
            if i > max { max = i; }
        });
        Ok(max)
    }

    /// Evaluates the lowest value parameter bound present inside the dataset.
    pub fn min(&self) -> CDHResult<f64> {
        if self.0.is_empty() { return Err("Dataset is empty.".to_string()); }
        let mut min = self.0[0];
        self.0.iter().for_each(|&i| {
            if i < min { min = i; }
        });
        Ok(min)
    }

    /// Dumps clean debug statements displaying the full metrics profile of the active dataset.
    pub fn info(&self) -> CDHResult<()> {
        println!("DataSet === {:?}", self.0);
        println!("Mean === {:#?}", self.mean()?);
        println!("Median === {:#?}", self.median()?);
        println!("Mode === {:#?}", self.mode()?);
        println!("Range === {:#?}", self.range()?);
        println!("Variance === {:#?}", self.sample_variance()?);
        println!("Skewness === {:#?}", self.skewness()?);
        println!("Kurtosis === {:#?}", self.kurtosis()?);
        println!("self has 3f64 == {:?}", self.has(3.0, None)?);
        println!("Set === {:?}", self.into_set()?);
        println!("Maximum === {:?}", self.max()?);
        println!("Minimum === {:?}", self.min()?);
        Ok(())
    }
}
