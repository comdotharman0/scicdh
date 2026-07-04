use crate::probability::CDHResult;
use std::collections::HashMap;
use std::ops::{Add,Sub, Mul, Div, Index, IndexMut,
AddAssign, SubAssign, MulAssign, DivAssign};
use std::fmt;


/// A high-performance structure encapsulating a linear array of values 
/// designed for descriptive and bivariate statistical data mining.
#[derive(Debug, Clone)]
pub struct DataSet<T>{
pub data: Vec<T>
}



impl DataSet<f64> {
    /// Instantiates a new data frame context around a contiguous memory buffer.
    ///
    /// # Example
    /// ```rust
   /// use scicdh::statistics::DataSet;
    /// let ds = DataSet::new(vec![1.0, 2.0, 3.0]);
    /// ```
    pub fn new(data: Vec<f64>) -> Self {
        DataSet{data}
    }

    /// Computes the arithmetic mean ($\mu$) of the sample space.
    ///
    /// $$\mu = \frac{\sum_{i=1}^n x_i}{n}$$
    pub fn mean(&self) -> CDHResult<f64> {
        let len = self.data.len() as f64;
        if len == 0.0 {
            return Err("Cannot compute mean of an empty dataset.".to_string());
        }
        Ok(self.data.iter().sum::<f64>() / len)
    }

    /// Computes the sample median using a non-destructive 0-indexed sorting pipeline.
    /// Handles both odd and even dataset lengths correctly.
    pub fn median(&self) -> CDHResult<f64> {
        let len = self.data.len();
        if len == 0 {
            return Err("Cannot compute median of an empty dataset.".to_string());
        }
        let mut sorted_arr = self.data.clone();
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
        if self.data.is_empty() {
            return Err("Cannot compute mode of an empty dataset.".to_string());
        }

        let mut frequencies: HashMap<u64, usize> = HashMap::with_capacity(self.data.len());
        let mut max_count = 0;

        // Step A: Aggregate counts via float-to-bit maps in pure O(n)
        for &num in &self.data {
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
        if self.data.is_empty() {
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
            None => &self.data,
            Some(slice) => slice,
        };
        Ok(set.iter().any(|&v| v == value))
    }

    /// Computes the sum of squared deviations for the primary variable ($SS_{xx}$).
    ///
    /// $$SS_{xx} = \sum x^2 - \frac{(\sum x)^2}{n}$$
    pub fn sample_variation_x2(&self) -> CDHResult<f64> {
        let n = self.data.len() as f64;
        if n == 0.0 { return Err("Empty dataset context.".to_string()); }
        let xsqrsum = self.data.iter().fold(0.0, |acc, b| acc + b.powi(2));
        let xsum = self.data.iter().sum::<f64>();
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
        let len = self.data.len();
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
        let n = self.data.len() as f64;
        if n == 0.0 { return Err("Empty dataset context.".to_string()); }
        let mean = self.mean()?;
        Ok(self.data.iter().fold(0.0, |acc, b| acc + (b - mean).powf(k)) / n)
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
        if y.len() != self.data.len() {
            return Err("Length of y should be equal to length of dataset".to_string());
        }
        let n = y.len() as f64;
        let xysum = self.data.iter().zip(y.iter()).fold(0.0, |acc, b| acc + (b.0 * b.1));
        let xsum = self.data.iter().sum::<f64>();
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
        let mut set: Vec<f64> = Vec::with_capacity(self.data.len());
        for &a in self.data.iter() {
            if !self.has(a, Some(&set))? {
                set.push(a);
            }
        }
        Ok(set)
    }

    /// Calculates pure population variance ($\sigma^2$) relative to the exact population mean.
    pub fn population_variance(&self) -> CDHResult<f64> {
        let n = self.data.len() as f64;
        if n == 0.0 { return Err("Empty dataset context.".to_string()); }
        let mean = self.mean()?;
        Ok(self.data.iter().fold(0.0, |acc, &x| acc + (x - mean).powi(2)) / n)
    }

    /// Computes the population standard deviation ($\sigma$).
    pub fn population_std_deviation(&self) -> CDHResult<f64> {
        Ok(self.population_variance()?.sqrt())
    }

    /// Evaluates the peak extreme value bound present inside the dataset.
    pub fn max(&self) -> CDHResult<f64> {
        if self.data.is_empty() { return Err("Dataset is empty.".to_string()); }
        let mut max = self.data[0];
        self.data.iter().for_each(|&i| {
            if i > max { max = i; }
        });
        Ok(max)
    }

    /// Evaluates the lowest value parameter bound present inside the dataset.
    pub fn min(&self) -> CDHResult<f64> {
        if self.data.is_empty() { return Err("Dataset is empty.".to_string()); }
        let mut min = self.data[0];
        self.data.iter().for_each(|&i| {
            if i < min { min = i; }
        });
        Ok(min)
    }

    /// Dumps clean debug statements displaying the full metrics profile of the active dataset.
 pub fn info(&self) -> CDHResult<String> {
      
 Ok(   format!(
"DataSet === {:?}
Mean === {:#?} 
Median === {:#?}
Mode === {:#?} 
Range === {:#?}
Variance === {:#?}
Skewness === {:#?}
Kurtosis === {:#?}
Maximum === {:?} 
Minimum === {:?} 
",self.data,self.mean()?,self.median()?,
self.mode()?,self.range()?,self.sample_variance()?,
self.skewness()?,self.kurtosis()?,
self.max()?,self.min()?)
)


    }


  #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.data.iter()
    }

    /// Returns a mutable iterator over the data points.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, f64> {
        self.data.iter_mut()
    }

}



impl fmt::Display for DataSet<f64> {
 fn fmt(&self, f: &mut fmt::Formatter)
 -> fmt::Result {
write!(f,
 " \u{001b}[34m DataSet Info \n\
 {} \n\
 {} \n\
{} \u{001b}[0m \n ","=".repeat(45),
 self.info().unwrap(),
"=".repeat(45))

}
}


impl Add<DataSet<f64>> for DataSet<f64>{
    type Output = CDHResult<DataSet<f64>>;

    fn add(self, other: DataSet<f64>) 
-> <Self as Add<DataSet<f64>>>::Output

 {
let len1 = self.data.len();
let len2 = self.data.len();
if len1 != len2{
return Err("Both Datasets should be of equal length".to_string());
}
let data: Vec<f64> = self.data.iter()
.zip(other.data.iter())
.map(|b| b.0+b.1).collect();
Ok(DataSet{data})
}

}



impl Add<f64> for DataSet<f64> {
    type Output = CDHResult<DataSet<f64>>;

    fn add(self, rhs: f64) -> Self::Output {
let data = self.data.iter()
.map(|i| i+rhs).collect();
Ok(DataSet{data})
}

}






// 1. Immutable Indexing: allows reading data via dataset[index]
impl Index<usize> for DataSet<f64> {
    type Output = f64;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        // Enforces normal Rust boundary panics if index is out of bounds
        &self.data[index]
    }
}

// 2. Mutable Indexing: allows writing data via dataset[index] = value
impl IndexMut<usize> for DataSet<f64> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}


// =========================================================================
// 1. DATASET - DATASET (Vector Subtraction)
// =========================================================================
impl Sub<DataSet<f64>> for DataSet<f64> {
    type Output = CDHResult<DataSet<f64>>;

    fn sub(self, other: DataSet<f64>) -> Self::Output {
        let len1 = self.data.len();
        let len2 = other.data.len();
        
        if len1 != len2 {
            return Err("Both Datasets should be of equal length".to_string());
        }

        // Optimized zero-allocation zip operation
        let data: Vec<f64> = self.data.iter()
            .zip(other.data.iter())
            .map(|(a, b)| a - b)
            .collect();

        Ok(DataSet{data})
    }
}

// =========================================================================
// 2. DATASET - SCALAR (Shifting every element down)
// =========================================================================
impl Sub<f64> for DataSet<f64> {
    type Output = CDHResult<DataSet<f64>>;

    fn sub(self, rhs: f64) -> Self::Output {
        let data: Vec<f64> = self.data.iter()
            .map(|&x| x - rhs)
            .collect();

        Ok(DataSet{data})
    }
}

// =========================================================================
// 3. IN-PLACE SUBTRACTION (Performance Optimization)
// =========================================================================
impl SubAssign<f64> for DataSet<f64> {
    fn sub_assign(&mut self, rhs: f64) {
        // Modifies the vector directly in memory without allocating a new Vec!
        self.data.iter_mut().for_each(|x| *x -= rhs);
    }
}



// =========================================================================
// DATASET * DATASET (Element-wise Hadamard Product)
// =========================================================================
impl Mul<DataSet<f64>> for DataSet<f64> {
    type Output = CDHResult<DataSet<f64>>;

    fn mul(self, other: DataSet<f64>) -> Self::Output {
        if self.data.len() != other.data.len() {
          return Err("Both Datasets should be of equal length".to_string());
        }

        let data: Vec<f64> = self.data.iter()
            .zip(other.data.iter())
            .map(|(a, b)| a * b)
            .collect();

        Ok(DataSet { data })
    }
}

// =========================================================================
// DATASET * f64 (Scalar Scaling)
// =========================================================================
impl Mul<f64> for DataSet<f64> {
    type Output = DataSet<f64>;

    fn mul(self, rhs: f64) -> Self::Output {
        let data: Vec<f64> = self.data.iter()
            .map(|&x| x * rhs)
            .collect();

        DataSet { data }
    }
}

// =========================================================================
// DATASET *= f64 (In-place Scaling - Zero Heap Allocations)
// =========================================================================
impl MulAssign<f64> for DataSet<f64> {
    fn mul_assign(&mut self, rhs: f64) {
        self.data.iter_mut().for_each(|x| *x *= rhs);
    }
}


// =========================================================================
// DATASET / DATASET (Vector Division)
// =========================================================================
impl Div<DataSet<f64>> for DataSet<f64> {
    type Output = DataSet<f64>;

    fn div(self, other: DataSet<f64>) -> Self::Output {
        assert_eq!(
            self.data.len(),
            other.data.len(),
            "Both Datasets should be of equal length"
        );

        let data: Vec<f64> = self.data.iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| a / b)
            .collect();

        DataSet { data }
    }
}

// =========================================================================
// DATASET / f64 (Scalar Division)
// =========================================================================
impl Div<f64> for DataSet<f64> {
    type Output = DataSet<f64>;

    fn div(self, rhs: f64) -> Self::Output {
        let factor = 1.0 / rhs; // Reciprocal multiplication is optimized for CPU execution
        let data: Vec<f64> = self.data.iter()
            .map(|&x| x * factor)
            .collect();

        DataSet { data }
    }
}


// Vector Assignments (Requires dimension verification)
impl AddAssign<DataSet<f64>> for DataSet<f64> {
    fn add_assign(&mut self, rhs: DataSet<f64>) {
        assert_eq!(
            self.data.len(),
            rhs.data.len(),
            "Both Datasets should be of equal length"
        );
        self.data.iter_mut().zip(rhs.data.iter()).for_each(|(a, &b)| *a += b);
    }
}

impl SubAssign<DataSet<f64>> for DataSet<f64> {
    fn sub_assign(&mut self, rhs: DataSet<f64>) {
        assert_eq!(
            self.data.len(),
            rhs.data.len(),
            "Both Datasets should be of equal length"
        );
        self.data.iter_mut().zip(rhs.data.iter()).for_each(|(a, &b)| *a -= b);
    }
}

impl MulAssign<DataSet<f64>> for DataSet<f64> {
    fn mul_assign(&mut self, rhs: DataSet<f64>) {
        assert_eq!(
            self.data.len(),
            rhs.data.len(),
            "Both Datasets should be of equal length"
        );
        self.data.iter_mut().zip(rhs.data.iter()).for_each(|(a, &b)| *a *= b);
    }
}

impl DivAssign<DataSet<f64>> for DataSet<f64> {
    fn div_assign(&mut self, rhs: DataSet<f64>) {
        assert_eq!(
            self.data.len(),
            rhs.data.len(),
            "Both Datasets should be of equal length"
        );
        self.data.iter_mut().zip(rhs.data.iter()).for_each(|(a, &b)| *a /= b);
    }
}

impl PartialEq for DataSet<f64> {
    fn eq(&self, other: &Self) -> bool {
        if self.data.len() != other.data.len() {
            return false;
        }
        
        // Exact binary matching across vectors
        self.data.iter().zip(other.data.iter()).all(|(a, b)| a == b)
    }
}


