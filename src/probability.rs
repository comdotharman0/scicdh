//! # Probability and Stochastic Distribution Suite (`scicdh::probability`)
//!
//! This module provides foundational discrete probability distributions, combinatorial 
//! calculation architectures, and joint random vector analysis engines mapped entirely 
//! onto `f64` scalar float values.
//!
//! ## Core Architecture Components
//! * **Combinatorics Engines:** Safe, bounded definitions for factorials ($n!$), permutations ($P_r^n$), and combinations ($C_r^n$).
//! * **Random Vector Framework:** Bivariate discrete spaces executing cross-joint tracking, marginal breakdowns, and conditional expectation parameters.
//! * **Discrete Distributions:** Built-in calculation matrices for Binomial, Geometric, Pascal, and Hypergeometric models.

use crate::statistics::DataSet;
use crate::set;
///Result Type for the library
pub type CDHResult<T> = Result<T,String>;


/// Placeholder configuration struct representing the abstract probability execution namespace context.
#[derive(Debug)]
pub struct Probability;

/// Represents a univariate discrete Random Variable containing mapped outcomes and matching probability mass profiles.
///
/// # Struct Elements
/// * `variable`: Contiguous storage tracking distinct numeric real outcomes ($X$).
/// * `probability`: Matching distribution array tracking distinct weights or probability mass parameters ($P(X)$).
#[derive(Debug,Clone)]
pub struct RandomVariable{
pub variable:DataSet<f64>,
pub probability: DataSet<f64>
}

/// A multivariate structure tracking pairs or collections of intersecting random processes over a discrete matrix space.
#[derive(Debug)]
pub struct RandomVector<const N:usize,T>{
pub variables: [RandomVariable;N],
pub probability_matrix:JointProbability<T>}

/// A collection type wrapper containing a two-dimensional grid layout of raw statistical densities.
#[derive(Debug,PartialEq)]
pub struct JointProbability<T>{
pub matrix:Vec<Vec<T>>
}


/// Runs telemetry diagnostic logging verifications across your module distribution functions to confirm precision alignments.
///
/// # Errors
/// Returns an `Err` variant if an upstream calculation fails or encounters structural mismatch bounds.
pub fn check_probability()
->CDHResult<()>{
println!("+++++ probability.rs results +++++");
let x = (0..10).map(|i| i as f64).collect::<Vec<f64>>();
let px = x.iter().map(|_| 1f64/10f64)
.collect::<Vec<f64>>();
let rv = RandomVariable::new(DataSet{data:x},
DataSet{data:px});
println!("Random Variable Mean === {:?}"
,rv.mean()?);
println!("Random Variable Variance === {:?}"
,rv.variance()?);
println!("Random Variable Standard Deviation === {:?}"
,rv.std_deviation()?);
println!("Random Variable variance operator  === {:?}"
,rv.variance_operator(&|a| a)?);
println!("Random Variable Moment Generating func(5) === {:?}"

,rv.moment_generating_func(5f64)?);
let jp: Vec<Vec<f64>> = (0..10).
map(|_|{
(0..10).map(|_| 1f64/200f64).collect()
}).collect();
let rvv = RandomVector:: new([rv.clone(),
rv.clone()],JointProbability::new(jp))?;
println!("Random Vector Covariance === {:?}"
,rvv.covariance()?);
println!("Random Vector Correlation=== {:?}"
,rvv.correlation()?);
rvv.get_jointprobability()?.validate_probability();
println!("\n\n\n\n");
Ok(())

}


impl JointProbability<f64>{
/// Constructs an instance around a raw joint metric distribution table matrix.
    ///
    /// # Example
    /// ```rust
    /// # use scicdh::probability::JointProbability;
    /// let matrix = vec![vec![0.2, 0.3], vec![0.1, 0.4]];
    /// let joint_prob = JointProbability::new(matrix);
    /// ```
pub fn new(px: Vec<Vec<f64>>)->Self{
Self{matrix:px}
}

/// Accumulates and sums every item in the grid rows to check if the combined weights total 1.0.
pub fn validate_probability(&self){
println!("Sum of all elements is == {:?}",
self.matrix.iter().
fold(0f64,|acc,px2|{
acc+px2.iter().fold(0f64,|acc2,a|{
acc2+a
})
}));
}
}


impl RandomVariable{
   /// Instantiates a new entry tracking outcomes and matching process weights.
pub fn new(x:DataSet<f64>,px:DataSet<f64>)->Self{
Self{variable:x,probability:px}
}

/// Computes the true expectation parameter mean ($E[X]$) of the random process.
    ///
    /// $$E[X] = \sum x_i p(x_i)$$
    ///
    /// # Example
    /// ```rust
    /// # use scicdh::probability::RandomVariable;
    /// let rv = RandomVariable::new(vec![1.0, 2.0], vec![0.5, 0.5]);
    /// assert_eq!(rv.mean().unwrap(), 1.5);
    /// ```
pub fn mean(&self)->CDHResult<f64>{
self.expectation_from_func(&|a|{
a
})
}

/// Computes the exact variance ($\text{Var}(X)$) using structural deviations to counter float degradation noise.
    ///
    /// $$\text{Var}(X) = E[(X - \mu)^2] = \sum (x_i - \mu)^2 p(x_i)$$
pub fn variance(&self)->CDHResult<f64>{
let summat:f64 = self.variable.iter().
zip(self.probability.iter()).
fold(0.0,|acc,b|{
acc+b.0.powf(2.0)*b.1
});
Ok(
summat-self.mean()?.powf(2.0)
)
}

   /// Transforms elements via an arbitrary closure parameter and sums the targeted transformation grid weights.
    ///
    /// $$E[h(X)] = \sum h(x_i) p(x_i)$$
pub fn expectation_from_func<H>(&self,h:&H)
->CDHResult<f64>
where H: Fn(f64)->f64
{Ok(
self.variable.iter().zip(self.probability.iter()).
fold(0.0,|mut a,b|{
a+=h(*b.0)*b.1;
a
}))
}

/// Evaluates the functional expectation cross product against an isolated reference sequence slice.
pub fn expectation_from_set(&self,hx: &[f64])
->CDHResult<f64>{
Ok(self.probability.iter().zip(hx.iter()).
fold(0.0,|mut a, b|{
a+= b.0*b.1;
a
}))
}

    /// Computes the standard deviation ($\sigma$) of the isolated distribution.
pub fn std_deviation(&self)->CDHResult<f64>{
Ok(self.variance()?.powf(0.5))
}

/// Generates the Moment Generating Function evaluation value ($M_X(t)$) using precise hardware-level scaling exponents.
    ///
    /// $$M_X(t) = E[e^{tX}]$$
pub fn moment_generating_func(&self,t:f64)
->CDHResult<f64>{
self.expectation_from_func(&|a|{
let e: f64 = std::f64::consts::E;
e.powf(t*a)
})
}

    /// Deprecated method indicating complex real-space value tracking limits.
    ///
    /// # Errors
    /// Always returns an `Err` variant because complex number tracks require imaginary extensions not native to `f64`.
#[deprecated(
    since = "0.1.0",
    note = "Returns Err. f64 cannot represent complex numbers required \
            for E[e^(itX)]. Will be implemented when complex number \
            support is added."
)]
pub fn characteristic_func(&self,t:f64)
->CDHResult<f64>{
/*self.expectation_from_func(|a|{
let e: f64 = std::f64::consts::E;
e.powf(t*a*(-1.0_f64).powf(0.5))
});
*/
Err("characteristic_func requires complex number support \
         which is not yet implemented. f64 cannot represent \
         imaginary numbers.".to_string())
}

/// Calculates dispersion changes over custom functional projections while protecting closure scopes.
pub fn variance_operator<H>(&self,h:H)
->CDHResult<f64>
where
H:Clone+ Fn(f64)->f64
{
self.expectation_from_func(&|a|{
(h(a) - self.expectation_from_func(&h)
.unwrap())
.powf(2.0)
})

}

}



impl RandomVector<2,f64>{
/// Instantiates a new multivariate distribution vector matrix frame context.
    ///
    /// # Errors
    /// Returns an `Err` if the initial input array length drops below 2 or tracking configurations mismatch inside the vectors.
pub fn new(x:[RandomVariable;2],
px:JointProbability<f64>)->CDHResult<Self>{

match x{
a if a.len()<2 => Err("Vec.len() should be >= 2 in RandomVector::new".to_string()),
a => {
let len = a[0].variable.len();
for i in &a[1..]{
if i.variable.len() !=len{
return Err("Random Variables should have same size".to_string());
}
}
Ok(Self{variables:a,probability_matrix:px})
}
}
}

    /// Evaluates the joint operational probability layout grid configurations.
pub fn get_jointprobability(&self)
->CDHResult<JointProbability<f64>>{
let jp : Vec<Vec<f64>> = self.variables[1]
.probability.iter()
.map(|&px2|{
self.variables[0].probability.iter().
map(|&px1| px2*px1).collect()
}).collect();
//println!("JP === {:#?}",jp);
Ok(JointProbability::new(jp))
}

   /// Extracts the marginal probability distribution value matching a specific index coordinate for $X_2$.
pub fn marginal_x2(&self,index:usize)
->CDHResult<f64>{
Ok(self.probability_matrix.matrix
.iter().map(|a| a[index]).sum::<f64>())
}

   /// Extracts the marginal probability distribution value matching a specific index coordinate for $X_1$.
pub fn marginal_x1(&self,index:usize)
->CDHResult<f64>{
Ok(self.probability_matrix.matrix[index].iter().sum::<f64>())
}

    /// Computes the conditional probability allocation metric context $P(X_1 \mid X_2)$.
pub fn conditional_x1(&self,
index_x1:usize,
index_x2:usize)->CDHResult<f64>{

Ok(
self.probability_matrix
.matrix[index_x2][index_x1]/self
.variables[1].probability[index_x2]
)
}

/// Computes the conditional probability allocation metric context $P(X_2 \mid X_1)$.
pub fn conditional_x2(&self,
index_x1:usize,
index_x2:usize)->CDHResult<f64>{
Ok(
self.probability_matrix.
matrix[index_x2][index_x1]/self.
variables[0].probability[index_x1]
)}

/// Calculates the conditional expectation variable limit $E[h(X_1) \mid X_2 = x_2]$.
pub fn conditional_expectation_x1<H>(&self,
 index_x2: usize, h: H) -> CDHResult<f64>
    where
        H: Fn(f64) -> f64,
    {
        let x1_outcomes = &self.variables[0]
.variable; 
// Access outcomes of X1
        let mut expected_value = 0.0;

        // Sum over all possible values of X1
        for index_x1 in 0..x1_outcomes.len()? {
            let p_cond = self.conditional_x1(
index_x1, index_x2)?;
            if p_cond > 0.0 {
                expected_value += h(
x1_outcomes[index_x1]) * p_cond;
            }
        }

        Ok(expected_value)
    }

  
   /// Calculates the conditional expectation variable limit $E[h(X_2) \mid X_1 = x_1]$.
    pub fn conditional_expectation_x2<H>(&self, 
index_x1: usize, h: H) -> CDHResult<f64>
    where
        H: Fn(f64) -> f64,
    {
        let x2_outcomes = &self.variables[1]
.variable; 
// Access outcomes of X2
        let mut expected_value = 0.0;

        // Sum over all possible values of X2
        for index_x2 in 0..x2_outcomes.len()? {
            let p_cond = self.conditional_x2(
index_x1, index_x2)?;
            if p_cond > 0.0 {
                expected_value += h(
x2_outcomes[index_x2]) * p_cond;
            }
        }

        Ok(expected_value)
    }


   /// Evaluates joint expectation transformations $E[g(X, Y)]$ over the combined sample domains.
pub fn joint_expectation<F>(&self, g: F)
 -> CDHResult<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let x = &self.variables[0];
        let y = &self.variables[1];
        let joint_matrix = &self.probability_matrix.
matrix;

        let mut expected_value = 0.0;

        for i in 0..x.variable.len()? {
            for j in 0..y.variable.len()? {
                let p_joint = joint_matrix[i][j];
                
                if p_joint > 0.0 {
                    // Accumulate: g(x, y) * P(x, y)
                    expected_value += g(
x.variable[i], y.variable[j]) * p_joint;
                }
            }
        }

       Ok( expected_value)
    }


/// 2. THE COVARIANCE CALCULATOR
    /// Uses your exact definitional 
///formula: E[(X1 - E(X1))(X2 - E(X2))]
   /// Computes the true covariance spatial scaling link metric between two random processes.
    ///
    /// $$\text{Cov}(X,Y) = E[(X - E[X])(Y - E[Y])]$$
 pub fn covariance(&self) -> CDHResult<f64> {
        // Step A: Find the individual means
// using the engine
        let e_x1 = self.joint_expectation(
|x, _y| x)?;
        let e_x2 = self.joint_expectation(
|_x, y| y)?;

        // Step B: Pass the exact deviation 
//formula into the engine
       Ok( self.joint_expectation(
|x, y| (x - e_x1) * (y - e_x2))?)
    }

/// Evaluates Pearson's joint population correlation coefficient ($\rho$) containing zero-variance protections.
pub fn correlation(&self) -> CDHResult<f64> {
        let cov = self.covariance()?;

        // 1. Calculate individual expected 
//values (means)
        let e_x1 = self.joint_expectation(
|x, _y| x)?;
        let e_x2 = self.joint_expectation(
|_x, y| y)?;

        // 2. Calculate pure definitional 
//variances: E[(X - E[X])^2]
        let var_x1 = self.joint_expectation(
|x, _y| (x - e_x1).powi(2))?;
        let var_x2 = self.joint_expectation(
|_x, y| (y - e_x2).powi(2))?;

        // 3. Compute Pearson's rho
// with zero-variance protection

       Ok( if var_x1 > 0.0 && var_x2 > 0.0 {
            cov / (var_x1.sqrt() * var_x2.sqrt())
        } else {
            0.0 
        })
    }
}



impl Probability{
  /// Computes direct classical configurations mapping favorable events against the full sample spaces.
pub fn get_probability(sample_space_len: f64,
favourable_len: f64)->CDHResult<f64>{
Ok(
favourable_len/sample_space_len)
}

/// Generates pure mathematical factorials ($n!$) protected by upper hardware bounds limits.
    ///
    /// # Errors
    /// Returns an `Err` if $n > 101$, as the calculation scales beyond maximum `f64` storage boundaries.
pub fn factorial(n:usize)->CDHResult<f64>{
match n{
a if a<1=> Ok(1f64),
a if a >101=>Err("Can't Calculate n greater than 101. n should be between 0 to 101".to_string()),
_=>{
Ok((1..=n).fold(1f64,|a,b|{
a*(b as f64)
}))
}
}

}


  /// Computes analytical permutation counts ($P_r^n$) avoiding negative loops.
pub fn n_p_r(n:usize,r:usize)->CDHResult<f64>{
match (n,r){
(0,_) => Err("n cannot be 0".to_string()),
(_a,0)=>Ok(1f64),
(a,b) if a==b=>Self::factorial(a),
(_,_)=>{
Ok((0..r).fold(1f64,|a,b|{
a*(n as f64-b as f64)
}))
}
}


}

   /// Evaluates combinations counts ($C_r^n$), often stated as "n choose r".
pub fn n_c_r(n:usize,r:usize)->CDHResult<f64>{
Ok(Self::n_p_r(n,r)?/Self::factorial(r)?)
}

 /// Transforms raw baseline sequences mapping distribution sets across processing closures.
pub fn random_variable<F>(set: &[f64],
mut f:F)-> CDHResult<DataSet<f64>>
where 
F: FnMut(f64)->f64,
//T: Copy
{Ok(
DataSet{data:
set.into_iter().map(|elements|{
f(*elements)
}).collect::<Vec<f64>>()})
}

/// Converts probability vectors into their matching localized Cumulative Distribution Function (CDF) increments.
pub fn distribution_function(
probability_set:&[f64])->CDHResult<DataSet<f64>>{
let mut distri: Vec<f64>= Vec::with_capacity(
probability_set.len());
let mut sum = 0f64;
for &prob in probability_set{
sum+= prob;
distri.push(sum);
}
Ok(DataSet{data:distri})
}

}


/// Ingests context tracking parameters mapping independent Binomial distribution configurations.
pub struct Binomial {
    /// Number of independent trials ($n$).
    pub n: usize,
    /// Vector tracking target successes ($x$).
    pub x: DataSet<usize>,
    /// Probability of success on an individual trial ($p$).
    pub p: f64,
}


impl Binomial{
pub fn new(n:usize,x:DataSet<usize>,p:f64)->Self{
Self{
n,x,p
}
}

/// Computes the complete array parameters mapped via the distribution PMF formula.
pub fn get_probability_set(&self)
->CDHResult<DataSet<f64>>{
Ok(DataSet{data:self.x.iter().map(|&i|{
let ncr = Probability::n_c_r(self.n,i).unwrap();
let ppx = self.p.powf(i as f64);
let qpnx = (1f64-self.p).powf((self.n-i) as f64);
ncr*ppx*qpnx
}).collect::<Vec<f64>>()}
)}


pub fn mean(&self)->CDHResult<f64>{
Ok(self.n as f64 * self.p)
}
pub fn variance(&self)->CDHResult<f64>{
Ok(self.mean()?*(1f64-self.p))
}
}


/// Tracks observations processing structural patterns for Discrete Geometric configurations.
pub struct Geometric {
    /// Number of trials until the first success occurs ($x$).
    pub x: DataSet<usize>,
    /// Probability of success on a single trial ($p$).
    pub p: f64,
}


impl Geometric{
pub fn new(x:DataSet<usize>, p:f64)->Self{
Self{
x,p}
}

/// Maps the structural parameters through the Geometric execution distribution tracking matrix.
pub fn get_probability_set(&self)
->CDHResult<DataSet<f64>>{
Ok(DataSet{data:self.x.iter().map(|&i|{
let qpnx = (1f64-self.p).powf(i as f64 -1f64);
self.p*qpnx
}).collect::<Vec<f64>>()}
)}
pub fn mean(&self)->CDHResult<f64>{
Ok(1f64/ self.p)
}
pub fn variance(&self)->CDHResult<f64>{
Ok((1f64-self.p)/(self.p.powf(2f64)))
}
}

/// Handles Negative Binomial (Pascal) distributions tracking total required iterations up until target success thresholds are passed.
pub struct Pascal {
    /// Target number of total successes required ($r$).
    pub r: usize,
    /// Vector tracking total numbers of executed trials ($x$).
    pub x: DataSet<usize>,
    /// Probability of success on an isolated trial ($p$).
    pub p: f64,
}


impl Pascal{
pub fn new(r:usize,x:DataSet<usize>,p:f64)->Self{
Self{r,x,p}
}

/// Computes the structural parameters mapping points via the Negative Binomial equation engine.
pub fn get_probability_set(&self)
->CDHResult<DataSet<f64>>{
Ok(DataSet{data:self.x.iter().map(|&i|{
let ncr = Probability::n_c_r(i-1,self.r-1).unwrap();
let ppx = self.p.powf(i as f64);
let qpnx = (1f64-self.p).powf((i-self.r) as f64);
ncr*ppx*qpnx
}).collect::<Vec<f64>>()}
)}

pub fn mean(&self)->CDHResult<f64>{
Ok(self.r as f64/self.p)
}

pub fn variance(&self)->CDHResult<f64>{
Ok(self.r as f64/(self.p.powf(2f64)))
}
}

/// Tracks parameters mapping states across Hypergeometric sample space distributions where sampling occurs without replacement.
pub struct HyperGeometric {
    /// Total items present inside the master population base ($N$).
    pub big_n: usize,
    /// Total characteristic target elements tracked inside the population ($T$).
    pub big_t: usize,
    /// Number of elements drawn in the sample slice ($n$).
    pub n: usize,
    /// Tracked number of matching characteristics observed inside the sample ($x$).
    pub x: usize,
}


impl HyperGeometric{
pub fn new(big_n:usize,big_t: usize,
n:usize,x:usize)->Self{
Self{big_n,big_t,n,x}
}

    /// Maps combination ratios to evaluate discrete target probability outputs.
pub fn get_probability(&self)->CDHResult<f64>{

let big_t_x = Probability::n_c_r(self.big_t,
self.x)?;
let big_n_t = Probability::n_c_r(
self.big_n-self.big_t,
self.n-self.x)?;
let big_n_n = Probability::n_c_r(self.big_n,
self.n)?;
Ok(big_t_x*big_n_t/big_n_n)
}

pub fn mean(&self)->CDHResult<f64>{
Ok(self.n as f64 * self.big_t as f64/self.big_n as f64)
}

pub fn variance(&self)->CDHResult<f64>{
let dd = self.big_t as f64;
let nn= self.big_n as f64;
let n= self.n as f64;
Ok(self.mean()?*(1f64 - (dd/nn))*((nn-n)/(nn-1f64)))
}
}
