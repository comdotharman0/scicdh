///Result Type for the library
pub type CDHResult<T> = Result<T,String>;
///Probability Struct
#[derive(Debug)]
pub struct Probability;

#[derive(Debug)]
pub struct RandomVariable(pub Vec<f64>,
pub Vec<f64>);

#[derive(Debug)]
pub struct RandomVector<T>(pub Vec<RandomVariable>,
pub JointProbability<T>);

#[derive(Debug)]
pub struct ProbabilityVector<T>(pub Vec<Vec<T>>);

#[derive(Debug)]
pub struct JointProbability<T>(pub Vec<Vec<T>>);


impl RandomVariable{
pub fn new(x:Vec<f64>,px:Vec<f64>)->Self{
Self(x,px)
}
pub fn mean(&self)->CDHResult<f64>{
self.expectation_from_func(|a|{
a
})
}

pub fn variance(&self)->CDHResult<f64>{
let summat:f64 = self.0.iter().zip(self.1.iter()).
fold(0.0,|mut acc,b|{
acc+= b.0.powf(2.0)*b.1;

acc
});
Ok(
summat-self.mean()?.powf(2.0)
)
}

pub fn expectation_from_func<H>(&self,h:H)
->CDHResult<f64>
where H: Fn(f64)->f64
{Ok(
self.0.iter().zip(self.1.iter()).
fold(0.0,|mut a,b|{
a+=h(*b.0)*b.1;
a
}))
}
pub fn expectation_from_set(&self,hx: &[f64])
->CDHResult<f64>{
Ok(self.1.iter().zip(hx.iter()).
fold(0.0,|mut a, b|{
a+= b.0*b.1;
a
}))
}
pub fn std_deviation(&self)->CDHResult<f64>{
Ok(self.variance()?.powf(0.5))
}
pub fn moment_generating_func(&self,t:f64)
->CDHResult<f64>{
self.expectation_from_func(|a|{
let e: f64 = std::f64::consts::E;
e.powf(t*a)
})
}
pub fn characteristic_func(&self,t:f64)
->CDHResult<f64>{
self.expectation_from_func(|a|{
let e: f64 = std::f64::consts::E;
e.powf(t*a*(-1.0_f64).powf(0.5))
})
}

pub fn variance_operator<H>(&self,h:H)
->CDHResult<f64>
where
H:Clone+ Fn(f64)->f64
{
self.expectation_from_func(|a|{
(h(a) - self.expectation_from_func(h.clone())
.unwrap())
.powf(2.0)
})

}

}



impl RandomVector<f64>{
pub fn new(x:Vec<RandomVariable>,
px:JointProbability<f64>)->Self{
Self(x,px)
}
pub fn marginal_x2(&self,index:usize)
->CDHResult<f64>{
Ok(self.1.0.iter().map(|a| a[index]).sum::<f64>())
}
pub fn marginal_x1(&self,index:usize)
->CDHResult<f64>{
Ok(self.1.0[index].iter().sum::<f64>())
}
pub fn conditional_x1(&self,
index_x1:usize,
index_x2:usize)->CDHResult<f64>{
Ok(
self.1.0[index_x2][index_x1]/self.0[1].1[index_x2]
)
}
pub fn conditional_x2(&self,
index_x1:usize,
index_x2:usize)->CDHResult<f64>{
Ok(
self.1.0[index_x2][index_x1]/self.0[1].1[index_x1]
)}
pub fn conditional_expectation_x1<H>(&self,
 index_x2: usize, h: H) -> CDHResult<f64>
    where
        H: Fn(f64) -> f64,
    {
        let x1_outcomes = &self.0[0].0; 
// Access outcomes of X1
        let mut expected_value = 0.0;

        // Sum over all possible values of X1
        for index_x1 in 0..x1_outcomes.len() {
            let p_cond = self.conditional_x1(
index_x1, index_x2)?;
            if p_cond > 0.0 {
                expected_value += h(
x1_outcomes[index_x1]) * p_cond;
            }
        }

        Ok(expected_value)
    }

    /// Equation (4-19): E[h(Y) | X1 = x1] = \sum_{x2} h(x2) * p_{X2 | X1}(x2)
    /// Iterates through all possible outcomes of X2 for a fixed X1 index.
    pub fn conditional_expectation_x2<H>(&self, 
index_x1: usize, h: H) -> CDHResult<f64>
    where
        H: Fn(f64) -> f64,
    {
        let x2_outcomes = &self.0[1].0; 
// Access outcomes of X2
        let mut expected_value = 0.0;

        // Sum over all possible values of X2
        for index_x2 in 0..x2_outcomes.len() {
            let p_cond = self.conditional_x2(
index_x1, index_x2)?;
            if p_cond > 0.0 {
                expected_value += h(
x2_outcomes[index_x2]) * p_cond;
            }
        }

        Ok(expected_value)
    }
pub fn joint_expectation<F>(&self, g: F)
 -> CDHResult<f64>
    where
        F: Fn(f64, f64) -> f64,
    {
        let x = &self.0[0];
        let y = &self.0[1];
        let joint_matrix = &self.1.0;

        let mut expected_value = 0.0;

        for i in 0..x.0.len() {
            for j in 0..y.0.len() {
                let p_joint = joint_matrix[i][j];
                
                if p_joint > 0.0 {
                    // Accumulate: g(x, y) * P(x, y)
                    expected_value += g(
x.0[i], y.0[j]) * p_joint;
                }
            }
        }

       Ok( expected_value)
    }
/// 2. THE COVARIANCE CALCULATOR
    /// Uses your exact definitional 
///formula: E[(X1 - E(X1))(X2 - E(X2))]
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
pub fn get_probability(sample_space_len: f64,
favourable_len: f64)->CDHResult<f64>{
Ok(
favourable_len/sample_space_len)
}

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

pub fn n_c_r(n:usize,r:usize)->CDHResult<f64>{
Ok(Self::n_p_r(n,r)?/Self::factorial(r)?)
}
pub fn random_variable<F>(set: &[f64],
mut f:F)-> CDHResult<Vec<f64>>
where 
F: FnMut(f64)->f64,
//T: Copy
{Ok(
set.into_iter().map(|elements|{
f(*elements)
}).collect::<Vec<f64>>())
}

pub fn distribution_function(
probability_set:&[f64])->CDHResult<Vec<f64>>{
let mut distri: Vec<f64>= Vec::with_capacity(
probability_set.len());
let mut sum = 0f64;
for &prob in probability_set{
sum+= prob;
distri.push(sum);
}
Ok(
distri)
}

}


pub struct Binomial{
pub n:usize,pub x:Vec<usize>,pub p:f64}
impl Binomial{
pub fn new(n:usize,x:Vec<usize>,p:f64)->Self{
Self{
n,x,p
}
}
pub fn get_probability_set(&self)
->CDHResult<Vec<f64>>{
Ok(self.x.iter().map(|&i|{
let ncr = Probability::n_c_r(self.n,i).unwrap();
let ppx = self.p.powf(i as f64);
let qpnx = (1f64-self.p).powf((self.n-i) as f64);
ncr*ppx*qpnx
}).collect::<Vec<f64>>()
)}
pub fn mean(&self)->CDHResult<f64>{
Ok(self.n as f64 * self.p)
}
pub fn variance(&self)->CDHResult<f64>{
Ok(self.mean()?*(1f64-self.p))
}
}

pub struct Geometric{
pub x:Vec<usize>,pub p:f64}
impl Geometric{
pub fn new(x:Vec<usize>, p:f64)->Self{
Self{
x,p}
}
pub fn get_probability_set(&self)
->CDHResult<Vec<f64>>{
Ok(self.x.iter().map(|&i|{
let qpnx = (1f64-self.p).powf(i as f64 -1f64);
self.p*qpnx
}).collect::<Vec<f64>>()
)}
pub fn mean(&self)->CDHResult<f64>{
Ok(1f64/ self.p)
}
pub fn variance(&self)->CDHResult<f64>{
Ok((1f64-self.p)/(self.p.powf(2f64)))
}
}

pub struct Pascal{
pub r:usize,pub x:Vec<usize>,pub p:f64}
impl Pascal{
pub fn new(r:usize,x:Vec<usize>,p:f64)->Self{
Self{r,x,p}
}
pub fn get_probability_set(&self)
->CDHResult<Vec<f64>>{
Ok(self.x.iter().map(|&i|{
let ncr = Probability::n_c_r(i-1,self.r-1).unwrap();
let ppx = self.p.powf(i as f64);
let qpnx = (1f64-self.p).powf((i-self.r) as f64);
ncr*ppx*qpnx
}).collect::<Vec<f64>>()
)}
pub fn mean(&self)->CDHResult<f64>{
Ok(self.r as f64/self.p)
}
pub fn variance(&self)->CDHResult<f64>{
Ok(self.r as f64/(self.p.powf(2f64)))
}
}

pub struct HyperGeometric{
pub big_n:usize,pub big_t: usize,
pub n:usize,pub x:usize}
impl HyperGeometric{
pub fn bew(big_n:usize,big_t: usize,
n:usize,x:usize)->Self{
Self{big_n,big_t,n,x}
}

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
