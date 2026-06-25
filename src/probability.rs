pub struct Probability;
pub enum RandomVariable{
Discrete,
Continuous,
Mixed
}
impl Probability{
pub fn get_probability(sample_space_len: usize,
favourable_len: usize)->f64{

favourable_len as f64/sample_space_len as f64
}

pub fn factorial(n:usize)->usize{
(1..=n).product()
}
pub fn n_p_r(n:usize,r:usize)->usize{
(0..r).fold(1,|a,b|{
a*(n-b)
})
}

pub fn n_c_r(n:usize,r:usize)->usize{
Self::n_p_r(n,r)/Self::factorial(r)
}
pub fn random_variable<F>(set: &[f64],
mut f:F)-> Vec<f64>
where 
F: FnMut(f64)->f64,
//T: Copy
{
set.into_iter().map(|elements|{
f(*elements)
}).collect::<Vec<f64>>()
}

pub fn distribution_function(
probability_set:&[f64])->Vec<f64>{
let mut distri: Vec<f64>= Vec::with_capacity(
probability_set.len());
let mut sum = 0f64;
for &prob in probability_set{
sum+= prob;
distri.push(sum);
}
distri
}
pub fn discrete_random_variable_mean(){}
pub fn continuous_random_variable_mean(){}
pub fn discrete_random_variable_variance(){}
pub fn continuous_random_variable_variance(){}
pub fn discrete_random_variable_expectation(){}
pub fn continuous_random_variable_expectation(){}
pub fn discrete_random_variable_variance_op(){}


}

impl RandomVariable{
pub fn mean(&self,x: &[f64],px:&[f64])->f64{
match self{
 &RandomVariable::Discrete=>{

},
 &RandomVariable::Continuous=>{},
&RandomVariable::Mixed=>{}
};
0f64
}
pub fn variance(&self,x: &[f64],px:&[f64])->f64{
match self{
 &RandomVariable::Discrete=>{

},
 &RandomVariable::Continuous=>{},
&RandomVariable::Mixed=>{}
};
0f64
}
pub fn expectation_op()->f64{0f64}
pub fn variance_op()->f64{0f64}





}
