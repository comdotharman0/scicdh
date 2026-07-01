pub mod statistics;
pub mod probability;
use crate::probability::*;
use std::collections::HashMap;
pub struct DataSet<T>(Vec<T>);
impl DataSet<f64>{
pub fn new(data:Vec<f64>)->Self{
DataSet(data)

}
pub fn mean(&self)->CDHResult<f64>{
let len = self.0.len() as f64;
Ok(self.0.iter().sum::<f64>()/len)

}
pub fn median(&self)->CDHResult<f64>{
let len = self.0.len();
let mut  sorted_arr= self.0.clone();
sorted_arr.sort_by(|a,b| a.total_cmp(b));
match len{
a if a%2==1=>Ok(sorted_arr[(a+1)/2 -1]),
_=>{
let center = len/2;
let center_elem= sorted_arr[center-1];
let next_elem= sorted_arr[center];
Ok((center_elem+next_elem)/2f64)
}
}
}
pub fn mode(&self)->CDHResult<Option<f64>>{
let mut frequencies: Vec<(f64,usize)> = vec![];
let set = self.into_set()?;
let mut max_count = 0;
set.iter()
.for_each(|&a|{
let i = self.0.iter().filter(|&&b| b==a).count();
if i>max_count{
max_count = i;
}
frequencies.push((a,i));
});
let mut  desired = frequencies
.iter().filter(|&a| a.1 == max_count);
if  desired.clone().count()>1{
Ok(None)
}else{
Ok(Some(desired.next().ok_or("")?.0))
}
}
pub fn range(&self)->CDHResult<f64>{
/*self.0.iter().max().unwrap_or_else(|e|{
format!("The error is {:#?}",e)})-self.0.iter().min().unwrap()
*/ Ok(self.max()?-self.min()?)
}
pub fn has(&self,value:f64,other:Option<&[f64]>)
->CDHResult<bool>{
let set = match other{
None=> &self.0,
_ => &other.ok_or("")?.to_vec()
};
Ok(set.iter().filter(|&&v| v==value).count()>0)
}

pub fn sample_variation_x2(&self)->CDHResult<f64>{
let n = self.0.len() as f64;
let xsqrsum= self.0.iter()
.fold(0.0,|acc,b| acc+b.powf(2f64));
let xsum= self.0.iter()
.fold(0.0,|acc,b| acc+b);
Ok(
xsqrsum-xsum.powf(2f64)/n
)
}
pub fn sample_variation_y2(&self,y:&[f64])
->CDHResult<f64>{
let n = y.len() as f64;
let ysqrsum= y.iter()
.fold(0.0,|acc,b| acc+b.powf(2f64));
let ysum= y.iter()
.fold(0.0,|acc,b| acc+b);
Ok(
ysqrsum-ysum.powf(2f64)/n
)
}
pub fn sample_variance(&self)->CDHResult<f64>{
Ok(self
.sample_variation_x2()?/((self.0.len()-1) as f64))
}
pub fn sample_std_deviation(&self)->CDHResult<f64>{
Ok(self.sample_variance()?.powf(0.5))
}

pub fn sample_coeff_of_variation(&self)
->CDHResult<f64>{
Ok(self.sample_std_deviation()?/self.mean()?)
}

pub fn kth_central_moment(&self,k:f64)
->CDHResult<f64>{
let n = self.0.len() as f64;
let mean = self.mean()?;
Ok(self.0.iter()
.fold(0f64,|acc,b| acc+(b-mean).powf(k))/n)
}

pub fn skewness(&self)->CDHResult<f64>{
Ok(self.kth_central_moment(3f64)?/self.
population_std_deviation()?.powf(3f64))
}

pub fn kurtosis(&self)->CDHResult<f64>{
Ok((self.kth_central_moment(4f64)?/self.
population_std_deviation()?.powf(4f64)) -3f64)
}

pub fn sample_variation_xy(&self,y:&[f64])
->CDHResult<f64>{
 match y.len(){
a if a ==self.0.len()=>{
let xysum= self.0.iter().zip(y.iter())
.fold(0.0,|acc,b| acc+(b.0*b.1) );
let xsum= self.0.iter()
.fold(0.0,|acc,b| acc+b);
let ysum= y.iter()
.fold(0.0,|acc,b| acc+b);
let std_variation_xy= xysum-(xsum*ysum/y
.len() as f64);

Ok(std_variation_xy)
},
_=> Err("Length of y should be equal to length of dataset".to_string())
}
}

pub fn pearson_correlation_coeff(&self,y:&[f64])
->CDHResult<f64>{
Ok(self.sample_variation_xy(y)?/((self
.sample_variation_x2()?*self.sample_variation_y2(y)?)
.powf(0.5)))
}

pub fn into_set(&self)->CDHResult<Vec<f64>>{
let mut set: Vec<f64> = vec![];
for &a in self.0.iter(){
if !self.has(a,Some(&set))?{
set.push(a);
}
}
Ok(set)
}
pub fn population_variance(&self)->CDHResult<f64>{
let n = self.0.len() as f64;
let mean = self.mean()?;
Ok(self.0.iter().fold(0.0,|acc,&x|{
acc+ (x-mean).powf(2f64)

})/n)
}
pub fn population_std_deviation(&self)
->CDHResult<f64>{
Ok(self.population_variance()?.powf(0.5))
}
pub fn max(&self)->CDHResult<f64>{
let mut max = self.0[0];
self.0.iter().
for_each(|&i| {
if i>max{
max = i;
}
});
Ok(max)
}
pub fn min(&self)->CDHResult<f64>{
let mut min = self.0[0];
self.0.iter().
for_each(|&i| {
if i<min{
min = i;
}
});
Ok(min)
}
pub fn info(&self)->CDHResult<()>{
let (mean,median,mode,range,variance)= (self.mean()?,
self.median()?,self.mode()?,self.range()?,
self.sample_variance()?);
println!("DataSet === {:?}",self.0);
println!("Mean === {:#?}",mean);
println!("Median === {:#?}",median);
println!("Mode === {:#?}",mode);
println!("Range === {:#?}",range);
println!("Variance === {:#?}",variance);
println!("Skewness === {:#?}",self.skewness()?);
println!("Kurtosis === {:#?}",self.kurtosis()?);
println!("self has 3f64 == {:?}",self.has(3f64,None)?);
println!("Set === {:?}",self.into_set()?);
println!("Maximum === {:?}",self.max()?);
println!("Minimum === {:?}",self.min()?);
Ok(())
}
}

fn main()->CDHResult<()>{
let _ = check_probability()?;
let d = DataSet:: new(vec![0.1,0.2,0.3,0.4,0.4,0.3,
0.2,0.1,0.0,0.4]);
let _ = d.info();
println!("Hello World {:#?}",Probability::factorial(100));
Probability::factorial(103)?;
Ok(())
}






#[cfg(test)]
mod tests{
use crate::statistics::*;
use crate:: probability::*;
#[test]
fn mean_test(){
let a = DataSet(vec![1.0,2.0,3.0,4.0]);
let p= Probability::get_probability(8,3);
let npr= Probability::n_p_r(4,2);
let f5= Probability::factorial(5);
let ncr = Probability::n_c_r(4,2);
let random_var = Probability::random_variable(
&[1.0,2.0,3.0],|a|{
a*2 as f64
});
println!("Random Variable = {:?}",random_var);
assert_eq!(a.mean(),2.5);
assert_eq!(p,0.375);
assert_eq!(npr,12);
assert_eq!(ncr,6);
assert_eq!(f5,120);
}
}










/*
fn main() {
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // 1. Hide the cursor so it doesn't flicker
    write!(handle, "\x1b[?25l").unwrap();

    for i in 0..20 {
        // 2. Clear the screen
        write!(handle, "\x1b[2J").unwrap();
        
        // 3. Draw a "graph" based on the loop index
        write!(handle, "\x1b[H").unwrap(); // Move cursor to top-left
        println!("--- Minimal TUI Engine ---");
        
        for y in 0..10 {
            if y == i % 10 {
                println!("* <--- Data Point at {}", i);
            } else {
                println!("|");
            }
        }
        
        handle.flush().unwrap();
        thread::sleep(time::Duration::from_millis(200));
    }

    // 4. Show the cursor again before exiting
    write!(handle, "\x1b[?25h").unwrap();
}

*/
