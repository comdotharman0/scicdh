pub mod statistics;
pub mod probability;
use crate::probability::*;

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
Ok(0f64)
}
pub fn mode(&self)->CDHResult<f64>{

Ok(0f64)
}
pub fn range(&self)->CDHResult<f64>{
/*self.0.iter().max().unwrap_or_else(|e|{
format!("The error is {:#?}",e)})-self.0.iter().min().unwrap()
*/ Ok(0f64)
}
pub fn variance(&self)->CDHResult<f64>{
let n = self.0.len() as f64;
let mean = self.mean()?;
Ok(self.0.iter().fold(0.0,|acc,&x|{
acc+ (x-mean).powf(2f64)

})/(n-1f64))
}
pub fn info(&self)->CDHResult<()>{
let (mean,median,mode,range,variance)= (self.mean()?,
self.median()?,self.mode()?,self.range()?,
self.variance()?);
println!("DataSet === {:#?}",self.0);
println!("Mean === {:#?}",mean);
println!("Median === {:#?}",median);
println!("Mode === {:#?}",mode);
println!("Range === {:#?}",range);
println!("Variance === {:#?}",variance);
Ok(())
}
}

fn main()->CDHResult<()>{
let d = DataSet:: new((0..11).map(|i| i as f64).collect::<Vec<f64>>());
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
