pub mod statistics;
pub mod probability;
use crate::probability::*;

pub struct DataSet<T>(Vec<T>);
impl DataSet<f64>{
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
}

fn main(){
println!("Hello World {:#?}",Probability::factorial(100));

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
