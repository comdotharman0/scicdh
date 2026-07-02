pub mod statistics;
pub mod probability;
pub mod series;
pub mod models;
use crate::models::regressions::*;
use crate::probability::*;
use crate::series::*;
 use crate::statistics::*;
fn main()->CDHResult<()>{
check_regressions()?;
let _ = check_probability()?;
let _ = check_series()?;
let d = DataSet:: new(vec![0.1,0.2,0.3,0.4,0.4,0.3,
0.2,0.1,0.0,0.4]);
let _ = d.info();
println!("Hello World {:#?}",Probability::factorial(100));
Ok(())
}






#[cfg(test)]
mod tests{
use crate::statistics::*;
use crate:: probability::*;
#[test]
fn mean_test(){
let a = DataSet(vec![1.0,2.0,3.0,4.0]);
let p= Probability::get_probability(8.0,3.0);
let npr= Probability::n_p_r(4,2);
let f5= Probability::factorial(5);
let ncr = Probability::n_c_r(4,2);
let random_var = Probability::random_variable(
&[1.0,2.0,3.0],|a|{
a*2 as f64
});
println!("Random Variable = {:?}",random_var);
assert_eq!(a.mean(),Ok(2.5f64));
assert_eq!(p,Ok(0.375f64));
assert_eq!(npr,Ok(12f64));
assert_eq!(ncr,Ok(6f64));
assert_eq!(f5,Ok(120f64));
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
