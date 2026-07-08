use scicdh::statistics::*;
use scicdh::probability::*;
use scicdh:: set;
use scicdh::traits::*;
use scicdh:: integrals::*;
fn main()->CDHResult<()>{
integrals_checking()?;
Ok(())
}






#[cfg(test)]
mod tests{
use scicdh::statistics::*;
use scicdh:: probability::*;
#[test]
fn mean_test(){
let a = DataSet{data:vec![1.0,2.0,3.0,4.0]};
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
