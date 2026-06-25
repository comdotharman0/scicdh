pub mod statistics;
pub mod probability;
use crate::probability::Probability;
pub struct RandomVector<T>(Vec<Vec<T>>);
pub struct ProbabilityVector<T>(Vec<Vec<T>>);
pub struct JointProbability<T>(Vec<Vec<T>>);

pub struct RandomVariable(Vec<f64>,Vec<f64>);
pub struct Matrix<T>(Vec<Vec<T>>);
impl Matrix<isize>{
pub fn get(&self,i:usize,j:usize)->isize{
self.0[i][j]
}
pub fn set(&mut self,i:usize,j:usize,value:isize)
->isize{
self.0[i][j]= value;
value
}
pub fn add_row(&mut self,row: Vec<isize>){
self.0.push(row);
}

}
impl RandomVariable{
pub fn mean(&self)->f64{

self.expectation_from_func(|a|{
a
})
}

pub fn variance(&self)->f64{
let summat:f64 = self.0.iter().zip(self.1.iter()).
fold(0.0,|mut acc,b|{
acc+= b.0.powf(2.0)*b.1;

acc
});

summat-self.mean().powf(2.0)

}

pub fn expectation_from_func<H>(&self,h:H)->f64
where H: Fn(f64)->f64
{
self.0.iter().zip(self.1.iter()).
fold(0.0,|mut a,b|{
a+=h(*b.0)*b.1;
a
})
}
pub fn expectation_from_set(&self,hx: &Vec<f64>)->f64{
self.1.iter().zip(hx.iter()).
fold(0.0,|mut a, b|{
a+= b.0*b.1;
a
})
}
pub fn std_deviation(&self)->f64{
self.variance().powf(0.5)
}
pub fn moment_generating_func(&self,t:f64)->f64{
self.expectation_from_func(|a|{
let e: f64 = 2.718;
e.powf(t*a)
})
}
pub fn characteristic_func(&self,t:f64)->f64{
self.expectation_from_func(|a|{
let e: f64 = 2.718;
e.powf(t*a*(-1.0 as f64).powf(0.5))
})
}

pub fn variance_operator<H>(&self,h:H)->f64
where
H:Clone+ Fn(f64)->f64
{
self.expectation_from_func(|a|{
(h(a) - self.expectation_from_func(h.clone()))
.powf(2.0)
})

}
/*
pub fn covariance(&self,other: RandomVariable)->f64{
let ex1=self.expectation_from_func(
|a|{

}

)

}*/
}
fn main(){
println!("Hello World");
let x = vec![0.0,1.0,2.0,3.0];
let px = vec![1,3,3,1].iter().
map(|&a|{

Probability::get_probability(8,a as usize)
}).collect::<Vec<f64>>();
println!("x= {:?} \n px= {:?}",x,px);
let rv = RandomVariable(x,px);
println!("Random Variable = {:#?}",
rv.moment_generating_func(1.0));
let npr= Probability::n_p_r(4,2);
println!("n_p_r={npr}");
let random_var = Probability::random_variable(
&[1.0,2.0,3.0],|a|{
a*2 as f64
});
println!("Random Variable = {:?}",random_var);
let pd = Probability::distribution_function(
&[0.1,0.2,0.3,0.4]);
println!("Probability Distribution = {:#?}",pd);
let mut a= [1.0,2.0,3.0,4.0];
let mut b = [0.1,0.2,0.3,0.4];
println!("a zip b = {:#?}",a.iter().zip(b.iter())
.scan(
1.0,|c,d|{
*c = d.0*d.1;
Some(*c)
}).sum::<f64>());
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
