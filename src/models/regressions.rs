use crate::statistics::DataSet;
use crate::probability::CDHResult;
pub fn check_regressions()->CDHResult<()>{
let x_set: Vec<f64> = (0..10)
.map(|i| i as f64).collect();
let mut lr = LinearRegression::new(
x_set.clone(),x_set
.iter().map(|&i| 2f64*i + 1f64).collect());
let pred = lr.predict(vec![11f64,12f64,13f64])?;
println!("Predictions == {:#?}",pred);
Ok(())
}

#[derive(Debug)]
pub struct LinearRegression<T>{
pub x_set: Vec<T>,
pub y_set: Vec<T>,
pub m: T,
pub b: T

}

impl LinearRegression<f64>{

pub fn new(x_set:Vec<f64>,y_set:Vec<f64>)-> Self{
let n = x_set.len() as f64;
let summat_xy = x_set.iter().zip(y_set.iter()).
fold(0f64,|acc,b|{
acc+(b.0*b.1)
});
let summat_x = x_set.iter().fold(0f64,|acc,b|{
acc+b
});
let summat_y = y_set.iter().fold(0f64,|acc,b|{
acc+b
});
let summat_xx = x_set.iter().fold(0f64,|acc,b|{
acc+b*b
});

let m = (n*(summat_xy) - 
summat_x*summat_y)/(n*summat_xx - summat_x.powi(2));
let b = (summat_y - m*summat_x)/n;
Self{
x_set,y_set,m,b
}
}

pub fn predict(&self,x_test:Vec<f64>)
->CDHResult<Vec<f64>>{
Ok(x_test.iter().
map(|&i| self.m*i + self.b).
collect::<Vec<f64>>())

}

}
