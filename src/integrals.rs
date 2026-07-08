use crate::statistics::DataSet;
use crate::probability:: CDHResult;
use crate::traits::{Numeric,Transform};
pub struct Partition<const N:usize,T>
where
T: Numeric{
interval:[T;N],
n: usize
}

pub struct Function<const N: usize,T,F>
where 
F: Fn(f64)->f64,
T: Numeric{
func: F,
partition:Partition<N,T>
} 

pub fn integrals_checking()
-> CDHResult<()>{
let d = DataSet::new(
vec![0f64,1f64,2f64]);
println!("d={:?}",d);
println!("d²= {:?}",d.transform(&|i| i*i));
Ok(())

}

impl<T:Numeric+std::cmp::PartialEq> Partition<2,T>{
pub fn new(interval:[T;2],n:usize)->CDHResult<Self>{
if n<2{
return Err("n cannot be less than 2".to_string());
}
if interval[1]<interval[0]{
return Err("a should be smaller than b".to_string());
}
Ok(Self{interval,n})

}

pub fn create(&self)->CDHResult<DataSet<f64>>{
let a= self.interval[0].to_f64();
let b= self.interval[0].to_f64();
let diff = b-a;
let steps = diff/(self.n as f64);
let mut ds = DataSet::new(vec![a]);
let mut i = a+steps;
while i<b{
ds.push(i);
i+= steps;
}
Ok(ds)
}

}


impl<T,F> Function<2,T,F>
where
F: Fn(f64)->f64,
T: Numeric+ std::cmp::PartialEq{
pub fn new(func:F,partition:Partition<2,T>)->CDHResult<Self>{
if partition.interval[1]<partition.interval[0]{
return Err("a should be smaller than b".to_string());
}
Ok(Self{func,partition})
}

pub fn infimum(&self)->CDHResult<f64>{
Ok(0f64)
}


}
