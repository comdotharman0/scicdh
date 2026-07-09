use crate::statistics::DataSet;
use crate::probability:: CDHResult;
use crate::traits::{Numeric,Transform};
#[derive(Debug)]
pub struct Partition<const N:usize,T>
where
T: Numeric{
interval:[T;N],
n: usize
}

#[derive(Debug)]
pub struct Function<const N: usize,T,F>
where 
F: Fn(f64)->f64,
T: Numeric{
func: F,
partition:Partition<N,T>
} 

pub fn integrals_checking()
-> CDHResult<()>{
let p = Partition::new([1f64,11f64],10)?;
let f = Function::new(|a| a*a,p)?;
//println!("f={:#?}",f.partition.create()?);
//println!("d²= {:?}",d.transform(&|i| i*i));
println!("infimum={:?}, supremum={:?}",f.infimum()?,
f.supremum()?);
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
let b= self.interval[1].to_f64();
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
self.partition.create()?
.to_f64()?.transform(&self.func)?
.min()
}

pub fn supremum(&self)->CDHResult<f64>{
self.partition.create()?
.to_f64()?.transform(&self.func)?
.max()
}


}
