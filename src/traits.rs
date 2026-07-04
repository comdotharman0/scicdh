use crate::probability::CDHResult;
trait Transform<T=Self>{
pub fn transform_to_collect<H>(&self,h: &H)
->CDHResult<Self>
where
H: Fn(f64)->f64;

pub fn transform_to_value<H>(&self,h: &H)
->CDHResult<Self>
where
H: Fn(&[f64])->f64;


}

trait Info{

pub fn info(&self)->CDHResult<()>;


}
