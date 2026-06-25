pub struct DataSet(pub Vec<f64>);
impl DataSet{
pub fn mean(&self)->f64{
let len: f64 = self.0.len() as f64;
let sum: f64 = self.0.iter().sum();
sum/len
}


}


