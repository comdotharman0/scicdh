pub struct Matrix<T>(pub Vec<Vec<T>>);

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