use std::borrow::Cow;
use crate::probability::CDHResult;
pub enum ArithmeticProgression<T>{
FromSet{
set:Vec<T>
},
FromValues{
a:T,d:T
}

}
pub fn check_series()->CDHResult<()>{
let a = ArithmeticProgression::FromSet{
set : vec![0f64,1f64,2f64]
};
if let ArithmeticProgression::FromSet{set:n}= a{
println!("n ==== {:?}",n);
}
Ok(())
}

type AP<T> = ArithmeticProgression<T>;

impl ArithmeticProgression<f64>{
pub fn get_setvalues<'a>(&'a self)
->CDHResult<Cow<'a,[f64]>>{
match self{
AP::FromSet{set} => Ok(Cow::Borrowed(set.as_slice())),
AP::FromValues{a,d}=>{
Ok(Cow::Owned((1..=10).map(|n| a+(n as f64-1f64)*d)
.collect::<Vec<f64>>()))
}
}
}


}
