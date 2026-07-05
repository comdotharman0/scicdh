use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
use crate::probability::CDHResult;

// =========================================================================
// 1. DATA TRANSLATION AND RENDERING TRAITS
// =========================================================================

pub trait Transform<T> {
    /// Maps a transformation closure element-wise over the container.
    fn transform_to_collect<H>(&self, h: &H) -> CDHResult<Self>
    where
        H: Fn(T) -> T,
        Self: Sized;
}

pub trait Info {
    /// Renders a comprehensive statistical summary profile string.
    fn info(&self) -> CDHResult<String>;
}

// =========================================================================
// 2. THE ABSOLUTE NUMERIC ABSTRACTION INTERFACE
// =========================================================================

pub trait Numeric: 
    Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> +
    AddAssign + SubAssign + MulAssign + DivAssign +
    Copy + PartialEq + PartialOrd 
{
    fn zero() -> Self;
    fn one() -> Self;
    fn from_usize(val: usize) -> Self;
    fn from_f64(val: f64) -> Self;
    fn sqrt(self) -> Self;
    fn powi(self, exp: i32) -> Self;
    fn powf(self, exp: Self) -> Self;
    fn to_bits_u64(self) -> u64;
    fn from_bits_u64(bits: u64) -> Self;
}

// =========================================================================
// 3. HARDWARE-NATIVE BARE METAL IMPLEMENTATIONS
// =========================================================================

impl Numeric for f64 {
    #[inline] fn zero() -> Self { 0.0 }
    #[inline] fn one() -> Self { 1.0 }
    #[inline] fn from_usize(val: usize) -> Self { val as f64 }
    #[inline] fn from_f64(val: f64) -> Self { val }
    #[inline] fn sqrt(self) -> Self { self.sqrt() }
    #[inline] fn powi(self, exp: i32) -> Self { self.powi(exp) }
    #[inline] fn powf(self, exp: Self) -> Self { self.powf(exp) }
    #[inline] fn to_bits_u64(self) -> u64 { self.to_bits() }
    #[inline] fn from_bits_u64(bits: u64) -> Self { f64::from_bits(bits) }
}

impl Numeric for f32 {
    #[inline] fn zero() -> Self { 0.0f32 }
    #[inline] fn one() -> Self { 1.0f32 }
    #[inline] fn from_usize(val: usize) -> Self { val as f32 }
    #[inline] fn from_f64(val: f64) -> Self { val as f32 }
    #[inline] fn sqrt(self) -> Self { self.sqrt() }
    #[inline] fn powi(self, exp: i32) -> Self { self.powi(exp) }
    #[inline] fn powf(self, exp: Self) -> Self { self.powf(exp) }
    #[inline] fn to_bits_u64(self) -> u64 { self.to_bits() as u64 }
    #[inline] fn from_bits_u64(bits: u64) -> Self { f32::from_bits(bits as u32) }
}

// =========================================================================
// 4. FLOATING-POINT EMULATED INTEGER IMPLEMENTATIONS
// =========================================================================

impl Numeric for usize {
    #[inline] fn zero() -> Self { 0 }
    #[inline] fn one() -> Self { 1 }
    #[inline] fn from_usize(val: usize) -> Self { val }
    #[inline] fn from_f64(val: f64) -> Self { val as usize }
    #[inline] fn sqrt(self) -> Self { (self as f64).sqrt() as usize }
    #[inline] fn powi(self, exp: i32) -> Self { (self as f64).powi(exp) as usize }
    #[inline] fn powf(self, exp: Self) -> Self { (self as f64).powf(exp as f64) as usize }
    #[inline] fn to_bits_u64(self) -> u64 { self as u64 }
    #[inline] fn from_bits_u64(bits: u64) -> Self { bits as usize }
}

impl Numeric for isize {
    #[inline] fn zero() -> Self { 0 }
    #[inline] fn one() -> Self { 1 }
    #[inline] fn from_usize(val: usize) -> Self { val as isize }
    #[inline] fn from_f64(val: f64) -> Self { val as isize }
    #[inline] fn sqrt(self) -> Self { (self as f64).sqrt() as isize }
    #[inline] fn powi(self, exp: i32) -> Self { (self as f64).powi(exp) as isize }
    #[inline] fn powf(self, exp: Self) -> Self { (self as f64).powf(exp as f64) as isize }
    #[inline] fn to_bits_u64(self) -> u64 { self as u64 }
    #[inline] fn from_bits_u64(bits: u64) -> Self { bits as isize }
}
