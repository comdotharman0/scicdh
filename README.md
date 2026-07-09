# SciCDH Statistics & Probability Library

A pure Rust library for descriptive statistics and 
discrete probability theory. Built with **zero external 
dependencies** using only Rust's standard library.

[![DOI](https://zenodo.org/badge/1280434144.svg)](https://doi.org/10.5281/zenodo.21279300)

## Author
Harman Singh  
ORCID: [0009-0004-5218-9471](https://orcid.org/0009-0004-5218-9471)
## Features

### Statistics (`statistics.rs`)
- Mean, Median, Mode
- Range, Min, Max
- Sample & Population Variance
- Sample & Population Standard Deviation
- Coefficient of Variation
- Pearson Correlation Coefficient
- Skewness & Kurtosis
- k-th Central Moment
- Unique set extraction

### Probability (`probability.rs`)
- Classical Probability
- Factorial, Permutations (nPr), Combinations (nCr)
- CDF (Cumulative Distribution Function)
- Random Variable transformations

### Discrete Distributions
| Distribution | PMF | Mean | Variance |
|---|---|---|---|
| Binomial | C(n,x)·pˣ·(1-p)^(n-x) | np | np(1-p) |
| Geometric | p·(1-p)^(x-1) | 1/p | (1-p)/p² |
| Pascal | C(x-1,r-1)·pˣ·(1-p)^(x-r) | r/p | r/p² |
| HyperGeometric | C(T,x)·C(N-T,n-x)/C(N,n) | nT/N | see docs |

### Random Variables
- Discrete Random Variable: E[X], Var(X), σ(X)
- Moment Generating Function M(t) = E[e^(tX)]
- Expectation of transformations E[h(X)]
- Joint Probability Distributions
- Covariance and Correlation
- Conditional Expectations E[X1|X2], E[X2|X1]

## Usage

```rust
use probability::{Probability, RandomVariable, Binomial};

// Combinations
let c = Probability::n_c_r(5, 2).unwrap();
assert_eq!(c, 10.0);

// Random Variable — fair die
let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
let px = vec![1.0/6.0; 6];
let die = RandomVariable::new(x, px);
println!("E[die] = {}", die.mean().unwrap()); // 3.5

// Binomial — 10 coin flips
let b = Binomial::new(10, vec![0,1,2,3,4,5], 0.5);
println!("Mean = {}", b.mean().unwrap()); // 5.0
```

## Testing

```bash
cargo test
```

32 tests, 0 failures. All expected values mathematically
verified by hand before being written into tests.

## Documentation

```bash
cargo doc --no-deps --target-dir doc
```

Opens `doc/` folder. HTML documentation is generated
from inline `///` doc comments on every public function.

## Design Decisions

- **No external dependencies** — pure `std` only.
  Every algorithm is implemented from mathematical
  definitions, not wrapped from another crate.
- **`CDHResult<T>`** — all functions return
  `Result<T, String>` so errors are explicit and
  recoverable, never panics on valid input.
- **`f64` throughout** — sufficient precision for
  scientific computing at this scale.
- **`characteristic_func` deprecated** — the
  characteristic function E[e^(itX)] requires complex
  numbers which `f64` cannot represent. Previous
  implementation silently returned `NaN`. Now returns
  `Err` with a clear message until complex support
  is added.

## Limitations

- `factorial(n)` only supports n ≤ 101 due to f64
  precision limits beyond ~15 significant digits
- Complex number support not yet implemented
  (affects characteristic function)
- Currently supports discrete distributions only —
  continuous distributions (Normal, Poisson, etc.)
  are planned

## References

The mathematical specifications, distributions, and foundational statistical algorithms implemented in this library are sourced from the following academic literature:

* Axler, S. (2020). *Measure, Integration & Real Analysis* (Graduate Texts in Mathematics, Vol. 280). Springer.
* Hines, W. W., Montgomery, D. C., Goldsman, D. M., & Borror, C. M. (2003). *Probability and Statistics in Engineering* (4th ed.). John Wiley & Sons.


## License

This project is licensed under the MIT License. See the full [LICENSE](https://github.com/comdotharman0/scicdh/blob/main/LICENSE) file for copyright details and permissions.

