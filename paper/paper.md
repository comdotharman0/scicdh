---
title: 'scicdh: A Zero-Dependency Pure Rust Library for Descriptive Statistics and Discrete Probability'
tags:
  - Rust
  - statistics
  - discrete probability
  - mathematics
  - no-external-dependencies
authors:
  - name: Harman Singh
    orcid: 0009-0004-5218-9471
    affiliation: Independent Researcher
date: 23 July 2026
bibliography: paper.bib
---

# Summary

`scicdh` is an open-source, lightweight Rust library designed for fast descriptive statistics and discrete probability modeling without relying on external third-party crates. The library provides self-contained mathematical abstractions for computing central tendencies, dispersion, and discrete distribution properties—including Binomial, Geometric, and Hypergeometric distributions—alongside custom probability mass functions (PMFs) and cumulative distribution functions (CDFs). `scicdh` is built to operate cleanly in resource-constrained environments, such as mobile terminal environments (e.g., Termux) and embedded systems, where minimal compile-time overhead and predictable binary footprints are essential.

# Statement of Need

In scientific computing and data science pipelines implemented in Rust, statistical operations are frequently offloaded to crates such as `statrs` or `nalgebra`. While feature-rich, these crates bring substantial transitively linked dependency trees (e.g., `num-traits`, `rand`, or platform-specific BLAS bindings). In low-resource development setups, offline environments, or embedded target architectures, large dependency trees increase compilation latency, expand binary bloat, and introduce external supply-chain vulnerabilities.

`scicdh` fills this gap by delivering a self-contained, pure-Rust implementation of foundational statistical and discrete probability routines with zero third-party dependencies. It allows researchers and software engineers to calculate key statistical metrics and evaluate discrete probability distributions natively using only Rust's standard numerical primitives.

# State of the Field

Existing numerical libraries in the Rust ecosystem typically fall into two categories: high-level numerical framework wrappers or heavy linear algebra suites. 

1. **`statrs`**: Provides comprehensive statistical distribution implementations but requires external traits (`num-traits`) and random sampling crates (`rand`), which introduce transitively linked dependencies.
2. **`nalgebra` / `ndarray`**: Optimized for multidimensional vector and matrix computations, but overly heavy when a project requires only discrete probability distribution calculations or basic descriptive statistics.

`scicdh` distinguishes itself by maintaining a strictly zero-dependency footprint while offering intuitive API access to core discrete distribution equations.

# Software Design & Mathematics

`scicdh` models discrete distributions using pure functional and iterative paradigms based on classical applied statistics [@hines2003probability] and rigorous analytical principles [@axler2020measure].

For a discrete random variable $X$ taking values in a finite or countable set with probability mass function $P(X = x)$, the expected value $E[X]$ and variance $\operatorname{Var}(X)$ are computed directly:

$$E[X] = \sum_{x} x \cdot P(X = x)$$

$$\operatorname{Var}(X) = E[X^2] - (E[X])^2$$

The Binomial, Geometric, and Hypergeometric modules provide exact evaluation for probability mass and cumulative distribution calculations using Rust's standard floating-point and integer math, avoiding external trait abstractions.

# Research Impact Statement

`scicdh` offers a lightweight, verifiable tool for researchers and developers needing deterministic statistical routines without external crate dependencies. By providing an offline-ready, standalone codebase, it enables educational modeling, embedded statistical tracking, and reproducible computations on minimal hardware configurations.

# AI Usage Disclosure

Generative AI tools were used solely to assist with code optimization during development. The core architecture, algorithms, and primary mathematical implementations—comprising over 98% of the library codebase—were authored independently.
# References
