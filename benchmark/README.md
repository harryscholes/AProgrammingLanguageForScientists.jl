# Benchmarking

## Installation

Julia provides a simple way to install reproducible virtual environments via its package manager and `Project.toml` and `Manifest.toml` files.

To install a virtual environment to run the array-summation benchmarks, change to this directory in the terminal, start the Julia REPL, enter the package manager using `]` key and then follow these steps:

```julia-repl
julia> ]

(v1.0) pkg> activate .

(benchmark) pkg> instantiate
```

## Running the benchmarks

To run the benchmarks:

```julia-repl
julia> include("speed.jl")

```
