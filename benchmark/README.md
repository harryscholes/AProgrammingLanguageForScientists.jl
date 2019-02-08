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

Run `benchmark.jl` and `benchmark.py`, then plot the results with `plot.jl`. Example benchmarking results are in `benchmark_julia.json` and `benchmark_python.json`, which are plotted in `timings.pdf`.
