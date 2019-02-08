#=
Julia array-summation benchmarks.

Refer to `README.md` for installation instructions for the requirments.
=#

cd(@__DIR__)

using AProgrammingLanguageForScientists, PyCall

@pyimport numpy as np

include("src.jl")

# Set benchmarking defaults to ensure 10,000 samples are always collected
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10_000

results = Dict()

A = rand(1_000_000)

results["sum"] = ExecutionTime(A, sum)
results["naivesum"] = ExecutionTime(A, naivesum)
results["improvedsum"] = ExecutionTime(A, improvedsum)
results["numpy.sum"] = ExecutionTime(A, np.sum)

write("benchmarks_julia.json", JSON.json(results))
