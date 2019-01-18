#=
Array-summation microbenchmarks.

Refer to `README.md` for installation instructions for the requirments.
=#

cd(@__DIR__)

using AProgrammingLanguageForScientists, BenchmarkTools, Statistics, PyCall, Plots

@pyimport numpy as np

# Use plotly backend
plotlyjs()

# Set benchmarking defaults to ensure 10,000 samples are always collected
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10_000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 999

# Storing execution time results
struct ExecutionTime
    mean::Float64
    std::Float64
end

ExecutionTime(x::BenchmarkTools.Trial) = ExecutionTime(mean(x.times)/1000,std(x.times)/1000)

macro executiontime(ex)
    quote
        ExecutionTime(@benchmark ($ex)($A))
    end
end

# Benchmarking

A = rand(1000,1000)

bsum = @executiontime sum
bnaivesum = @executiontime naivesum
bimprovedsum = @executiontime improvedsum
bnumpysum = @executiontime np.sum

timings = [bsum, bnaivesum, bimprovedsum, bnumpysum]

bar(["sum", "naivesum", "improvedsum", "numpy.sum"], getfield.(timings, :mean),
    yerror=getfield.(timings, :std), marker=stroke(2, :black),
    linewidth=0, c=:grey80,
    ylabel="Mean execution time ± SD (μs)", guidefontsize=12,
    tickfontsize=10, xtickfont=font(12, "Courier New"),
    grid=:y, legend=false, size=(600,400), dpi=300)

savefig("timings.png")
