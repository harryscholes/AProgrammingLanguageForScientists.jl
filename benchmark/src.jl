#=
Benchmarking code.
=#

using Statistics, BenchmarkTools, JSON

# Store execution time results
struct ExecutionTime
    mean::Float64
    std::Float64
end
function ExecutionTime(A::AbstractArray, f)
    x = @benchmark ($f)($A)
    ExecutionTime(mean(x.times)/10^9, std(x.times)/10^9)
end

plotter!(labels, timings, label, c) =
    bar!(labels, getfield.(timings, :mean), label=label, c=c,
        yerror=getfield.(timings, :std), marker=stroke(1.5, :black), α=.75)
