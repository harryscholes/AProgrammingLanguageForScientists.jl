#=
Plot Julia and Python array-summation benchmarks.

Refer to `README.md` for installation instructions for the requirments.
=#

cd(@__DIR__)

using Plots, JSON

include("src.jl")

plotlyjs() # Use plotly backend

j = JSON.parsefile("benchmarks_julia.json")
p = JSON.parsefile("benchmarks_python.json")

timings = map(x->ExecutionTime(x["mean"], x["std"]),
    [j["sum"], j["naivesum"], j["improvedsum"], p["sum"], p["naivesum"], p["numpy.sum"]])

plot(linewidth=0, marker=stroke(2, :black),
    ylabel="Mean execution time ± SD (s)", guidefontsize=12,
    tickfontsize=10, xtickfont=font(12, "Courier New"),
    grid=:y, legend=true, size=(400,350), dpi=1000, xrot=45, yscale=:log10, yticks=(10.).^(-4:0),
    thickness_scaling=2)

plotter!(["sum", "naivesum", "improvedsum"], timings[1:3], "Julia", :grey75)
plotter!([" sum", " naivesum", "numpy.sum"], timings[4:end], "Python", :grey25)

savefig("timings.pdf")
