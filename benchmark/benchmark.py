"""
Python array-summation benchmarks.

Refer to `README.md` for installation instructions for the requirments.
"""
import json
import numpy as np


def saveresult(f, res):
    results[f] = {"mean": res.average, "std": res.stdev}


results = {}

A = np.random.rand(1_000_000)

# numpy.sum
%timeit -o -n 1 -r 10_000 np.sum(A)
saveresult("numpy.sum", _)

A = list(A)

# sum
%timeit -o -n 1 -r 1_000 sum(A)
saveresult("sum", _)


# naivesumpy
def naivesum(A):
    accumulator = 0.

    for element in A:
        accumulator += element

    return accumulator


%timeit -o -n 1 -r 1_000 naivesum(A)
saveresult("naivesum", _)


with open("benchmarks_python.json", "w") as f:
    json.dump(results, f)
