#=
Speed

Array-summation methods. Refer to ../benchmark/speed.jl for microbenchmarks.
=#

"""
    naivesum(A)

Returns the sum of all elements in `A`.
"""
function naivesum(A::AbstractArray{T}) where T<:Number
    accumulator = zero(T)

    for element in A
        accumulator += element
    end

    return accumulator
end

"""
    improvedsum(A)

Returns the sum of all elements in `A`.
"""
function improvedsum(A::AbstractArray{T}) where T<:Number
    accumulator = zero(T)

    @inbounds @simd for element in A
        accumulator += element
    end

    return accumulator

    # @inbounds @simd for i = eachindex(A)
    #     accumulator += A[i]
    # end
    #
    # return accumulator
end
