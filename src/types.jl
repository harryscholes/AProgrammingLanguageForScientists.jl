#=
Types

Examples of how to define abstract and concrete types using biological sequences. Where appropriate, metaprogramming is used to generate repetitive code.
=#

#=
Abstract types.
=#

abstract type AbstractSequence end
abstract type AbstractNucleicAcid <: AbstractSequence end
abstract type AbstractRNA <: AbstractNucleicAcid end

#=
Concrete types.
=#

# Define the composite type `Protein` explicitly

"""
    Protein(sequence[, name, identifier])

A type representing a protein.
"""
struct Protein <: AbstractSequence
    sequence::String
    name::Union{String,Nothing}
    identifier::Union{String,Nothing}

    # Inner constructor to check whether a valid sequence was provided
    function Protein(s::String, n::Union{String,Nothing}, i::Union{String,Nothing})
        s = uppercase(s)
        Set(s) ≤ Protein_alphabet || throw(SequenceError())
        new(s, n, i)
    end
end

# Define nucleic acid composite types using metaprogramming and the `@sequence` macro

"""
    @sequence(T, U, Σ)

Define a sequence type `T` as a subtype of `U`, with alphabet `Σ`.
"""
macro sequence(T, U, Σ)
    local t = "$(Meta.quot(T))"[2:end]
    quote # Returns code expression for a new composite type (complete with docstring!)
        """
            $($t)(sequence [, name, identifier])

        $($t) sequence type.
        """
        struct $T <: $U
            sequence::String
            name::Union{String,Nothing}
            identifier::Union{String,Nothing}
            function ($T)(s::String, n::Union{String,Nothing}, i::Union{String,Nothing})
                s = uppercase(s)
                Set(s) ≤ $Σ || throw(SequenceError())
                new(s, n, i)
            end
        end
    end
end

# Define DNA and RNA types

@sequence DNA AbstractNucleicAcid DNA_alphabet

for T = (:mRNA, :tRNA, :rRNA)
    @eval @sequence $T AbstractRNA RNA_alphabet
end

#=
Methods on subtypes of `AbstractSequence`.
=#

for T = (:Protein, :DNA, :mRNA, :tRNA, :rRNA)
    @eval begin
        # Methods to get fields (in Julia, it is better to access fields using such methods,
        # rather than directly accessing the field)
        sequence(x::$T) = x.sequence
        name(x::$T) = x.name
        identifier(x::$T) = x.identifier

        # Alternative constructor that only accepts a sequence and sets `name` and
        # `identifier` to `nothing`
        ($T)(sequence::String) = ($T)(sequence, nothing, nothing)

        # String literal
        """
            $($T)\"<seq>\"

        $($T) string literal.
        """
        macro $(Symbol("$(T)_str"))(seq)
            quote
                ($$T)($seq)
            end
        end
    end
end

#=
Built-in functions can be extended with new methods. The syntax is either
```
import Base: foo
foo(...) = ...
```
or
```
Base.foo(...) = ...
```
=#

Base.length(x::AbstractSequence) = length(sequence(x))
Base.getindex(x::AbstractSequence, i) = getindex(sequence(x), i)
Base.lastindex(x::AbstractSequence) = lastindex(sequence(x))
Base.string(x::AbstractSequence) = sequence(x)
Base.count(x::AbstractSequence, b::Char) = count(i->i == b, sequence(x))
Base.count(x::AbstractSequence) = Dict(b=>count(x, b) for b = alphabet(typeof(x)))

#=
The type annotation `Type{}` means that `alphabet` accepts a type, rather than an instance.

Usage:
```
alphabet(Protein) # The type, rather than
alphabet(Protein("ACGHT")) # An instance
```
=#

"""
    alphabet(::Type{<:AbstractSequence})

The alphabet of a sequence type.
"""
alphabet(::Type{Protein}) = Protein_alphabet
alphabet(::Type{DNA}) = DNA_alphabet
alphabet(::Type{<:AbstractRNA}) = RNA_alphabet

"""
    alphabetsize(::Type{<:AbstractSequence})

The alphabet size of a sequence type.
"""
alphabetsize(x::Type{<:AbstractNucleicAcid}) = 4  # DNA/RNA
alphabetsize(x::Type{Protein}) = 20

"""
    gccontent(x::AbstractNucleicAcid)

GC content of a nucleic acid sequence.
"""
gccontent(x::AbstractNucleicAcid) = (count(x, 'G') + count(x, 'C')) / length(x)

"""
    complementtable(::Type{AbstractNucleicAcid})

Complement table of a nucleic acid type.
"""
complementtable(::Type{DNA}) = DNA_complement
complementtable(::Type{<:AbstractRNA}) = RNA_complement

# Example of how to use types to achieve certain behaviours

abstract type ComplementDirection end
struct ForwardComplement <: ComplementDirection end
struct ReverseComplement <: ComplementDirection end

function _complement(dir::Type{<:ComplementDirection}, x::T) where T<:AbstractNucleicAcid
    d = complementtable(T)
    range = dir == ForwardComplement ? (1:length(x)) : (length(x):-1:1)
    join([d[x[i]] for i = range])
end

"""
    complement(x::AbstractNucleicAcid)

Complement of a nucleic acid sequence.
"""
function complement(x::T) where T<:AbstractNucleicAcid
    T(_complement(ForwardComplement, x), name(x), identifier(x))
end

"""
    reversecomplement(x::AbstractNucleicAcid)

Reverse complement of a nucleic acid sequence.
"""
function reversecomplement(x::T) where T<:AbstractNucleicAcid
    T(_complement(ReverseComplement, x), name(x), identifier(x))
end

"""
    transcribe([T::Type{<:AbstractRNA},] x::DNA)

Transcribe a `DNA` sequence `x` to a type `T<:AbstractRNA`. If `T` is not provided, `mRNA`
is used by default.
"""
function transcribe(T::Type{<:AbstractRNA}, x::DNA)
    T(replace(sequence(x), "T"=>"U"), name(x), identifier(x))
end

# Method that defaults to transcribing `DNA` to `mRNA`
transcribe(x::DNA) = transcribe(mRNA, x)

"""
    translate(x::mRNA)
    translate(x::DNA)

Translate a sequence `x` to a `Protein` sequence. If `x` is a `DNA` type, the sequence is
first transcribed to `mRNA`.
"""
function translate(x::mRNA)
    l = length(x)
    l -= l % 3  # Full codons only
    seq = join(map(i->codontable[x[i:i+2]], 1:3:l))

    # Check for stop codons
    stop = findfirst(x->x == '*', seq)

    if stop !== nothing
        seq = seq[1:stop-1]
    end

    Protein(seq, name(x), identifier(x))
end

# Transcribe `DNA` to `mRNA` automagically before translating
translate(x::DNA) = translate(transcribe(x))

#=
Constants
=#

struct SequenceError <: Exception end

Base.show(io::IO, e::SequenceError) = print(io, "SequenceError: invalid sequence")

# Alphabet sets
const Protein_alphabet = Set("ACDEFGHIKLMNPQRSTVWY")
const DNA_alphabet = Set("ACGT")
const RNA_alphabet = Set("ACGU")

# Mapping tables
const DNA_complement = Dict('A'=>'T','T'=>'A','C'=>'G','G'=>'C')
const RNA_complement = Dict('A'=>'U','U'=>'A','C'=>'G','G'=>'C')

const codontable = Dict(
    "UUU" => "F", "CUU" => "L", "AUU" => "I", "GUU" => "V", "UUC" => "F", "CUC" => "L",
    "AUC" => "I", "GUC" => "V", "UUA" => "L", "CUA" => "L", "AUA" => "I", "GUA" => "V",
    "UUG" => "L", "CUG" => "L", "AUG" => "M", "GUG" => "V", "UCU" => "S", "CCU" => "P",
    "ACU" => "T", "GCU" => "A", "UCC" => "S", "CCC" => "P", "ACC" => "T", "GCC" => "A",
    "UCA" => "S", "CCA" => "P", "ACA" => "T", "GCA" => "A", "UCG" => "S", "CCG" => "P",
    "ACG" => "T", "GCG" => "A", "UAU" => "Y", "CAU" => "H", "AAU" => "N", "GAU" => "D",
    "UAC" => "Y", "CAC" => "H", "AAC" => "N", "GAC" => "D", "UAA" => "*", "CAA" => "Q",
    "AAA" => "K", "GAA" => "E", "UAG" => "*", "CAG" => "Q", "AAG" => "K", "GAG" => "E",
    "UGU" => "C", "CGU" => "R", "AGU" => "S", "GGU" => "G", "UGC" => "C", "CGC" => "R",
    "AGC" => "S", "GGC" => "G", "UGA" => "*", "CGA" => "R", "AGA" => "R", "GGA" => "G",
    "UGG" => "W", "CGG" => "R", "AGG" => "R", "GGG" => "G")
