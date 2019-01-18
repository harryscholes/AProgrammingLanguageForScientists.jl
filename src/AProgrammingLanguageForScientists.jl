module AProgrammingLanguageForScientists

using DataStructures

include("types.jl")

export
    AbstractSequence,
    AbstractNucleicAcid,
    AbstractRNA,
    Protein,
    DNA,
    mRNA,
    tRNA,
    rRNA,
    @Protein_str,
    @DNA_str,
    @mRNA_str,
    @tRNA_str,
    @rRNA_str,
    @sequence,
    sequence,
    name,
    identifier,
    alphabet,
    alphabetsize,
    gccontent,
    complement,
    reversecomplement,
    transcribe,
    translate,
    Protein_alphabet,
    DNA_alphabet,
    RNA_alphabet

include("arrays.jl")

export
    KmerVector,
    sequence,
    kay,
    hashes,
    seed,
    MinHashSketch,
    BottomKSketch,
    jaccard

include("speed.jl")

export
    naivesum,
    improvedsum

end
