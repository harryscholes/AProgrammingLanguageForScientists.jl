using AProgrammingLanguageForScientists, Test
using AProgrammingLanguageForScientists: SequenceError, DNA_complement, RNA_complement, complementtable

const n = "name"
const id = "id"

@testset "$T" for T = (Protein, DNA, mRNA, rRNA, tRNA)
    @testset for l = (10, 100, 1000), repeat = 1:10
        seq = join(rand(alphabet(T), l))
        x = T(seq, n, id)

        @testset "Getters" begin
            @test sequence(x) == seq
            @test name(x) == n
            @test identifier(x) == id
        end
        @testset "Equality" begin
            @test x == x
            @test string(x) == seq
            @test length(x) == l == lastindex(x)
        end
        @testset "Slicing" begin
            for idx = (l+1, 0, -1)
                @test_throws BoundsError x[idx]
            end
            for i = rand(1:l-5, 10), j = rand(1:5, 10)
                @test x[i:i+j] == seq[i:i+j]
            end
        end
        @testset "Counting" begin
            @test count(x, 'A') == count(x -> x == 'A', seq)
            @test sum(values(count(x))) == l
            d = count(x)
            for k = keys(d)
                k ∈ seq && (@test d[k] == count(x -> x == k, seq))
                k ∉ seq && (@test d[k] == 0)
            end
        end
    end
    @testset "SequenceError" begin
        @test_throws SequenceError T("zzz")
        @test_throws SequenceError T("zzz", n, id)
    end

    @test Protein"ACDEFGHIKLMNPQRSTVWY" == Protein("ACDEFGHIKLMNPQRSTVWY")
    @test DNA"ACTG" == DNA("ACTG")
    @test rRNA"ACUG" == rRNA("ACUG")
    @test mRNA"ACUG" == mRNA("ACUG")
    @test tRNA"ACUG" == tRNA("ACUG")
end

@testset "Alphabet" begin
    @testset "Protein" begin
        @test alphabet(Protein) == Protein_alphabet
        @test alphabetsize(Protein) == 20
    end
    @testset "DNA" begin
        @test alphabet(DNA) == DNA_alphabet
        d = complementtable(DNA)
        @test d == DNA_complement
        for (a,b) = (('A', 'T'), ('C', 'G'))
            @test d[a] == b
            @test d[b] == a
        end
        @test alphabetsize(DNA) == 4
        @test gccontent(DNA("GCAGCTCGACGT")) == 2/3
    end
    @testset "RNA" for T = (mRNA, rRNA, tRNA)
        @test alphabet(T) == RNA_alphabet
        d = complementtable(T)
        @test d == RNA_complement
        for (a,b) = (('A', 'U'), ('C', 'G'))
            @test d[a] == b
            @test d[b] == a
        end
        @test alphabetsize(T) == 4
    end
end

@testset "Complement" begin
    @testset for T = (DNA,mRNA,tRNA,rRNA), (seq,comp) = (("ACG","TGC"),("GCACAG","CGTGTC"))
        if T <: AbstractRNA
            seq = replace(seq, "T" => "U")
            comp = replace(comp, "T" => "U")
        end

        x = T(seq)

        @testset "Complement" begin
            y = complement(x)
            @test y isa T
            @test sequence(y) == comp
            @test complement(y) == x
            @test length(x) == length(y) == length(seq)
        end
        @testset "Reverse complement" begin
            y = reversecomplement(x)
            @test y isa T
            @test sequence(y) == reverse(comp)
            @test reversecomplement(y) == x
            @test T(reverse(sequence(complement(x)))) == y
            @test length(x) == length(y) == length(seq)
        end
    end
end

@testset "Transcribe" begin
    @testset for T = (mRNA, tRNA, rRNA), l = (10, 100, 1000), repeat = 1:10
        seq = join(rand(alphabet(DNA), l))
        x = DNA(seq)
        y = transcribe(T, x)
        @test y isa T
        @test sequence(y) == replace(seq, "T" => "U")
        @test y == T(replace(seq, "T" => "U"))
        @test length(x) == length(y) == l
    end
end

@testset "Translate" begin
    for (dna,protein) = (("GGTGGCACC","GGT"), ("TGGGTCCAG","WVQ"), # Regular
                         ("GGTGGCACCAA","GGT"), ("TGGGTCCAGA","WVQ"), # Partial codons
                         ("GCGGTTTGA","AV"), # Stop at end
                         ("GCGGTTTGAGCG", "AV") # Stop in middle
                        )
        x = DNA(dna, n, id)
        y = translate(x)
        @test sequence(y) == protein
        @test name(y) == n
        @test identifier(y) == id
    end
end

x = DNA"ACGT"
kv = KmerVector(x, 3)
IndexStyle(kv)

@testset "Arrays" begin
    @testset for T = (Protein, DNA, mRNA, rRNA, tRNA)
        @testset "Length $l $k" for l = (10, 100, 1000), k = rand(4:8, 10)
            seq = join(rand(alphabet(T), l))
            x = T(seq)
            y = KmerVector(x, k)
            @test IndexStyle(y) == IndexLinear()

            @testset "KmerVector" begin
                len = length(x) - k + 1
                @test y == y
                @test length(y) == len
                @test size(y) == (len,)
                @test sequence(y) == x
                @test kay(y) == k
                for i = rand(1:len, 10), j = rand(1:len, 10)
                    @test y[i] == x[i:i+k-1]
                    @test y[[i,j]] == [x[i:i+k-1],x[j:j+k-1]]
                end
            end

            @testset "MinHashSketch" begin
                z = MinHashSketch(y, k)
                @test z == z
                a = MinHashSketch(y, k)
                @test z == a
                for i = rand(1:k, 10), j = rand(1:k, 10)
                    @test z[i] == a[i]
                    @test z[[i,j]] == a[[i,j]]
                end
                @test length(z) == kay(z) == k == length(hashes(z))
                @test size(z) == (k,)
                @test IndexStyle(z) == IndexLinear()
            end

            @testset "BottomKSketch" begin
                s = 123
                if k > length(y)
                    @test_throws ArgumentError BottomKSketch(y, l)
                else
                    z = BottomKSketch(y, k, s)
                    @test z == z
                    a = BottomKSketch(y, k, s)
                    @test z == a
                    for i = rand(1:k, 10), j = rand(1:k, 10)
                        @test z[i] == a[i]
                        @test z[[i,j]] == a[[i,j]]
                    end
                    @test length(z) == kay(z) == k
                    @test seed(z) == s
                    @test_throws ArgumentError BottomKSketch(y, l+1)
                    @test length(z) == kay(z) == k == length(hashes(z))
                    @test size(z) == (k,)
                    @test IndexStyle(z) == IndexLinear()
                end
            end
        end
    end

    @testset "Jaccard" begin
        @test jaccard([],[]) == 1.0

        for _ = 1:10
            x = rand(100)
            @test jaccard(x,x) == 1.0

            y, z = rand(1:50, 100), rand(51:100, 100)
            @test jaccard(y, z) == 0.

            y, z = Set(y), Set(z)
            l = min(length(y), length(z))
            y, z = collect(y)[1:l], collect(z)[1:l]
            @test jaccard(y, [y;z]) == 0.5
        end
    end
end

@testset "Speed" begin
    for i = rand(1:100, 100), j = rand(1:100, 10)
        A = rand(i,j)
        t = sum(A)
        for f = (naivesum, improvedsum)
            u = f(A)
            @test t ≈ u
            @test typeof(u) == eltype(A)
        end
    end
end
