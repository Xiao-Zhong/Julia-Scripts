using Glob
using DataFrames
using CSV

df = DataFrame( READ_PAIR_OK = [], 
                TOTAL_READ_PAIRS = [], 
                DROPPED_PAIR_STRANDS_MISMATCH = [],
                StrandTest_frFirstStrand = [],
                StrandTest_frSecondStrand = [],
                StrandTest_ambig_genesFountOnBothStrands  = [],
                StrandTest_ambig_noGenes = [],
                StrandTest_ambig_other = [],
                iLab_ID = [],
                Project = [],
                Sample = []
)

fs = glob("*/*/*/*.summary.txt")
for f in fs 
    fo = open(f, "r")

    numbers = []
    for line in readlines(fo)
        columns = split(line)

        columns[1] == "READ_PAIR_OK" && push!(numbers, columns[2])
        columns[1] == "TOTAL_READ_PAIRS" && push!(numbers, columns[2])
        columns[1] == "DROPPED_PAIR_STRANDS_MISMATCH" && push!(numbers, columns[2])
        columns[1] == "StrandTest_frFirstStrand" && push!(numbers, columns[2])
        columns[1] == "StrandTest_frSecondStrand" && push!(numbers, columns[2])
        columns[1] == "StrandTest_ambig_genesFountOnBothStrands" && push!(numbers, columns[2])
        columns[1] == "StrandTest_ambig_noGenes" && push!(numbers, columns[2])
        columns[1] == "StrandTest_ambig_other" && push!(numbers, columns[2])
    end
    close(fo)

    append!(numbers, split(f, "/")[1:3])
    #println(numbers)
    push!(df, numbers)
end

#vscodedisplay(df)
CSV.write("strandedness_summary.csv", df)