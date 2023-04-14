using DataFrames, CSV

function convert(asv_in)
    df = CSV.File(asv_in) |> DataFrame
    #names(df)

    fname = splitext(asv_in)[1]
    open("$(fname).fa", "w") do file
        for (i, seq) in enumerate(names(df)[2:end])
            write(file, ">ASV-$(i)\n$(seq)\n")
        end
    end

    df2 = permutedims(df)
    df2 = rename(df2, Symbol.(Vector(df2[1, :])))[2:end, :]

    asv_out = hcat(["ASV-$x" for x in 1:ncol(df)-1], df2)
    CSV.write("$(fname)_simplified.txt", asv_out; delim = "\t")
end

## input dada2 ASV table ##
convert("Immuno_8/dada2_asv_Immuno_8.csv")
convert("Immuno_9/dada2_asv_Immuno_9.csv")
convert("SCP/dada2_asv_SCP.csv")
