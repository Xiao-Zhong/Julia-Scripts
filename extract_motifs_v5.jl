using DataFrames
using CSV

motifs_df = DataFrame(CSV.File("motif_for_1486_target2.csv", header=false, ignorerepeated=true))
#vscodedisplay(motifs_df)

# group by protein id
gdf = groupby(motifs_df, :Column1)

# extract last two motifs
gdf2 = [last(df, 2) for df in gdf]

# print out
for df in gdf2
    ids = df.Column1
    seqs = df.Column5
    for c in 1:length(ids)
        println(">", ids[c], "-$c", "\n", seqs[c])
    end
end