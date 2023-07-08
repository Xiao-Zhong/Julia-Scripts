using DataFrames, CSV

otu = CSV.File("Immuno_8/otu.csv") |> DataFrame
#otu[:,5] |> sum
samples = CSV.File("Immuno_8/samples2_Immuno8_v3.csv") |> DataFrame
otu2 = leftjoin(otu, samples[:, [:Sample, :Treatment]], on = :Column1 => :Sample)
otu2_gdf = groupby(otu2, :Treatment)

##sum reads per ASV per group/treatment
#otu2_gdf_sum = [combine(otu2_gdf, name => sum) for name in names(otu2_gdf)[2:end-1]]
otu_sum = combine(otu2_gdf, names(otu2_gdf)[2:end-1] .=> sum) 
#otu_sum[:, 5] |> sum

##edit multiple column names by pattern match
rename_foo(x) = replace(x, "_sum" => "")
rename!(rename_foo, otu_sum)
#names(otu_sum)

#otu_sum.Treatment = replace.(otu_sum.Treatment, r" $" => "")
otu_sum.Treatment = strip.(otu_sum.Treatment)
#otu_sum.Treatment
CSV.write("otu_sum_test.csv", otu_sum)
