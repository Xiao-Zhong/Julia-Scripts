using DataFrames, CSV
using RCall
#using StatsBase

##read dada2 output files, sample table
function input2df(asv_file, taxa_file, sample_file)
    asv = CSV.File(asv_file) |> DataFrame;
    taxa = CSV.File(taxa_file) |> DataFrame;
    taxa = permutedims(taxa, 1)
    sample = CSV.File(sample_file) |> DataFrame;
    sample.Treatment = strip.(sample.Treatment)
    return asv, taxa, sample
end

##merge these table above to be a metadata table
function df_merger(asv, taxa, sample)
    asv_taxa = vcat(asv, taxa)
    new_asv_id = ["ASV-$(i-1)" for i in 2:ncol(asv_taxa)]
    rename!(asv_taxa, names(asv_taxa)[2:end] .=> new_asv_id)
    leftjoin!(asv_taxa, sample, on = "Column1" => "Sample Name")
    #unique(asv_taxa.Treatment)
    #remove taxa rows or with "missing" in the last few columns
    df_taxa = filter(:Treatment => ismissing, asv_taxa)
    #extract the rows with a treatment targeted.
    filter!(:Treatment => !ismissing, asv_taxa)
    #rename!(asv_taxa, Dict(:Column1 => "Sample"))
    return asv_taxa, df_taxa
end

##filter ASV with few read support and existing in few samples
function asv_filter(asv, treatment)
    filter!(:Treatment => x -> x ∈ treatment, asv)
    #CSV.write("asv.csv", asv)
    group_label = length(treatment) > 1 ? "Treatment" : "Column1"
    samples_cutoff = length(treatment) > 1 ? 2 : 1
    println("groupby: $group_label")

    ##split-apply-combine: filter ASVs within a treatment
    gdfs = groupby(asv, "$(group_label)")
    group_ordered = []; gdfs_filter = []
    for gdf in gdfs
        push!(group_ordered, unique(gdf[:, "$(group_label)"]))
        #println("before\t", size(gdf))
        select_c = [1,]
        for n in 2:ncol(gdf)-3
            c = (gdf[:, n] .> 5 ) |> count # more than 5 reads per sample, two samples at least per treatment
            (c >= samples_cutoff) && push!(select_c, n)
        end
        gdf = select(gdf, select_c)
        #println("after\t", size(gdf))
        push!(gdfs_filter, gdf)
    end

    return gdfs_filter, group_ordered
end

##extract good ASV based on any taxonomy classifcation and term
function asv_taxa_extract(asv_in, taxa_in, sample_in, treatment)
    (asv_taxa, df_taxa) = df_merger(asv_in, taxa_in, sample_in)
    #vscodedisplay(df_taxa)
    #CSV.write("taxa.csv", df_taxa)

    (asv_filter_dfs, groups) = asv_filter(asv_taxa, treatment)
    println("#Groups:\n", groups)
    println("#Number\tASVs")
    #vscodedisplay(asv_filter_dfs[1])
    #CSV.write("asv_filter.csv", asv_filter_dfs[1])

    gdfs_filter_term = [vcat(gdf, df_taxa, cols=:intersect) for gdf in asv_filter_dfs]
    #vscodedisplay(gdfs_filter_term[1])
    #CSV.write("gdfs_filter_term.csv", gdfs_filter_term[1])

    ## to get the index of one rank
    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    rank_i = findfirst(rank .== ranks) #faster and efficient

    ##extract ASVs assinged to a taxonomy term.
    asv_list = Vector{String}[]; #asv_list_c = Vector{String}[]
    for gdf in gdfs_filter_term
        asv_detected = String[]; #asv_detected_c = String[]
        taxa_detected = collect(gdf[end-7+rank_i, :])
        # for n in 2:length(taxa_detected)
        #     if term == taxa_detected[n]
        #         push!(asv_detected, names(gdf)[n])
        #         #[append!(asv_detected_c, repeat([names(gdf)[n],], reads)) for reads in gdf[1:end-7, n]]
        #     end
        # end
        asv_detected = ifelse.(taxa_detected[2:end] .== term, names(gdf)[2:end], missing)
        asv_detected = filter(!ismissing, asv_detected)
        println(length(asv_detected), "\t", asv_detected)
        
        push!(asv_list, asv_detected)
        #push!(asv_list_c, asv_detected_c)
    end
    #asv_list
    return asv_list, groups
end

function venn_plot(asv, category, category_label_size, output_file)
    ##plot veen diagram using R
    R"""
    library(ggVennDiagram); library(ggplot2)

    #png("mtcars.png")
    ggVennDiagram($asv, 
        label_alpha = 0, 
        label_percent_digit = 1,
        category.names = $category,
        edge_size = 0.5,
        set_size = $category_label_size,
    ) + 
    ggplot2::scale_fill_gradient(low="white",high = "red") +
    scale_x_continuous(expand = expansion(mult = c(.11)))

    ggsave($output_file)
    #dev.off()
    """
end

## to check ASV at one taxonomy level and with a term
rank = "Kingdom"
term = "Bacteria"

## input files
asv_in_file = "Immuno_8/dada2_asv_Immuno_8.csv"
taxa_in_file = "Immuno_8/dada2_taxa_names_Immuno_8.csv"
sample_in_file =  "Immuno_8/samples2_Immuno8_v3.csv"
## run functions
(asv_in, taxa_in, sample_in) = input2df(asv_in_file, taxa_in_file, sample_in_file)

#treatment = ["FM control", "BSF 10", "BSF 30"]
#treatment = ["FM control + Immuno", "BSF 10 + Immuno", "BSF 30 + Immuno"]
#treatment = ["FM control"]
#treatment = ["BSF 10"]
#treatment = ["BSF 30"]
#treatment = ["FM control + Immuno"]
#treatment = ["BSF 10 + Immuno"]
treatment = ["BSF 30 + Immuno"]

(asv_p, category_p) = asv_taxa_extract(asv_in, taxa_in, sample_in, treatment);

if length(asv_p) >= 2
    venn_plot(asv_p, category_p, 3, "BSF30_Immuno-samples-only.pdf")
else
    println("Cannot find overlapping ASV between any two groups!")
end
