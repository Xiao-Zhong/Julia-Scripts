using DataFrames
using Glob
using CSV

# extract the number of reads assigned to a taxon in Kraken Report file
function taxa_count(taxa, kraken_file)
    count = 0
    file = open(kraken_file, "r")
    for line in readlines(file)
        c = split(line)
        #println(c[6])
        c[6] == taxa && (count = c[2])
    end
    return count
end
#taxa_count("Vibrio", "1-Larvae-1421_S1_L001_kraken2_report.txt")

#generate a taxon, e.g., species, genus, family.
function taxa_term(terms)
    names = []
    for term in terms
        taxa = split(term, ";")
        push!(names, taxa[1], taxa[2], join(taxa[end-1:end], " "))
    end
    return names
end

species = ["Pseudoalteromonadaceae;Pseudoalteromonas;spongiae", 
"Vibrionaceae;Vibrio;scophthalmi",
"Shewanellaceae;Shewanella;algae"]

#extract the number of reads assgined to a specific taxon in each Kraken file, and output a DataFrame table
df = DataFrame()
df[!, :Sample] = [replace(f, "_L001_kraken2_report.txt" => "") for f in glob("*_L001_kraken2_report.txt")]
for name in taxa_term(species)
    #println(name)
    counts = []
    for f in glob("*_L001_kraken2_report.txt")
        sc = taxa_count(name, f)
        push!(counts, sc)
    end

    df[!, name] = counts
end

CSV.write("kraken2_summary.csv", df)