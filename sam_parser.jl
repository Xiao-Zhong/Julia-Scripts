using DataFrames
using Glob
using CSV

# count the number of lines in SAM file with a taxon
function taxa_count(taxa, sam_file)
    file = open(sam_file, "r")
    count = 0
    for line in readlines(file)
        startswith(line, "@SQ") && continue
        c = split(line)        
        if occursin(taxa, c[3]) 
            #println(s, "\t", line)
            count +=1
        end
    end
    return count
end

#generate a taxon, e.g., species, genus, family.
function taxa_term(term)
    names = []
    for t in term
        taxa = split(t, ";")
        push!(names, taxa[1], taxa[2], join(taxa[end-1:end], ";"))
    end
    return names
end

species = ["Pseudoalteromonadaceae;Pseudoalteromonas;spongiae", 
"Vibrionaceae;Vibrio;scophthalmi",
"Shewanellaceae;Shewanella;algae"]

#count the lines with a specific taxon in each SAM file, and output a DataFrame table
df = DataFrame()
df[!, :Sample] = [replace(f, "_L001.bwa.sam" => "") for f in glob("*.bwa.sam")]
for name in taxa_term(species)
    #println(name)
    counts = []
    for f in glob("*.bwa.sam")
        sc = taxa_count(name, f)
        push!(counts, sc)
    end

    name = replace(name, ";" => " ")
    df[!, name] = counts
    #println(name)
end

CSV.write("alignment_summary.csv", df)
#vscodedisplay(df)