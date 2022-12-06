using DataFrames, CSV
using FASTX

df = DataFrame(CSV.File("Uploaded_sample_set-matrix-all-221206.tsv"))
df = df[!, Not(r"^X|^Zymo")]

replace!(df[!,5], "NA" => "0")
df[!,5] = parse.(Int, df[!,5])
sort!(df, :5, rev=true)
filter!(:"neg.control_S8.cladeReads" => x -> x > 0, df)

first_name = [split(x)[1] for x in df.name]
first_name = replace.(first_name, r"\[|\]" => "")
#vscodedisplay(first_name)
#vscodedisplay(df)

# reader = open(FASTA.Reader, "silva_nr99_v138.1_wSpecies_train_set2.fa")
# for record in reader
#     seq  = FASTA.sequence(record)
#     id = FASTA.identifier(record)
#     #println(">$num;$id\n$seq")
#     taxas = split(id, ";")
#     tmp = 0
#     for taxa in taxas
#         if taxa in first_name
#             #println(taxa, "\t", first_name)
#             tmp +=1
#             break
#         end
#     end
#     tmp == 0 && println(">$id\n$seq")
# end
# close(reader)

black_list = ["Coxiella","uncultured","Faecalibacterium","Cutibacterium","Acinetobacter","Methylophaga","Streptococcus","Lactobacillus","Corynebacterium"]

reader = open(FASTA.Reader, "silva_nr99_v138.1_wSpecies_train_set2.fa")
for record in reader
    seq  = FASTA.sequence(record)
    id = FASTA.identifier(record)
    #println(">$num;$id\n$seq")
    taxas = split(id, ";")
    tmp = 0
    for taxa in taxas
        if taxa in black_list
            #println(taxa, "\t", first_name)
            tmp +=1
            break
        end
    end
    tmp == 0 && println(">$id\n$seq")
end
close(reader)
