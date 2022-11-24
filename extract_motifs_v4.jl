## read "motif_for_1486_target.csv"
motif_dict = Dict()
function motif_read(motif_input)
    f = open(motif_input, "r")
    lines = readlines(f)
    for line in lines
        columns = split(line)
        if haskey(motif_dict, columns[1])
            push!(motif_dict[columns[1]], columns)
        else
            motif_dict[columns[1]] = [columns]
        end
    end
    #println(motif_dict[columns[1]])
end
#motif_read("motif_for_1486_target2.csv")

## extract sequences
function fasta_output()
    for gene in keys(motif_dict)
        #println(motif_dict[gene])
        for motif in keys(motif_dict[gene][end-1:end])
            start = parse(Int64, motif_dict[gene][motif][2])
            stop = parse(Int64, motif_dict[gene][motif][3])
            seq = motif_dict[gene][motif][5]
            println(">$gene-$motif\t$start..$stop\n$seq")
        end
    end
end

motif_read("motif_for_1486_target2.csv")
fasta_output()