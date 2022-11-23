using FASTX
using DataStructures

#input_motifs = DataFrame(CSV.File(IOBuffer(replace(read("motif_for_1486_target2.csv"), UInt8('\t') => UInt8(' '))), delim=" ", ignorerepeated=true))
## read "motif_for_1486_target.csv"
motif_dict = DefaultDict(Dict)
function motif_read(motif_input)
    f = open(motif_input, "r")
    lines = readlines(f)
    for line in lines
        columns = split(line)
        #println(target[1])
        current = [columns]
        if haskey(motif_dict, columns[1])
            current = motif_dict[columns[1]]
            push!(current, columns)
        end
        motif_dict[columns[1]] = current
    end
    #println(motif_dict)
end
#motif_read("motif_for_1486_target2.csv")

## read "targets.fa"
seq_dict = Dict{String,String}()
function fasta_read(fasta_input)
    reader = open(FASTA.Reader, fasta_input)
    for record in reader
        sequence = FASTA.sequence(record)
        id = FASTA.identifier(record)
        seq_dict[id] = sequence
        #println(">", id, "\n", seq_dict[id])
    end
    close(reader)
end
#fasta_read("targets.fa")

## extract sequences
function fasta_output()
    for gene in keys(motif_dict)
        #println(motif)
        #println(motif_dict[motif][end-1:end])
        num = 1;
        for motif in keys(motif_dict[gene][end-1:end])
            gene = motif_dict[gene][motif][1]
            start = parse(Int64, motif_dict[gene][motif][2])
            stop = parse(Int64, motif_dict[gene][motif][3])
            seq = SubString(seq_dict[gene], start, stop)
            println(">$gene-$num\t$start..$stop\n$seq")
            num +=1
        end
    end
end

function extract_seq(input_fasta, input_motif)
    fasta_read(input_fasta)
    motif_read(input_motif)
    fasta_output()
end

extract_seq(ARGS[1], ARGS[2])