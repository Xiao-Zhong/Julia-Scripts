using FASTX # read fasta
using BioSequences # reverse_complement; 
#using DataStructures

seq_dict = Dict{String,String}()
function fasta_read(fasta_input)
    reader = open(FASTA.Reader, fasta_input)
    for record in reader
        sequence = FASTA.sequence(record)
        id = FASTA.identifier(record)
        #println(id)
        seq_dict[id] = sequence
    end
    close(reader)
end

gene_dict = Dict()
function gene_read(gene_input)
    f = open(gene_input, "r")
    lines = readlines(f)

    c = 1
    for line in lines
        #println(line)
        columns = split(line)
        #println(columns)
        gene_dict["copy_$c"] = columns
        c += 1
    end
end

function fasta_output()
    for gene in keys(gene_dict)
        #println(gene)
        seq_id = gene_dict[gene][2]
        start = parse.(Int64, gene_dict[gene][9])
        stop = parse.(Int64, gene_dict[gene][10])
        strand = (start < stop ) ? "+" : "-"
        start > stop && ((start, stop) = (stop, start))

        seq = SubString(seq_dict[seq_id], start, stop)
        strand == "-" && (seq = reverse_complement(LongDNA{4}(seq)))
        println(">$gene\t$seq_id|$start..$stop|$strand\t$(length(seq))\n$seq")
    end
end

function extract_seq(input_fasta, input_gene)
    fasta_read(input_fasta)
    gene_read(input_gene)
    fasta_output()
end

#extract_seq(ARGS[1], ARGS[2])