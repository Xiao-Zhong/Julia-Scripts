using FASTX

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

motif_dict = Dict{String,Vector}()
function motif_read(motif_input)
    f = open(motif_input, "r")
    lines = readlines(f)[2:end]
    num = 1
    for line in lines
        target = split(line)[1]
        motifs = split(split(line)[3], ",")
        for motif in motifs
            start = split(motif, "..")[1]
            stop = split(motif, "..")[2]
            columns = [target, parse.(Int64,start), parse.(Int64,stop)]
            motif_dict["motif_$num"] = columns
            num +=1
        end
    end
    #println(motif_dict)
end
#gene_read("motifs.txt")

function fasta_output()
    for motif in keys(motif_dict)
        genome_id = motif_dict[motif][1]
        start = motif_dict[motif][2]
        stop = motif_dict[motif][3]
        seq = SubString(seq_dict[genome_id], start, stop)
        println(">$motif\t$genome_id\t$start..$stop\n$seq")
    end
end

function extract_seq(input_fasta, input_motif)
    fasta_read(input_fasta)
    motif_read(input_motif)
    fasta_output()
end

#extract_seq("targets.fa", "motifs.txt")
extract_seq(ARGS[1], ARGS[2])
