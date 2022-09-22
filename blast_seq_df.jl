using FASTX # read fasta
using CSV
using DataFrames
using BioSequences # reverse_complement; 

function fasta_read(infasta)
    reader = open(FASTA.Reader, infasta)
    #ids = String[]; seqs = String[]
    df = DataFrame(id = String[], seq = String[])
    for record in reader
        sequence = FASTA.sequence(record)
        id = FASTA.identifier(record)
        #push!(ids, id) 
        #push!(seqs, sequence)
        push!(df, [id, sequence])
    end
    close(reader)

    return df
end

function gene_read(ingene)
    return DataFrame(CSV.File(ingene, header=false))
end

function fasta_output(infasta_file, ingene_file)
    fasta_df = fasta_read(infasta_file)
    gene_df = gene_read(ingene_file)

    for n in 1:nrow(gene_df)
        geneid = "copy$n"
        start = gene_df[n, 9]
        stop = gene_df[n, 10]
        strand = (start < stop ) ? "+" : "-"
        start > stop && ((start, stop) = (stop, start))

        seqdf = filter(:id => x -> x == gene_df[n, 2], fasta_df)
        (seq, seqid) = (seqdf.seq[1], seqdf.id[1])

        subseq = LongDNA{4}(SubString(seq, start, stop))
        strand == "-" && (reverse_complement!(subseq))
        println(">$geneid\t$seqid|$start..$stop|$strand\t$(length(subseq))\n$subseq")
    end
end

fasta_output(ARGS[1], ARGS[2])