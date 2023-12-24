#Run `using Pkg; Pkg.status("XX")` to check version.
#Run `import Pkg; Pkg.add("XX")` to install the XX package.
using FASTX
using DataFrames

#Check if there are enough command-line arguments
if length(ARGS) < 3
    println("Usage: julia split_fasta_gtf.jl <fasta_file> <gtf_file> <fragment_size>")
    exit(1)  # Exit the script with a non-zero status to indicate an error
end

infasta = ARGS[1]
ingtf = ARGS[2]
cutoff = parse(Int, ARGS[3])
#cutoff = 300000000 #The BAI index format can handle individual chromosomes up to 512 Mbp (2^29 bases) in length
#TODO: have variable cutoffs, to generate more uniform fragment size

###### Read GTF
reader = open(ingtf, "r")

gdf = DataFrame(seqid = [], seqstart = [], seqend = [], attributes = [])
largeseqlen = Dict(); #genomeversion = ""
for line in readlines(reader)
    columns = split(line)
    #println(columns[4])
    if startswith(line, "##sequence-region") && parse(Int, columns[4]) > cutoff
        largeseqlen[columns[2]] = parse(Int, columns[4])
    # elseif startswith(line, "#!genome-version")
    #     genomeversion = columns[2]
    elseif startswith(line, "#") 
        continue
    elseif haskey(largeseqlen, columns[1]) && columns[3] == "gene"
        gid = replace(columns[10], "\"" => "", ";" => "")
        push!(gdf, (columns[1], parse(Int, columns[4]),  parse(Int, columns[5]), gid))
    end
end
close(reader)
#largeseqlen
#genomeversion

#sort genes on a chromosome according to their positions
sort!(gdf, [:seqid, :seqstart])

#process genes chromosome by chromosome
ggdfs = groupby(gdf, :seqid);
#first(ggdfs[1], 10)

###### Identify breakpoints
BP_df = DataFrame(Seqid = [], IRstart = [], IRend = [], Breakpoint = [])
for ggdf in ggdfs
    cutoffs = cutoff
    cutoffs_n = Int[]

    #sequence with 2 genes at least
    #BP located between genes
    if nrow(ggdf) > 1
        for n in 1:nrow(ggdf)-1
            while ggdf[n, :seqend] < cutoffs < ggdf[n+1, :seqstart]
                push!(BP_df, (ggdf[n, :seqid], ggdf[n, :seqend], ggdf[n+1, :seqstart], cutoffs))
                push!(cutoffs_n, cutoffs/cutoff)
                cutoffs += cutoff
            end
        end
    end

    # println( cutoffs |> typeof)
    # println( ggdf[1, :seqstart] |> typeof)
    # println( largeseqlen[ggdf[end, :seqid]] |> typeof)
    # println( cutoffs |> typeof)
    
    #BP located before the first gene
    while cutoffs < ggdf[1, :seqstart]
        push!(BP_df, (ggdf[1, :seqid], 1, ggdf[1, :seqstart], cutoffs))
        push!(cutoffs_n, cutoffs/cutoff)
        cutoffs += cutoff
    end

    #BP located after the last gene
    while ggdf[end, :seqend] < cutoffs < largeseqlen[ggdf[end, :seqid]]
        push!(BP_df, (ggdf[end, :seqid], ggdf[end, :seqend], largeseqlen[ggdf[end, :seqid]], cutoffs))
        push!(cutoffs_n, cutoffs/cutoff)
        cutoffs += cutoff
    end
    
    #if large sequence doesn't have or miss any BP
    if cutoffs_n == [] || cutoffs_n[end] != length(cutoffs_n)
        println("Warning: no breakpoint is found due to no gene, or missed due to the expected located within a gene!")
    end
end
#BP_df

###### Break large sequences based on the BP identifed above
outfasta = splitext(basename(infasta))[1]
isfile("$(outfasta).split.fa") && rm("$(outfasta).split.fa")
writer = open("$(outfasta).split.fa", "a")

reader2 = open(FASTA.Reader, infasta)

newseq_df = DataFrame(Seqid = [], New_Seqid = [], Length = [])
for record in reader2
    #split large sequences
    if seqsize(record) > cutoff
        seqdf = filter(:Seqid => x -> x == identifier(record), BP_df)
        for (n, bp) in enumerate(seqdf.Breakpoint)
            #println(identifier(record), "\t", bp)
            start = bp-cutoff+1
            mark = ('a':'z')[n]
            write(writer, ">", identifier(record), "$(mark)", "\t$(start)_$(bp)\n", sequence(record)[start:bp], "\n")
            
            #original seq id, new id, and new sequence length to a dataframe
            push!(newseq_df, (identifier(record), identifier(record) * "$mark", cutoff))

            #the last fragment that starts from the last BP to sequence end
            if n == nrow(seqdf)
                mark = ('a':'z')[n+1]
                write(writer, ">", identifier(record), "$(mark)", "\t$(bp+1)_$(seqsize(record))\n", sequence(record)[bp+1:seqsize(record)], "\n")
                push!(newseq_df, (identifier(record), identifier(record) * "$(mark)", seqsize(record) - bp))
            end
        end
    else
        #write the small seuences without any change
        write(writer, ">", identifier(record),"\n", sequence(record), "\n")
        push!(newseq_df, (identifier(record), identifier(record), seqsize(record)))
    end
end
close(reader2)
close(writer)
#newseq_df

###### Transfer annotation
outgtf = splitext(basename(ingtf))[1]
isfile("$(outgtf).split.gtf") && rm("$(outgtf).split.gtf")
writer2 = open("$(outgtf).split.gtf", "a")

f = open(ingtf, "r")
seqs_id_used = String[]
for line in readlines(f)
    if startswith(line, "##gtf-version")
        #println(line)
        write(writer2, line, "\n")
        # for key in sort(collect(keys(new_largeseqlen)))
        #     #println("##sequence-region $key 1 $(new_largeseqlen[key])")
        #     write(writer2, "##sequence-region $key 1 $(new_largeseqlen[key])", "\n")
        # end
        for row in eachrow(newseq_df)
            write(writer2, "##sequence-region $(row.New_Seqid) 1 $(row.Length)", "\n")
        end
    elseif startswith(line, "##sequence-region")
        continue
    elseif startswith(line, "#")
        #println(line)
        write(writer2, line, "\n")
    else
        #CAUTION: use "\t" as delimiter rather than "\s" in case the attributes column has "\s"! 
        columns = split(line, "\t")
        
        columns[3] == "region" && continue
        
        #BP dataframe
        seqdf = filter(:Seqid => x -> x == columns[1], BP_df)
        #new length dataframe
        seqdf_id = filter(:Seqid => x -> x == columns[1], newseq_df)
        
        for (n, bp) in enumerate(sort(seqdf.Breakpoint, rev=true))
            fstart = parse(Int, columns[4])
            fend = parse(Int, columns[5])

            # #the line with "region"
            # if columns[3] == "region"
            #     columns[9] = replace(columns[9], "$(columns[1])" => seqdf_id.New_Seqid[1])
            #     columns[1] = seqdf_id.New_Seqid[1]
            #     columns[4] = "1"
            #     columns[5] = "$(seqdf_id.Length[1])"
            #     push!(seqs_id_used, seqdf_id.New_Seqid[1])
            # end
            
            #the first fragment (the last after reverse)with the name "XXbp1", but their coordinates don't change
            if columns[3] != "region" && n == length(seqdf.Breakpoint) #&& columns[1] in seqs_id_used
                columns[1] = seqdf_id.New_Seqid[1]
            end

            # #generate a new "region" line
            # if columns[3] != "region" && n != length(seqdf.Breakpoint) && seqdf_id.New_Seqid[1] ∉ seqs_id_used
            #     write(writer2, seqdf_id.New_Seqid[1], "\t", genomeversion, "\tregion\t1\t", "end", "\t.\t.\t.\tgene_id \"region:$(columns[1])\"; Alias \"XX\"; ID \"region:$(columns[1])\";\n")
            #     push!(seqs_id_used, seqdf_id.New_Seqid[n])
            # end

            #the following fragments with the names "XXbp2, XXbp3, ...,", and coordinates have changed!
            if bp < fstart
                columns[1] = sort(seqdf_id.New_Seqid, rev=true)[n]
                columns[4] = "$(fstart - bp)"
                columns[5] = "$(fend - bp)"
                break
            end
        end

        # #copy with the line with "region"
        # if columns[1] ∉ seqs_id_used
        #     seqdf_id2 = filter(:New_Seqid => x -> x == columns[1], newseq_df)
        #     #println(columns[1], "\t", seqdf_id2.Length[1], "\t", seqs_id_used)
        #     write(writer2, columns[1], "\t", genomeversion, "\tregion\t1\t", "$(seqdf_id2.Length[1])", "\t.\t.\t.\tgene_id \"region:$(columns[1])\"; Alias \"X$(columns[1])\"; ID \"region:$(columns[1])\";\n")
        #     push!(seqs_id_used, columns[1])
        # end
        
        #write all lines including those without BP
        #println(join(columns, "\t"))
        write(writer2, join(columns, "\t"), "\n")
    end
end
close(f)
close(writer2)