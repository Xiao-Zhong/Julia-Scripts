using Glob
using XAM
using DataFrames
using CSV
#using BenchmarkTools

function extract_names(files::Vector{String})
    df = DataFrame(Sample = String[], Reads = Vector{String}[])
    for file in files
        sample = split(basename(file), ".")[1]
        sample = split(sample, "_")[1] # Note: check if the same sample names can be extracted.
        println(sample)
        
        read_name = String[]
        reader = open(BAM.Reader, file)
        # for record in reader
        #     push!(read_name, replace(BAM.tempname(record), r"_\w+$" => ""))
        # end
        record = BAM.Record()
        while !eof(reader)
            empty!(record)
            read!(reader, record)
            push!(read_name, replace(BAM.tempname(record), r"_\w+$" => ""))
        end
        close(reader)

        push!(df, [sample, read_name])
    end
    return df
end

function overlapping(umi::Vector{String}, agent::Vector{String}, p::String)
    df_merge = leftjoin(extract_names(umi), extract_names(agent), on = :Sample, makeunique=true)
    df_merge.Overlapping = intersect.(df_merge.Reads, df_merge.Reads_1)

    df_merge.Umi = df_merge.Reads .|> length
    df_merge.UmiUnique = df_merge.Reads .|> unique .|> length
    df_merge.Agent = df_merge.Reads_1 .|> length
    df_merge.AgentUnique = df_merge.Reads_1 .|> unique .|> length
    df_merge.OverlappingSize = df_merge.Overlapping .|> length
    #println(df_merge)
    CSV.write("$(p)_overlappping_reads_df.csv", df_merge[!, [1, 5, 6, 7, 8, 9]])
end

##### set project, bam file path #####

project = "S94-RD-50"
umi_files = glob("./UMI-tools/aligned-deduplicated/*.bam")
agent_files = glob("./aligned-deduplicated/*.bam")

# Usage: julia -t 46 overlapping_reads_df.jl
# Run it in parallel to save time! 
#####

println(umi_files, "\n", agent_files)
# extract_names(umi_files) |> println
# extract_names(agent_files) |> println
overlapping(umi_files, agent_files, project)
