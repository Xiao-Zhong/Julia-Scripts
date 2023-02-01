#samtools view -F 4  PCB-09-PDX_S1_L001.Aligned.sortedByCoord.out.bam |cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' |sort |uniq -c|awk '{print $2"\t"$1}' > mapped_read_length_count.txt
using XAM
using StatsBase

# Open a BAM file.
reader = open(BAM.Reader, "PCB-09-PDX_S1_L001.Aligned.sortedByCoord.out.bam")

l = []
# Iterate over BAM records.
for record in reader
    # `record` is a BAM.Record object.
    if BAM.ismapped(record)
        #println(BAM.sequence(record), "\t", BAM.seqlength(record))
        push!(l, BAM.seqlength(record))
    end
end
# Close the BAM file.
close(reader)

d = countmap(l) |> sort
