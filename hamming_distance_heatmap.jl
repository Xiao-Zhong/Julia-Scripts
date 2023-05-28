using BioAlignments
using Plots
plotlyjs()

#input setting
seqs = """
AACTGCAA
AACGCATT
AGTCGCGA
AGTCACTA
AACGCTTA
CTGTAGCC
GCTCGGTA
ACACGACC
GGAGAACA
CATCAAGT"""
project = "GWA-index-check"
indexes_split = true
split_cutoff = 3

seqs = replace.(split(seqs, "\n"), "\t" => "")

function dist_plot(scores, seqs)
    nrow, ncol = size(scores)
    p = heatmap(scores,
        colorbar = false,
        c = palette([:white, :green], 10),
        xticks = (1:nrow, seqs),
        xrot = 90,
        yticks = (1:ncol, seqs),
        size = (400, 400),
        linewidth = 5,
        #title = project * " barcode hamming distance"
    )
    fontsize = 10
    annotate!([(i, j, text(scores[i,j], fontsize, :black)) for i in 1:nrow for j in 1:ncol])
    return p
    #savefig("test.pdf")
end

#"scores" is "mismatches"
distance(seqs) = [score(pairalign(HammingDistance(), i, j)) for i in seqs, j in seqs] 
dist_plot(distance(seqs), seqs)
savefig("test1.html")

if indexes_split == true
    bad = [pair[1] for pair in findall(x-> 0 < x < split_cutoff, distance(seqs))] |> unique |> sort
    dist_plot(distance(seqs[bad]), seqs[bad])
    savefig("test2.pdf")

    seqs_good = setdiff(seqs, seqs[bad])
    dist_plot(distance(seqs_good), seqs_good)
    savefig("test3.pdf")
end
