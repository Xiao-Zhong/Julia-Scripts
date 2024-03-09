#import Pkg; Pkg.add(["BioAlignments", "PlotlyJS", "PlotlyJS"]) #install the pacakges
using BioAlignments
using Plots
plotlyjs()

#input settings
indexes = """
O.T.2.Amp	B02	TAGAGCTC			GWA-KQ-971
O.T.3.Amp	C02	ATGTCAAG			GWA-KQ-971
O.T.4.Amp	D02	GCATCATA			GWA-KQ-971
O.T.5.Amp	E02	GACTTGAC			GWA-KQ-971
O.T.6.Amp	F02	CTACAATG			GWA-KQ-971
O.T.7.Amp	G02	TCTCAGCA			GWA-KQ-971
""" #paste indexes here
project = "index_check.html" #change file name, and extension (e.g., png, pdf and svg, etc.)
plot_size = (700, 600) #adjust when required, very likely when there're too many indexes
font_size = 10 #adjust when required, ...

function dist_plot(scores, seqs, plot_s, font_s)
    nrow, ncol = size(scores)
    p = heatmap(scores,
        colorbar = false,
        c = palette([:white, :green], 10),
        xticks = (1:nrow, seqs),
        xrot = 90,
        yticks = (1:ncol, seqs),
        size = plot_s, #adjust when required, very likely when there're too many indexes
        linewidth = 5,
        #title = project * " barcode hamming distance"
    )
    fontsize = font_s #adjust when required, ...
    annotate!([(i, j, text(scores[i,j], fontsize, :black)) for i in 1:nrow for j in 1:ncol])
    return p
end

function read_index(idxs)
    sample_dict = Dict()
    for line in split(idxs, "\n")[1:end-1]
        cols = split(line, "\t")[[1, 3, 5]] #extract sample, index1, index2
        sample_dict[cols[2] * cols[3]] = cols[1] * "-" * cols[2] * cols[3]
    end
    return sample_dict |> keys |> collect, sample_dict |> values |> collect
end

#one line function for hamming distance
distance(x) = [score(pairalign(HammingDistance(), i, j)) for i in x, j in x] 

(seqs, seqs_samples) = read_index(indexes)
dist_plot(distance(seqs), seqs_samples, plot_size, font_size)
savefig(project) 

# indexes_split = true
# split_cutoff = 3
# #if indexes_split == true
#     bad = [pair[1] for pair in findall(x-> 0 < x < split_cutoff, distance(seqs))] |> unique |> sort
#     dist_plot(distance(seqs[bad]), seqs[bad])
#     savefig("index_heatmap_good.pdf")

#     seqs_good = setdiff(seqs, seqs[bad])
#     dist_plot(distance(seqs_good), seqs_good)
#     savefig("index_heatmap_bad.pdf")
