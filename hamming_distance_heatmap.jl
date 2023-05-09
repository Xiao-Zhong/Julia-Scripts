using BioAlignments
using Plots
plotlyjs()

#input
project = "YY"
#seqs = ("AAGACTAC", "TGCCGTTA", "TTGGATCT", "TCCTCCAA", "CGAGTCGA")
seqs = """
AAGACTAC
TGCCGTTA
TTGGATCT
TCCTCCAA
TCCTCCAA
""" |> split

#"scores" is "mismatches"
scores = [score(pairalign(HammingDistance(), i, j)) for i in seqs, j in seqs]
nrow, ncol = size(scores)

heatmap(scores,
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
savefig("$project.html")
