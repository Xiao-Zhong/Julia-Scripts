using BioAlignments
using Plots
plotlyjs()

#input
project = "gwa"
#seqs = ("AAGACTAC", "TGCCGTTA", "TTGGATCT", "TCCTCCAA", "CGAGTCGA")
seqs = """
CGCTTCCA	CAGTTCGA
TATTCCTG	TGGAAGCG
CAAGTTAC	CAGGAATA
CACATAAACA	GTTTGAAGTA
ACCAATCGTG	CCCGTGGAGA
"""
seqs = replace.(split(seqs, "\n"), "\t" => "")

#"scores" is "mismatches"
scores = [score(pairalign(HammingDistance(), i, j)) for i in seqs, j in seqs]
nrow, ncol = size(scores)

heatmap(scores,
    colorbar = false,
    c = palette([:white, :green], 10),
    xticks = (1:nrow, seqs),
    xrot = 90,
    yticks = (1:ncol, seqs),
    size = (2000, 2000),
    linewidth = 5,
    #title = project * " barcode hamming distance"
)

fontsize = 10
annotate!([(i, j, text(scores[i,j], fontsize, :black)) for i in 1:nrow for j in 1:ncol])
savefig("$project.html")
