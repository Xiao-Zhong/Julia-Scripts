# seqs = """
# AGCTGAAG
# ATTCCGTG
# TATGCCGC
# TCAGCTCA
# ATTAGGAG
# CAGCAATA
# GCCAAGCT
# TCCGTTAA
# GTGCAACG
# AGTAACGC
# CATAGCCA
# CACTAGTA
# TTAGTGCG
# TCGATACA
# ATAGTGAC
# CAAGGTGA
# TAGACCAA
# CGGTAGAG
# TCAGCATC
# AGAAGCAA
# GCAGGTTC
# AAGTGTCT
# CTACCGAA
# TAGAGCTC
# ATGTCAAG
# GACTTGAC
# CTACAATG
# TCTCAGCA
# AGACACAC
# AATACGCG
# GCACACAT
# GCACCTAA
# TGCTGCTC
# TGGCACCA
# AGATGGAT
# GAATTGTG
# GAGCACTG
# GTTGCGGA
# AATGGAAC
# TCAGAGGT
# GCAACAAT
# GTCGATCG
# ATGGTAGC
# CGCCAATT
# GACAATTG
# ATATTCCG
# TCTACCTC
# TCGTCGTG
# ATGAGAAC
# AATGACCA
# CAGACGCT
# TCGAACTG
# CGCTTCCA
# TATTCCTG
# CAAGTTAC
# CACATAAA
# ACCAATCG
# AGATAGTG
# AGAGGTTA
# CTGACCGT
# GCATGGAG
# GCGTCACT
# TCACCACG
# AGACCTGA
# GCCGATAT
# CGATACCT
# CTCGACAT
# GAGATCGC
# CGGTCTCT
# TAACTCAC
# CACAATGA
# GACTGACG
# CTTAAGAC
# GAGTGTAG
# TGCACATC
# CGATGTCG
# AACACCGA
# GATCCATG
# AGCTACAT
# CGCTGTAA
# CACTACCG
# TGGCTTAG
# TCCAGACG
# AGTGGCAT
# TGTACCGA
# AAGACTAC
# ATCATTCC
# ACCACTGT
# ACTATGCA
# AACTCACC
# CTGTAGCC
# GCTCGGTA
# ACACGACC
# GGAGAACA
# CATCAAGT"""

using BioAlignments
using Plots
plotlyjs()

#input setting
project = "GWA-index-check"
indexes_split = true
split_cutoff = 1
seqs = """
AACTGCAA
AACGCATT
AGTCGCGA
CAGGTCTG
CTTGCATA
ATCCTCTT
GTCCTATA
CTGCCTTA
GCTCACGA
CTGGCATA
ACCTCCAA
GCTAACGA
CAGATCTG
ATCCTGTA
AGTCACTA
AACGCTTA
GCGAGTAA
GCGATTAC
GCCACATA
GCATCATA
CGGATTGC
CTTATTGC"""
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

scores = distance(seqs)
dist_plot(scores, seqs)

if indexes_split == true
    bad = [pair[1] for pair in findall(x-> x == split_cutoff, scores)] |> unique
    dist_plot(distance(seqs[bad]), seqs[bad])

    seqs_good = setdiff(seqs, seqs[bad])
    dist_plot(distance(seqs_good), seqs_good)
end
