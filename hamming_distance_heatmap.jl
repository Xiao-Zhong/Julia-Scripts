using BioAlignments
using Plots
plotlyjs()

#input setting
project = "GWA_index1-check"
remove_collided_indexes = true
#seqs = ("AAGACTAC", "TGCCGTTA", "TTGGATCT", "TCCTCCAA", "CGAGTCGA")
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

function dist_plot(scores)
    heatmap(scores,
        colorbar = false,
        c = palette([:white, :green], 10),
        xticks = (1:nrow, seqs),
        xrot = 90,
        yticks = (1:ncol, seqs),
        size = (700, 500),
        linewidth = 5,
        #title = project * " barcode hamming distance"
    )
    fontsize = 10
    annotate!([(i, j, text(scores[i,j], fontsize, :black)) for i in 1:nrow for j in 1:ncol])
    #savefig("$project-Collision.html")
end

#"scores" is "mismatches"
scores = [score(pairalign(HammingDistance(), i, j)) for i in seqs, j in seqs]
nrow, ncol = size(scores)
dist_plot(scores)

#if remove_collided_indexes == true
    bad = [pair[1] for pair in findall(x-> x == 1, scores)] |> unique
    scores_bad = scores[bad, bad]
    nrow, ncol = size(scores_bad)
    dist_plot(scores_bad)

    good = [pair[1] for pair in findall(x-> x > 1, scores)] |> unique
    scores_good = scores[good, good]
    nrow, ncol = size(scores_good)
    dist_plot(scores_good)
#end
