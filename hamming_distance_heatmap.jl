using BioAlignments
using Plots
#plotlyjs()

#input setting
seqs = """
TCATCCTT
AACACTCT
CACCTAGA
AGTTCATG
GTTGGTGT
GCTACGCA
TCAACTGC
AAGCGAAT
GTGTTACA
CAAGCCAT
CTCTCGTG
TCGACAAC
TCGATGTT
CAAGGAAG
ATTGATGC
TCGCAGAT
GCAGAGAC
CTGCGAGA
CAACCAAC
ATCATGCG
TCTGAGTC
TCGCCTGT
GCGCAATT
AGACGCCT
CAGGCAGA
TCCGCGAT
CTCGTACG
CACACATA
CGTCAAGA
TTCGCGCA
CGACTACG
GAAGGTAT
TTGGCATG
CGAATTCA
TTAGTTGC
GATGCCAA
AGTTGCCG
GTCCACCT
ATCAAGGT
GAACCAGA
CATGTTCT
TCACTGTG
ATTGAGCT
GATAGAGA
TCTAGAGC
GAATCGCA
CTTCACGT
CTCCGGTT
TGTGACTA
GCTTCCAG
CATCCTGT
GTAATACG
GCCAACAA
CATGACAC
TGCAATGC
CACATTCG
CAATCCGA
CATCGACG
GTGCGCTT
ATAGCGTT
GAGTAAGA
CTGACACA
ATACGTGT
GACCGAGT
GCAGTTAG
CGTTCGTC
CGTTAACG
TCGAGCAT
GCCGTAAC
GAGCTGTA
AGGAAGAT
CTAACAAG
CGCTCAGA
TAACGACA
ATGCCTAA
AACGTGAT
GACTAGTA
ATTGGCTC
GATGAATC
AGCAGGAA"""
project = "GWA-index-check"

seqs = replace.(split(seqs, "\n"), "\t" => "")

function dist_plot(scores, seqs)
    nrow, ncol = size(scores)
    p = heatmap(scores,
        colorbar = false,
        c = palette([:white, :green], 10),
        xticks = (1:nrow, seqs),
        xrot = 90,
        yticks = (1:ncol, seqs),
        size = (1800, 1800),
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
savefig("index_heatmap_v2.pdf")

indexes_split = true
split_cutoff = 3
#if indexes_split == true
    bad = [pair[1] for pair in findall(x-> 0 < x < split_cutoff, distance(seqs))] |> unique |> sort
    dist_plot(distance(seqs[bad]), seqs[bad])
    savefig("index_heatmap_good.pdf")

    seqs_good = setdiff(seqs, seqs[bad])
    dist_plot(distance(seqs_good), seqs_good)
    savefig("index_heatmap_bad.pdf")
