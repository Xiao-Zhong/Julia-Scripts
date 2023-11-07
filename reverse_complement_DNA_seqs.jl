using BioSequences

seqs = """
GAACTGAGCG
CCGTATGTTC
AAAAAAAAAA
AAAAATTTTT
"""

for seq in split(seqs)
    println(reverse_complement(LongDNA{2}(seq)))
end
