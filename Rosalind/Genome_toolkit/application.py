from genome_toolkit import genomeToolkit

gt = genomeToolkit()

seq = 'AAAGAAAATTGA'
kmer = 'AA'

print(gt.count_kmer(seq, kmer))