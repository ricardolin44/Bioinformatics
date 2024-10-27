class genomeToolkit:
    def __init__(self):
        print('Genome Toolkit initiated')

    def count_kmer(self, sequence, kmer):
        kmer_count = 0
        for position in range(len(sequence)-(len(kmer)-1)):
            if sequence[position:position+len(kmer)] == kmer:
                kmer_count += 1
        return kmer_count