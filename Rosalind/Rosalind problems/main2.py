from bio_seq import bio_seq
from utilities import read_FASTA, readTextFile, writeTextFile
from time import perf_counter

start_time = perf_counter()

test_dna = bio_seq()
test_dna.generate_rnd_seq(40, 'RNA')

print(test_dna.get_seq_info())
print(test_dna.get_seq_biotype())
print(test_dna.gc_content())
print(test_dna.gc_content_subsec())
print(test_dna.translate_seq())

for rf in test_dna.gen_reading_frames():
    print(rf)

print(test_dna.all_proteins_from_orfs())

writeTextFile('test.txt', test_dna.seq)
for rf in test_dna.gen_reading_frames():
    writeTextFile('test.txt', str(rf), 'a')

fasta = read_FASTA('Rosalind problems/fasta_sample.txt')
print(fasta)

elapsed_time = perf_counter()
execution_time = elapsed_time - start_time