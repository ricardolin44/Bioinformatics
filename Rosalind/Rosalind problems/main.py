from dna_toolkit import *
import random

randDNAStr = ''.join([random.choice(nucleotides) for nuc in range(333)])
# DNAStr = validseq(randDNAStr)
DNAStr = 'AUGGACACUACCGUCCGACAGCGCUGUUCAAUCUCGAUCGAAGCGUCUUGUAGCGUGCUU'.replace('U','T')

print(f'\nSequence: {DNAStr}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide Frequency: {countseq(DNAStr)}\n')
print(f'[3] + DNA/RNA Transcription: {transcription(DNAStr)}\n')
print(f"[4] + DNA String + Complement + Reverse Complement: \n5' {DNAStr} 3'")
print(f"   {''.join(['|'for c in range(len(DNAStr))])}")
print(f"3' {revtrans1(DNAStr)[::-1]} 5' [Complement]")
print(f"5' {revtrans1(DNAStr)} 3' [Rev. Complement]")
print(f'[5] + GC Content: {gccontent(DNAStr)}%\n')
print((''.join(translation(DNAStr))).strip('_'))
print(codon_usage(DNAStr, 'A'))