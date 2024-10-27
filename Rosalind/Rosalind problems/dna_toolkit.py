import collections
nucleotides = ['A','T','G','C']
DNA_Complement = {'A':'T', 'T':'A','C':'G','G':'C'}

DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}


#Validating and counting nucleotides
def validseq(seq):
    tmp_seq = seq.upper()
    for nuc in tmp_seq:
        if nuc not in nucleotides:
            return False
    return tmp_seq

def countseq(seq):
    cntseq = {'A':0,'C':0,'G':0,'T':0}
    for nuc in seq:
        cntseq[nuc] += 1
    return cntseq
    '''return collections.Counter(seq)'''
    '''print(' '.join(str(val) for key,val in result.items()))'''

#Transcribing DNA into RNA
def transcription_old(seq):
    rna_seq = []
    dna_seq = seq.upper()
    for nuc in dna_seq:
        if nuc == 'T':
            nuc = 'U'
            rna_seq.append(nuc)
        else:
            rna_seq.append(nuc)
    return ''.join(rna_seq)

def transcription(seq):
    dna_seq = seq.upper()
    return dna_seq.replace('T','U')

#Reverse Transcription
def revtrans1(seq):
    return ''.join(DNA_Complement[nuc] for nuc in seq)[::-1]

def revtrans2(seq):
    mapping = str.maketrans('ATGC','TACG') 
    return seq.translate(mapping)[::-1]

#GC Contents
def gccontent(seq):
    return round((seq.count('C')+seq.count('G'))/len(seq)*100)

#Hamming Distances
def h_d_loop(str1, str2):
    h_distance = 0
    for itr in range(len(str1)):
        if str1[itr] != str2[itr]:
            h_distance += 1
        return h_distance

def h_d_set(str1, str2):
    nucleotide_set_1 = set([(x,y) for x,y in enumerate(str1)])
    nucleotide_set_2 = set([(x,y) for x,y in enumerate(str2)])
    return len(nucleotide_set_1.difference(nucleotide_set_2))

def h_d_zip(str1, str2):
    zipped_dna = zip(str1, str2)
    h_distance = [(nuc1,nuc2)for nuc1, nuc2 in zipped_dna if nuc1!=nuc2]
    return len(h_distance)

def codon_usage(seq, aminoacid):
    codon = {}
    for i in range(0, len(seq)-2, 3):
        if DNA_Codons[seq[i:i+3]] == aminoacid:
            if seq[i:i+3] not in codon:
                codon[seq[i:i+3]] = 1
            else: codon[seq[i:i+3]] += 1

    freq = sum(codon.values())
    codonfreq = {key: value/freq for key, value in codon.items()}
    return codonfreq

def translation(seq, init_pos=0):
    # seq.replace('U','T')
    return [DNA_Codons[seq[pos:pos+3]] for pos in range(init_pos, len(seq)-2, 3)]

def gen_reading_frames(seq):
    frames = []
    frames.append(translation(seq, 0))
    frames.append(translation(seq, 1))
    frames.append(translation(seq, 2))
    frames.append(translation(revtrans1(seq),0))
    frames.append(translation(revtrans1(seq),1))
    frames.append(translation(revtrans1(seq),2))

def proteins_from_rf(aa_seq):
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == '_':
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == 'M':
                current_prot.append('')
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    if endReadPos>startReadPos:
        rfs = gen_reading_frames(seq[startReadPos:endReadPos])
    else:
        rfs = gen_reading_frames(seq)
    res = []
    for rf in rfs:
        proteins = proteins_from_rf(rf)
        for p in proteins:
            res.append(p)

    if ordered:
        return sorted(res, key=len, reverse=True)
    return res