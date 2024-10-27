from collections import Counter
import math
import requests
import re
import itertools

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
protein_mRNA = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATA', 'ATT', 'ATC'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'], 'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'], 'Y': ['TAT', 'TAC'],
    '_': ['TAA', 'TAG', 'TGA']
    }

amino_acid_complement = {
    'A' : 'T',
    'T' : 'A',
    'G' : 'C',
    'C' : 'G'
}
# for key,value in DNA_Codons.items():
#     if value in protein_mRNA:
#         protein_mRNA[value].append(key)
#     else:
#         protein_mRNA[value] = [key]
# print(protein_mRNA)

seq = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
dic_seq = {'A':0,'C':0,'G':0,'T':0}
def count_nuc(seq):
    for nuc in seq:
        dic_seq[nuc] += 1
    return dic_seq

def countseq(seq):
    cntseq = {'A':0,'C':0,'G':0,'T':0}
    for nuc in seq:
        cntseq[nuc] += 1
    return cntseq

def translation_all(seq, init_pos=0):
    res = [DNA_Codons[seq[x:x+3]] for x in range(init_pos, len(seq)-2, 3)]
    return ''.join(res[:-1])

def complement(seq):
    return ''.join([amino_acid_complement[amino_acid] for amino_acid in seq[::-1]])

# ans = list(count_nuc(seq).values())
# print((an for an in ans), sep=' ')
# print(countseq(seq))
# print(' '.join(str(val) for keyss,val in countseq(seq).items()))

# print(countseq(seq).items())
# print(' '.join(str(val) for keyss,val in countseq(seq).items()))

text1 = '''>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT
>Rosalind_8
ATGGCACT
'''

def readFastaintoDict(text):
    text = '\n'.join(line.strip() for line in text.splitlines() if line.strip())
    listA = text.split('\n')
    listdna = ['' for index in listA if index[0]=='>']
    indexnum = -1
    for index in listA:
        if index[0] != '>':
            listdna[indexnum] += index
        else:
            indexnum += 1
    listname = [index[1:] for index in listA if index[0]=='>']
    length = len(min(listdna, key=len))
    res = {listname[num]: listdna[num].upper() for num in range(len(listdna))}
    return [res, length]

def motif(text, subtext):
    positions = []
    position = text.find(subtext)
    while position != -1:
        positions.append(position+1)     #plus 1 to make it easy to count
        position = text.find(subtext, position+1)   #plus 1 to find the subtext after the position
    return positions

def consensus(text):
    dictdna, length = readFastaintoDict(text)[0], readFastaintoDict(text)[1]
    counting = {'A':[0] * length, 'C':[0] * length, 'G':[0] * length, 'T':[0] * length}
    for value in dictdna.values():
        for num, nuc in enumerate(value):
            counting[nuc][num] += 1
    cons =  [max(counting, key=lambda keys: counting[keys][index]) for index in range(length)]  #all keys to be compared
    matrix = [counting['A'],counting['C'],counting['G'],counting['T']]
    name = ['A', 'C', 'G', 'T']
    print(''.join(cons))
    for num, row in enumerate(matrix):
        res = f'''{name[num]}: {' '.join(map(str,row))}'''
        print(res)
    return [counting, cons]

text2 = 'CCACTGTTGTAGCGCTGTTGTCTGTTGTTCCTGTTGTGTCCTGTTGTCCTGTTGTGGCCTCTGTTGTCGCATTCTGTTGTATACTGTTGTTGATTCAACTGTTGTGACTGTTGTCTGTTGTACTGTTGTTCGGCAACTACTGTTGTGCTGTTGTCCTCCTGTTGTGCTGTTGTCTGTTGTGGCTGTTGTTCTGTTGTCTGTTGTACGTCTGTTGTACTGTTGTCATCTGTTGTTAGCTGTTGTGGAGCGCTGTTGTCGGGCGAGACTGTTGTCTGTTGTCTGTTGTGTTCGCTGTTGTGCCTGTTGTGCTGTTGTACATCTGTTGTCTGTTGTTCTCTGTTGTTCTGTTGTGAGACTGTTGTGGATTCACTGTTGTGTAGTACTGTTGTTTAGGATAGACGCTTAGCTCTGTTGTAACCCTGTTGTCTGTTGTCGTACCTGTTGTCCCTGTTGTACTGTTGTGAAGATGCTGTTGTCCACCTGTTGTAGGCTGTTGTCTGTTGTTCCTGTTGTACTGTTGTCTGTTGTTGCCTGTTGTAATGTGCACCTCCTGTTGTAGAGACTGTTGTCTGTTGTCTGTTGTCCTGTTGTGGGCCTGTTGTCTGTTGTCTGTTGTATAACTGTTGTTGCTGTTGTTTGGTTGCTGTTGTGCCCGTAACCTGTTGTCATGACTGTTGTTCTGTTGTCTGTTGTTACGATCACTACCACGCTGTTGTACTGTTGTCTGTTGTAATGACTCACACCTGTTGTCCTGTTGTTGATGCTGTTGTCTGTTGTCTGTTGTGTTGCACTGTTGTCTGTTGTCTGTTGTTCCTACAACCTCTGTTGTCACTGTTGTAGTCTGTTGTCCTGTTGTCTGTTGTCTGTTGTAAAA'
text3 = 'CTGTTGTCT'

consensus(text1)
print(' '.join(map(str,motif(text2,text3))))

# students = {'Alice': 85, 'Bob': 92, 'Charlie': 78, 'David': 95}
# highest_scorer = max(students, key=students.get)
# print(highest_scorer)  # Output: David

# counting = {'A': [6, 1, 0, 0, 5, 6, 0, 0],
#             'C': [0, 0, 1, 4, 3, 0, 7, 1],
#             'G': [1, 1, 7, 4, 0, 1, 0, 0],
#             'T': [1, 6, 0, 0, 0, 1, 1, 7]}

# max_values_per_position = [max(values) for values in zip(*counting.values())]

# print("Max values for each position:")
# print(max_values_per_position)

def mendelsprob(k,m,n):
    sample = ['AA']*k + ['Aa']*m + ['aa']*n
    combs = []
    for num, index1 in enumerate(sample):
        temp = sample[num+1:]
        for index2 in temp:
            combs += [a+b for a in index1 for b in index2]
    prob = (1 - combs.count('aa')/len(combs))
    return prob

print(mendelsprob(30,15,20))

def overlap(text, k):
    dictdna, length = readFastaintoDict(text)[0], readFastaintoDict(text)[1]
    listname, listdna = list(dictdna.keys()), list(dictdna.values())
    res = []
    for num1, index1 in enumerate(listdna):
        temp = listdna[:num1] + listdna[num1+1:]
        for num2, index2 in enumerate(temp):
            if index1[-k:]==index2[:k]:
                if num2>=num1:
                    num3 = num2+1
                else:
                    num3 = num2
                combination = f'{listname[num1]} {listname[num3]}'
                res.append(combination)
    return res
text4 = '''>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG
'''
print('\n'.join(overlap(text4, 3)))

text5 = '''>Rosalind_1
GATTACA
>Rosalind_2
TAGACCA
>Rosalind_3
ATACA
'''

def sharedmotif(text):
    dictdna, length = readFastaintoDict(text)[0], readFastaintoDict(text)[1]
    listname, listdna = list(dictdna.keys()), list(dictdna.values())
    shortest = min(listdna, key=len)
    temp = listdna
    temp.remove(shortest)
    res = []
    for num in range(length,1,-1):
        for times in range(0,length-num+1,1):
            motif = shortest[times:times+num]
            check = all(motif in index2 for index2 in temp)
            if check:
                res.append(motif)
        if check:
            break
    return res[0]   #only output the first one

# text_list = ["apple", "apricot", "application", "apartment"]
#Find a common substring
# common_substring = text_list[0]
# for text in text_list[1:]:
#     while text.find(common_substring) == -1:
#         common_substring = common_substring[:-1]

# print("Common substring:", common_substring)

print(sharedmotif(text5))

def independent(generation, N):
    #probability for Aa and Aa parent to breed Aa child is 0.5
    #probability is independant for Aa and Bb so AaBb child has the probability of 0.5*0.5 = 0.25
    #this problem is focusing on 2 pairs of alleles(4independant allele which has 16combinations, half of them has the probability 0.5 as stated above)
    #but we can tinkle the problem by changing the alleles into 3pairs(6**2=36combinations)
    #and also change the focus, not Aa but AA or another pair
    prob = 1
    samplenum = 2**generation
    for num in range(N):
        prob -= math.comb(samplenum, num)*(0.25)**num*(0.75)**(samplenum-num)
    return round(prob,3)

print(independent(2,1))

def fetch_fasta(uniprot_id):
    try:
        url = f'http://www.uniprot.org/uniprot/{uniprot_id[:6]}.fasta'
        response = requests.get(url)
        response.raise_for_status()  # Check for any errors in the HTTP response
        fasta_content = response.text
        # print("FASTA file fetched successfully.")
        return fasta_content
    except requests.exceptions.RequestException as e:
        print(f"Error fetching FASTA file: {e}")

def fetch_fasta_into_dict(text):
    protein_list = [line for line in text.splitlines()]
    all_protein = ''
    for protein in protein_list:
        info = fetch_fasta(protein)
        all_protein += info
    protein_seq = list(readFastaintoDict(all_protein)[0].values())
    protein_dict = dict(zip(protein_list,protein_seq))
    return protein_dict

def find_motif(text):

    return None

text7 = '''P07725_CD8A_RAT
Q90304_C166_CARAU
Q13VE3
P02760_HC_HUMAN
P07306_LECH_HUMAN
A1TJ10
A5WBR3
Q8WW18
P06870_KLK1_HUMAN
P01878_ALC_MOUSE
P01047_KNL2_BOVIN
Q3B391
'''

def find_motif(text):   #N-glycosylation motif
    pattern = re.compile(r'(?=(N[^P][ST][^P]))')   #(?=) returns overlapping span, but the matches will be 0 length
    protein_dict = fetch_fasta_into_dict(text)   #you can also use overlapped=True but must be used with findall (new Python regex)
    protein_name = list(protein_dict.keys())
    res = []
    for num,protein in enumerate(protein_dict.values()):
        matches = pattern.finditer(protein)
        temp_range = [match.span(0) for match in matches]
        temp_init = [index[0]+1 for index in temp_range]
        if temp_range:
            res.append(protein_name[num])
            res.append(' '.join(map(str,temp_init)))
    return '\n'.join(map(str,res))

# print(find_motif(text7))

def reverse_trans_rna_possibilities(protein_seq):
    protein_seq += '_'
    combs = 1
    for amino_acid in protein_seq:
        combs = combs*len(protein_mRNA[amino_acid])
    return combs%1000000

text8 = 'MTMNGVSINKIKAFFQKGMYMPIPHAPEMLVFTIDPIHQAS'
print(reverse_trans_rna_possibilities(text8))

def ORF(text):
    dna_dict = readFastaintoDict(text)[0]
    dna_name, dna_list = list(dna_dict.keys()), list(dna_dict.values())
    # pattern = re.compile(r'(?=(ATG(\w*)(TAG|TGA|TAA)))')
    dna_read = dna_list
    for dna in dna_list:
        dna_read = dna_read + [complement(dna)]
    res = []
    for dna in dna_read:
        # matches = pattern.finditer(dna_origin)
        # dna_seq = [match.group(1) for match in matches]
        starts = []
        for start_position in range(3):
            starts += [num for num in range(start_position, len(dna)-start_position-2, 3) if dna[num:num+3] == 'ATG']
        orfs = []
        for start in starts:
            temp1 = []
            for num in range(start, len(dna)-2, 3):
                temp1 += [dna[num:num+3]]
                if dna[num:num+3] == 'TGA' or dna[num:num+3] == 'TAG' or dna[num:num+3] == 'TAA':
                    break
            if 'TGA' in temp1 or 'TAG' in temp1 or 'TAA' in temp1:
                orfs.append(''.join(temp1))
        for orf in orfs:
            temp = translation_all(orf)
            if temp not in res:
                res.append(temp)
    return '\n'.join(map(str,res))

text9 = '''>Rosalind_1009
CTTAATTTCGATTGCAGC
'''

print(ORF(text9))

# list1 = ['AAAAA']
# list2 = list1
# list1 = list1 + ['AAAAAAA']
# print(list2)

# list3 = ['AAAAA']
# list4 = list3
# list3 += ['AAAAAAA']
# print(list4)

def palindrome(text):
    dna_dict = readFastaintoDict(text)[0]
    dna_name, dna_list = list(dna_dict.keys()), list(dna_dict.values())
    res = []
    for seq in dna_list:
        for i in range(4,13):
            seq2 = complement(seq)[::-1]
            for num1 in range(0,len(seq)-i+1,1):
                if seq[num1:num1+i] == seq2[num1:num1+i][::-1]:
                    res += [f'{num1+1} {i}']
    return '\n'.join(res)
#start 4-12,

text10 = '''>Rosalind_5375
TATAGATCGGATAGTCGTATATTGT'''
print(palindrome(text10))

def slicing(text):
    dna_dict = readFastaintoDict(text)[0]
    dna_name, dna_list = list(dna_dict.keys()), list(dna_dict.values())
    res = dna_list[0]
    for dna in dna_list[1::]:
        res = res.replace(dna,'')
    new_dna = f'>{dna_name[0]}\n{res}'
    res_list = ORF(new_dna).split('\n')
    return max(res_list, key=len)

text11 = '''>Rosalind_9452
ATGAGGACAGTGCAGGGTCAGGGAAAGTGCAGCATCGGGAGGTGGTTAGCACACAAATCG
GGACGTTGTCTATGAAGCGGCAACACAGGGATGCTTCGCAACGGTACTGGAGAGCTAAAA
CTTATTAGAGTTCTTGATATCGTTGAAGGACCACTAAAGGGGACGGATAGCCTTGAGACA
TGTTGTTACCAAATGCAACTGACGTCTTCACGCGTTCCATACCTAACCGTTTGTCCAACC
AACAGAAATCTAAACGAGTGGATTATTACAGATGGAAGCCACCTCCCCAGTGAAAGTTAC
CGCCTCG
>Rosalind_6301
GGACTTCACGCCTGGTACTTG
>Rosalind_6458
'''
print(slicing(text11))

def generate_combination(perm):
    print(math.factorial(perm))
    data = [num for num in range(1, perm+1)]
    for combination in itertools.permutations(data,len(data)):  #len(data) for how many lengths of combination 
        yield ' '.join(map(str,combination))

for combo in generate_combination(3):
    print(combo)
# print(generate_combination(3))   generate object, need to iterate

def probs_using_gc(text):
    dna = text.splitlines()[0]
    GC_content = countseq(dna)['C'] + countseq(dna)['G']
    probs_list = text.splitlines()[1].split(' ')
    res = []
    for prob in probs_list:
        temp = (float(prob)/2)**GC_content * ((1-float(prob))/2)**(len(dna)-GC_content)
        res.append(round(math.log10(temp),3))
    return ' '.join(map(str,res))

text12 = '''ACGATACAA
0.129 0.287 0.423 0.476 0.641 0.742 0.783
'''

print(probs_using_gc(text12))

def motif_search(text):
    dna_dict = readFastaintoDict(text)[0]
    dna_name, dna_list = list(dna_dict.keys()), list(dna_dict.values())
    print(dna_list[0])
    res = []
    for dna in dna_list:
        temp = ['0']
        num_temp = ['0']
        for amino_acid in dna[1:]:
            for index, no in enumerate(num_temp):
                if amino_acid == dna[int(no)]:
                    num_temp[index] = str(int(num_temp[index]) +1)
                    max_num = max(int(num1) for num1 in num_temp)
                    if no!='0' and amino_acid == dna[0] and len(num_temp) < max_num:  #finding another amino acid that has the same motif while finding the other
                        num_temp.append('0')
                elif amino_acid == dna[0]:
                    num_temp[index] = '1'
                else:
                    if num_temp != ['0']:
                        num_temp[index] = '0'
                    if all(int(num)==0 for num in num_temp):
                        num_temp = ['0']
                if index+1 == len(num_temp):
                    max_num = max(int(num1) for num1 in num_temp)
                    temp.append(str(max_num))
                elif num_temp == ['0']:
                    temp.append('0')
        res.append(temp)
    return ' '.join(map(str,temp))

text13 = '''>Rosalind_9718
AAACAAAAACTTCCGTCCTATGCGCACTGTGCATCACGACCTTTTAAGATGATTGGCTTC
'''

with open('result5.txt','w') as file:
    file.write(motif_search(text13))

# def longest_subs(text):
#     max_num = float(text.splitlines()[0])
#     order = [num for num in text.splitlines()[1] if num!= ' ']
#     inc_res = [[num] for num in order if int(num)<=max_num/2]
#     dec_res = [[num] for num in order if int(num)>=max_num/2]
#     for no1, index1 in enumerate(inc_res):
#         no_in_order = order.index(index1[0])
#         check = order[no_in_order:]
#         temp = []
#         for no2, index2 in enumerate(check):
#             if no2+1 != len(check):
#                 no_checking = no2+1
#             else:
#                 if index2 > temp[-1]:
#                     temp.append(index2)
#                 break
#             if temp == []:
#                 num_base = int(order[no_in_order])
#             else:
#                 num_base = int(temp[-1])
#             no_base = check.index(str(num_base))
#             no_range = no_checking - no_base
#             num_range = int(check[no_checking]) - num_base
#             degree = no_range + num_range
#             degree_temp = [degree]
#             while no_range<min(degree_temp) and int(index2)>=int(index1[0]):
#                 no_range = no_checking - no_base
#                 num_range = int(check[no_checking]) - num_base
#                 degree = no_range + num_range
#                 degree_temp.append(degree)
#                 no_checking += 1
#             if len(degree_temp) > 1:
#                 if degree_temp[0] == degree_temp[1]:
#                     degree_temp.remove(degree_temp[0])
#             min_no = degree_temp.index(min(degree_temp))
#             temp.append(check[(no2+1):][min_no])
#             if len(temp) > 2:
#                 if temp[-1] == temp[-2]:
#                     temp.remove(temp[-1])
#         for num in temp:
#             inc_res[no1] += f'{ num}'
#     return inc_res

# def longest_subs(text):
#     max_num = float(text.splitlines()[0])
#     text.replace('\n', ' ')
#     order = [num for num in text.split(' ')[1:]]
#     inc_res = [[num] for num in order if int(num)<=max_num/2]
#     dec_res = [[num] for num in order if int(num)>=max_num/2]
#     for no1, index1 in enumerate(inc_res):
#         no_in_order = order.index(index1[0])
#         check = order[no_in_order:]
#         temp = []
#         no_checking = 1
#         no_base = 0
#         while no_checking < len(check):
#             no_range = no_checking - no_base
#             num_check = check[no_checking]
#             num_range = int(check[no_checking]) - int(check[no_base])
#             degree = no_range + num_range
#             degree_temp = [degree]
#             while no_range<min(degree_temp):
#                 no_range = no_checking - no_base
#                 if check[no_checking] > check[no_base]:
#                     num_range = int(check[no_checking]) - int(check[no_base])
#                 else:
#                     num_range = 1000000
#                 degree = no_range + num_range
#                 degree_temp.append(degree)
#                 no_checking += 1
#                 if no_checking == len(check):
#                     break
#             if len(degree_temp) > 1:
#                 target_index = degree_temp[1:].index(min(degree_temp[1:]))+1+no_base
#                 temp.append(str(check[target_index]))
#             else:
#                 break
#             no_base = target_index
#             no_checking = no_base+1
#             if no_base == len(check)-1:
#                 break
#         inc_res[no1] = inc_res[no1] + temp
#     for no1, index1 in enumerate(dec_res):   #changed
#         no_in_order = order.index(index1[0])
#         check = order[no_in_order:]
#         temp = []
#         no_checking = 1
#         no_base = 0
#         while no_checking < len(check):
#             no_range = no_checking - no_base
#             num_range = int(check[no_base]) - int(check[no_checking])  #changed
#             degree = no_range + num_range
#             degree_temp = [degree]
#             while no_range<min(degree_temp):
#                 no_range = no_checking - no_base
#                 if check[no_checking] < check[no_base]:
#                     num_range = int(check[no_base]) - int(check[no_checking]) #changed
#                 else:
#                     num_range = 1000000
#                 degree = no_range + num_range
#                 degree_temp.append(degree)
#                 no_checking += 1
#                 if no_checking == len(check):
#                     break
#             if len(degree_temp) > 1:
#                 target_index = degree_temp[1:].index(min(degree_temp[1:]))+1+no_base
#                 temp.append(str(check[target_index]))
#             else:
#                 break
#             no_base = target_index
#             no_checking = no_base+1
#             if no_base == len(check)-1:
#                 break
#         dec_res[no1] = dec_res[no1] + temp
#     inc_ans = ' '.join(max(inc_res, key=len))
#     dec_ans = ' '.join(max(dec_res, key=len))
#     ans = inc_ans + '\n' + dec_ans
#     return ans

# def longest_subs(text):
#     max_num = float(text.splitlines()[0])
#     text.replace('\n', ' ')
#     order = [int(num) for num in text.split(' ')[1:]]
#     res = [[num] for num in order]
#     for i in range(len(order)-1, -1, -1):
#         for j in range(i+1, len(order)):
#             if order[i] < order[j]:
#                 res[i] = order[i] + max(res[i], res[i]+)


# with open('rosalind_lgis.txt','r') as file:
#     text15 = file.read()
# print(longest_subs(text15))

res = [[1] for num in range(5)]
print(res)