from dna_toolkit import *

test = {1:6, 2:5, 3:4}
result = test[max(test)]
print(result)

list1 = ['M','L']
list1[0] += 'M'
print(list1)

aa_seq = ['A', 'M', 'S', 'H', 'D', 'M', 'A', 'V', '_', 'X']
print(proteins_from_rf(aa_seq))
