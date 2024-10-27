def readFile(filePath):
    with open(filePath, 'r') as f:
        return [l.strip() for l in f.readlines()]
    
def gccontent(seq):
    return round((seq.count('C')+seq.count('G'))/len(seq)*100,3)
    
FASTAFile = readFile('Rosalind problems/test_data/gc_content.txt')
FASTADict = {}
FASTALabel = ''

for line in FASTAFile:
    if '>' in line:
        FASTALabel = line
        FASTADict[FASTALabel] = ""
    else:
        FASTADict[FASTALabel] += line

print(FASTAFile)
print(FASTADict)

ResultDict = {key: gccontent(value) for (key, value) in FASTADict.items()}
print(ResultDict)
MAXGCKey = max(ResultDict, key=ResultDict.get)
print(f'{MAXGCKey[1:]}\n{ResultDict[MAXGCKey]}')
