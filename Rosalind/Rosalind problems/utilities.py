def readTextFile(filePath):
    with open(filePath, 'r') as f:
        return ''.join([l.strip() for l in f.readlines()])
    
def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq + '\n')

def read_FASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAFILE = [l.strip() for l in f.readlines()]
    
    FASTADict = {}
    FASTALabel = ''

    for line in FASTAFILE:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ''
        else:
            FASTADict[FASTALabel] += line

    return FASTADict