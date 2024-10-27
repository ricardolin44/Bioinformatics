def rabrec(n,k):
    seq = [1,1]
    for times in range(n-2):
        Fn = seq[times]*k+seq[times+1]
        seq.append(Fn)
    return Fn

def rabrec2(n,k):
    child, parent = 1, 1
    for itr in range(n-1):
        child, parent = parent ,parent+ child*k
    return child

# print(rabrec(29,5))
print(rabrec2(29,5))

