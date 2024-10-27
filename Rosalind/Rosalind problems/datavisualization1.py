import numpy as np
import pandas as pd

data = np.genfromtxt("hg19.cage_peak_phase1and2combined_tpm_ann_decoded.osc.txt.gz.extract.tsv",
                     comments="#", usecols=range(2,73,1), names=True, dtype=object, delimiter="\t")
df = pd.DataFrame(data)
print("Number of genes: " + str(len(df)))
print(df.head())