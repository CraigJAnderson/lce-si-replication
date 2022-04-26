##take two discreetly phased trinucleotide counts, reverse complement those listed first and add them all together to get strand sensitive trinucleotide counts for each nodule listed
##combine_tnt_counts.py
##python 3
#use as python script.py inut1.tnt input2.tnt > output.tnt
import pandas as pd
import numpy as np
import sys

np.random.seed(1111)

##define function to produce complementary base
def revcom(dna):
 complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
 return ''.join([complement[base] for base in dna[::-1]])

#the following defines a specific order for trinucleotide output
bases=['G','T','A','C']
kmer=list() 
for i in bases:
 for j in bases:
  for k in bases:
   kmer.append(i+j+k)

a2n = pd.read_csv(sys.argv[1],sep="\t",header=[0],index_col=[0])
t2n = pd.read_csv(sys.argv[2],sep="\t",header=[0], index_col=[0])

##reverse complement all the names in the a2n df
new_names = []
for x in a2n.columns:
 new_names.append(str((revcom(x[2]))+(revcom(x[1]))+(revcom(x[0]))))
 
a2n.columns = new_names

##set up blank df based upon t2n order and add the same tnts from the same nodule together, storing them in the new df
df = t2n.replace(t2n, 0)
for x, row in a2n.iterrows():
 for z in row.index:
   df.at[x,z] = a2n.at[x,z] + t2n.at[x,z]

#correct for autosomes
from io import StringIO
output = StringIO()
df.to_csv(output,sep="\t",header=True, index=True)
output.seek(0)
print(output.read(),end="")
