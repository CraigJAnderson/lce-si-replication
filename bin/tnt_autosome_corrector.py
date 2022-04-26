##take two discreetly phased trinucleotide counts, reverse complement those listed first and add them all together to get strand sensitive trinucleotide counts for each nodule listed
##combine_tnt_counts.py
##python 3
#use as python script.py inut1.tnt input2.tnt > output.tnt
import sys
import pandas as pd

a2n = pd.read_csv(sys.argv[1],sep="\t",header=[0],index_col=[0])

#correct for autosomes
a2n = a2n*2

##set up blank df based upon t2n order and add the same tnts from the same nodule together, storing them in the new df
from io import StringIO
output = StringIO()
a2n.to_csv(output,sep="\t",header=True, index=True)
output.seek(0)
print(output.read(),end="")
