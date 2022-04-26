#python script to calculate mutation rate, given matching tables for mutation counts and trinucleotide counts
#use as: python matrix_mu_rate.py tnt.bed mut.bed > out.bed
import pandas as pd
import sys
import numpy as np

pd.options.mode.chained_assignment = None  # default='warn'

tnt = pd.read_csv(sys.argv[1],sep="\t",header=0,index_col=0)
rfd_bin = sys.argv[4]
tnt = tnt.filter(like=str(rfd_bin), axis=0)
mut = pd.read_csv(sys.argv[2],sep="\t",header=0,index_col=0)
mut = mut.filter(like=str(rfd_bin), axis=0)
norm = pd.read_csv(sys.argv[3],sep="\t",header=0,index_col=None)
#tnt = pd.read_csv("test.tri",sep="\t",header=0,index_col=None)
#mut = pd.read_csv("test.mut",sep="\t",header=0,index_col=None)

out = mut.iloc[:,0:1]
out['mu'] = ((mut.iloc[:,0:64]/tnt.iloc[:,0:64])*norm.iloc[0,0:64]).sum(axis=1)
#out['mu'] = ((((mut.iloc[:,4:69]/tnt.iloc[:,4:69])*norm.iloc[0,0:65]).sum(axis=1)*tnt.iloc[:,69])/(tnt.iloc[:,68]+tnt.iloc[:,69]))/tnt.iloc[:,70]
out['mu'].replace([np.inf], "NA", inplace=True)
#out = out[((out.end-out.start) > 990) & ((out.end-out.start) < 1010) & (tnt.XXX < 10)]
out = out.iloc[:,1:3]
from io import StringIO
output = StringIO()
out.to_csv(output,sep="\t",header=False, index=True)
output.seek(0)
print(output.read(),end="")

