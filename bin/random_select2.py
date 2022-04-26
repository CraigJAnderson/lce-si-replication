#python script to take a bed file and randomised index to create a new bed to bootstrap against
#use as: python random_select.py random_lines.txt test.bed out.bed
import pandas as pd #1.0.5
import numpy as np
import sys

in_bed = pd.read_csv(sys.argv[1], sep="\t",header=None, index_col=None)
out_bed_name = sys.argv[2]
num_variants = int(sys.argv[3])
out_bed = in_bed.iloc[np.random.randint(0,len(in_bed),size=num_variants)]

out_bed.to_csv(out_bed_name,sep="\t",header=False, index=False)

