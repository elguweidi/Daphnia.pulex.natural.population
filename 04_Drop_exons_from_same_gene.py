# This code was used to keep the longest exon from each gene
import pandas as pd
from Bio import SeqIO
myList = []
for seq_record in SeqIO.parse("scgFile.fa", "fasta"):
    myList.append([seq_record.id, str(seq_record.seq), len(seq_record)])

myList = pd.DataFrame(myList)

df = (myList.sort_values([0,2], ascending=[True, False])).drop_duplicates([0]).reset_index(drop=True)

scgFile = open("scgenes.fa", 'w')
for i in range(len(df)):
    print(">" +  str(df.loc[i,0]), file=scgFile)
    print(str(df.loc[i, 1]), file=scgFile)
scgFile.close()
