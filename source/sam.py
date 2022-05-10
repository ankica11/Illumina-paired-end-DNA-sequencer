from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import os

seq1= Seq("ATTTTT")
seq2 = Seq("GCCCCC")

seq3 = seq1[0:3] + seq2[3:6]

#print(type(seq3))
#path="D:/1master/1gi/projekat/GCF_000861905.1_ViralProj15445_genomic.fasta"

#print(os.path.split(path)[1])




            
           
