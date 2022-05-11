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



mapq_dict = {'aligned' : {'0-10':0, '10-20':1, '20-30':0, '30-40':0, '40-50':0, '50-60':0},\
                 'misaligned' : {'0-10':0, '10-20':0, '20-30':0, '30-40':2, '40-50':0, '50-60':0}}

print(mapq_dict['misaligned']['30-40'])
            
           
