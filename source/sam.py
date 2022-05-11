from cProfile import label
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


import matplotlib.pyplot as plt
import numpy as np

mapq_dict = {'aligned' : {'0-10':0, '10-20':1, '20-30':0, '30-40':0, '40-50':0, '50-60':0},\
                 'misaligned' : {'0-10':5, '10-20':0, '20-30':0, '30-40':2, '40-50':0, '50-60':0}}

print([(num/7*100) for num in mapq_dict['misaligned'].values()])




# make data:

labels = ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60']

bwa_mem_aligned =    [10, 10, 15, 10, 63, 434]
bwa_mem_misaligned = [400, 20, 50, 2, 5, 3]
bottoms_bwa_aligned = []
bottoms_bwa_misaligned = []
for i in range(len(bwa_mem_aligned)):
    item1=bwa_mem_aligned[i]
    item2=bwa_mem_misaligned[i]
    if item1 < item2:
        bottoms_bwa_misaligned.append(item1)
        bwa_mem_misaligned[i] -= item1
        bottoms_bwa_aligned.append(0)
    elif item1 > item2:
        bottoms_bwa_aligned.append(item2)
        bwa_mem_aligned[i] -= item2
        bottoms_bwa_misaligned.append(0)


bowtie = [4, 40, 50, 72, 80, 300]
label_loc = np.arange(len(labels))
width = 0.2

# plot
fig, ax = plt.subplots()


bwa_mem_aligned_rects = ax.bar(label_loc - width/2, bwa_mem_aligned, width=0.1, color="green", label = "BWA MEM aligned", bottom=bottoms_bwa_aligned)
bwa_mem_misaligned = ax.bar(label_loc - width/2, bwa_mem_misaligned, width=0.1, color="green", label = "BWA MEM misaligned", alpha=0.3, bottom=bottoms_bwa_misaligned)
bowtie_rects = ax.bar(label_loc + width/2, bowtie, width=0.1, color="orange", label = "Bowtie aligned")
ax.set_xlabel("Mapping quality")
ax.set_ylabel("Number of reads")
ax.set_title("Mapping quality distribution of aligned and misaligned reads")
ax.set_xticklabels(labels)
ax.set_xticks(label_loc)
ax.legend()



fig.tight_layout()

plt.show()



            
           
