
import matplotlib.pyplot as plt
import numpy as np
from comparator import compare


def get_accuracy_histogram(bwa_data, bowtie_data, genome_name):
    
    tools = ["BWA MEM", "Bowtie"]
   
    bwa_data_low = bwa_data[0]
    bwa_data_medium = bwa_data[1]
    bwa_data_high = bwa_data[2]
    bowtie_data_low = bowtie_data[0]
    bowtie_data_medium = bowtie_data[1]
    bowtie_data_high = bowtie_data[2]
    
    bwa_mem_aligned = [round(bwa_data_low[1]/(bwa_data_low[1]+bwa_data_low[2]+bwa_data_low[3])*100, 2),\
                       round(bwa_data_medium[1]/(bwa_data_medium[1]+bwa_data_medium[2]+bwa_data_medium[3])*100, 2),\
                       round(bwa_data_high[1]/(bwa_data_high[1]+bwa_data_high[2]+bwa_data_high[3])*100, 2)]
    bwa_mem_misaligned =  [round(bwa_data_low[2]/(bwa_data_low[1]+bwa_data_low[2]+bwa_data_low[3])*100, 2),\
                           round(bwa_data_medium[2]/(bwa_data_medium[1]+bwa_data_medium[2]+bwa_data_medium[3])*100, 2),\
                           round(bwa_data_high[2]/(bwa_data_high[1]+bwa_data_high[2]+bwa_data_high[3])*100, 2)]
    
    bwa_mem_unmapped = [round(bwa_data_low[3]/(bwa_data_low[1]+bwa_data_low[2]+bwa_data_low[3])*100, 2),\
                       round(bwa_data_medium[3]/(bwa_data_medium[1]+bwa_data_medium[2]+bwa_data_medium[3])*100, 2),\
                       round(bwa_data_high[3]/(bwa_data_high[1]+bwa_data_high[2]+bwa_data_high[3])*100, 2)]

    bowtie_aligned = [round(bowtie_data_low[1]/(bowtie_data_low[1]+bowtie_data_low[2]+bowtie_data_low[3])*100, 2),\
                       round(bowtie_data_medium[1]/(bowtie_data_medium[1]+bowtie_data_medium[2]+bowtie_data_medium[3])*100, 2),\
                       round(bowtie_data_high[1]/(bowtie_data_high[1]+bowtie_data_high[2]+bowtie_data_high[3])*100, 2)]
    bowtie_misaligned =  [round(bowtie_data_low[2]/(bowtie_data_low[1]+bowtie_data_low[2]+bowtie_data_low[3])*100, 2),\
                           round(bowtie_data_medium[2]/(bowtie_data_medium[1]+bowtie_data_medium[2]+bowtie_data_medium[3])*100, 2),\
                           round(bowtie_data_high[2]/(bowtie_data_high[1]+bowtie_data_high[2]+bowtie_data_high[3])*100, 2)]
    
    bowtie_unmapped = [round(bowtie_data_low[3]/(bowtie_data_low[1]+bowtie_data_low[2]+bowtie_data_low[3])*100, 2),\
                       round(bowtie_data_medium[3]/(bowtie_data_medium[1]+bowtie_data_medium[2]+bowtie_data_medium[3])*100, 2),\
                       round(bowtie_data_high[3]/(bowtie_data_high[1]+bowtie_data_high[2]+bowtie_data_high[3])*100, 2)]
    
    
    width = 0.08
    labels = ['errors set 1', 'errors set 2', 'errors set 3']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_aligned_rects = ax.bar(label_loc - 3*width, bwa_mem_aligned, width=0.02, color="green", label="BWA MEM aligned")
    bowtie_aligned_rects = ax.bar(label_loc - 2*width, bowtie_aligned, width=0.02, color="orange", label="Bowtie aligned")
   
    bwa_mem_misaligned_rects = ax.bar(label_loc - width, bwa_mem_misaligned, width=0.02, color="green", label="BWA MEM misaligned", alpha=0.6)
    bowtie_misaligned_rects = ax.bar(label_loc, bowtie_misaligned, width=0.02, color="orange", label="Bowtie misaligned", alpha=0.6)
  
    bwa_mem_unmapped_rects = ax.bar(label_loc + width, bwa_mem_unmapped, width=0.02, color="green", label="BWA MEM unmapped", alpha=0.2)
    bowtie_unmapped_rects = ax.bar(label_loc + 2*width, bowtie_unmapped, width=0.02, color="orange", label="Bowtie unmapped", alpha=0.2)
    ax.set_ylabel("Accuracy (%)", fontsize=14)
    ax.set_title("Percentage of aligned, misaligned and unmapped reads for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_aligned_rects, padding=3)
    ax.bar_label(bwa_mem_misaligned_rects, padding=3)
    ax.bar_label(bwa_mem_unmapped_rects, padding=3)
    ax.bar_label(bowtie_aligned_rects, padding=3)
    ax.bar_label(bowtie_misaligned_rects, padding=3)
    ax.bar_label(bowtie_unmapped_rects, padding=3)
    ax.legend()
    fig.tight_layout()
    plt.show()


def get_precision_histogram(bwa_data, bowtie_data, genome_name):
    

    true_positives_bwa_mem = (bwa_data[0][1], bwa_data[1][1], bwa_data[2][1])
    true_positives_bowtie = (bowtie_data[0][1], bowtie_data[1][1], bowtie_data[2][1]) # correctly aligned reads in 3 sets of errors

    false_positives_bwa_mem = (bwa_data[0][2], bwa_data[1][2], bwa_data[2][2])
    false_positives_bowtie = (bowtie_data[0][2], bowtie_data[1][2], bowtie_data[2][2]) # misaligned reads in 3 sets of errors


    #precision = TP/(TP+FP)
    bwa_mem_precision = [round(true_positives_bwa_mem[0]/(true_positives_bwa_mem[0]+false_positives_bwa_mem[0])*100, 2),\
                         round(true_positives_bwa_mem[1]/(true_positives_bwa_mem[1]+false_positives_bwa_mem[1])*100, 2),\
                         round(true_positives_bwa_mem[2]/(true_positives_bwa_mem[2]+false_positives_bwa_mem[2])*100, 2)]
    bowtie_precision = [round(true_positives_bowtie[0]/(true_positives_bowtie[0]+false_positives_bowtie[0])*100, 2),\
                        round(true_positives_bowtie[1]/(true_positives_bowtie[1]+false_positives_bowtie[1])*100, 2),\
                        round(true_positives_bowtie[2]/(true_positives_bowtie[2]+false_positives_bowtie[2])*100, 2)]
    
    width = 0.08
    labels = ['errors set 1', 'errors set 2', 'errors set 3']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_precision_rects = ax.bar(label_loc - width, bwa_mem_precision, width=0.06, color="green", label="BWA MEM")
    bowtie_precision_rects = ax.bar(label_loc + width, bowtie_precision, width=0.06, color="orange", label="Bowtie")
   
    ax.set_ylabel("Precision (%)", fontsize=14)
    ax.set_title("Precision of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_precision_rects, padding=3)
    ax.bar_label(bowtie_precision_rects, padding=3)
    ax.legend()
    fig.tight_layout()
    plt.show()



def get_recall_histogram(bwa_data, bowtie_data, genome_name):
    true_positives_bwa_mem = (bwa_data[0][1], bwa_data[1][1], bwa_data[2][1])
    true_positives_bowtie = (bowtie_data[0][1], bowtie_data[1][1], bowtie_data[2][1]) # correctly aligned reads in 3 sets of errors

    false_negatives_bwa_mem = (bwa_data[0][3], bwa_data[1][3], bwa_data[2][3])
    false_negatives_bowtie = (bowtie_data[0][3], bowtie_data[1][3], bowtie_data[2][3]) # unmapped reads in 3 sets of errors


    #recall = TP/(TP+FN)
    bwa_mem_recall = [round(true_positives_bwa_mem[0]/(true_positives_bwa_mem[0]+false_negatives_bwa_mem[0])*100, 2),\
                         round(true_positives_bwa_mem[1]/(true_positives_bwa_mem[1]+false_negatives_bwa_mem[1])*100, 2),\
                         round(true_positives_bwa_mem[2]/(true_positives_bwa_mem[2]+false_negatives_bwa_mem[2])*100, 2)]
    bowtie_recall = [round(true_positives_bowtie[0]/(true_positives_bowtie[0]+false_negatives_bowtie[0])*100, 2),\
                        round(true_positives_bowtie[1]/(true_positives_bowtie[1]+false_negatives_bowtie[1])*100, 2),\
                        round(true_positives_bowtie[2]/(true_positives_bowtie[2]+false_negatives_bowtie[2])*100, 2)]
    
    width = 0.08
    labels = ['errors set 1', 'errors set 2', 'errors set 3']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_recall_rects = ax.bar(label_loc - width, bwa_mem_recall, width=0.06, color="green", label="BWA MEM")
    bowtie_recall_rects = ax.bar(label_loc + width, bowtie_recall, width=0.06, color="orange", label="Bowtie")
   
    ax.set_ylabel("Recall (%)", fontsize=14)
    ax.set_title("Recall (sensitivity) of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_recall_rects, padding=3)
    ax.bar_label(bowtie_recall_rects, padding=3)
    ax.legend()
    fig.tight_layout()
    plt.show()
    


def get_fcore_histogram(bwa_data, bowtie_data, genome_name):
    true_positives_bwa_mem = (bwa_data[0][1], bwa_data[1][1], bwa_data[2][1])
    true_positives_bowtie = (bowtie_data[0][1], bowtie_data[1][1], bowtie_data[2][1]) # correctly aligned reads in 3 sets of errors

    false_negatives_bwa_mem = (bwa_data[0][3], bwa_data[1][3], bwa_data[2][3])
    false_negatives_bowtie = (bowtie_data[0][3], bowtie_data[1][3], bowtie_data[2][3]) # unmapped reads in 3 sets of errors

    false_positives_bwa_mem = (bwa_data[0][2], bwa_data[1][2], bwa_data[2][2])
    false_positives_bowtie = (bowtie_data[0][2], bowtie_data[1][2], bowtie_data[2][2]) # misaligned reads in 3 sets of errors



    #recall = TP/(TP+FN)
    bwa_mem_recall = [round(true_positives_bwa_mem[0]/(true_positives_bwa_mem[0]+false_negatives_bwa_mem[0])*100, 2),\
                         round(true_positives_bwa_mem[1]/(true_positives_bwa_mem[1]+false_negatives_bwa_mem[1])*100, 2),\
                         round(true_positives_bwa_mem[2]/(true_positives_bwa_mem[2]+false_negatives_bwa_mem[2])*100, 2)]
    bowtie_recall = [round(true_positives_bowtie[0]/(true_positives_bowtie[0]+false_negatives_bowtie[0])*100, 2),\
                        round(true_positives_bowtie[1]/(true_positives_bowtie[1]+false_negatives_bowtie[1])*100, 2),\
                        round(true_positives_bowtie[2]/(true_positives_bowtie[2]+false_negatives_bowtie[2])*100, 2)]

    bwa_mem_precision = [round(true_positives_bwa_mem[0]/(true_positives_bwa_mem[0]+false_positives_bwa_mem[0])*100, 2),\
                         round(true_positives_bwa_mem[1]/(true_positives_bwa_mem[1]+false_positives_bwa_mem[1])*100, 2),\
                         round(true_positives_bwa_mem[2]/(true_positives_bwa_mem[2]+false_positives_bwa_mem[2])*100, 2)]
    bowtie_precision = [round(true_positives_bowtie[0]/(true_positives_bowtie[0]+false_positives_bowtie[0])*100, 2),\
                        round(true_positives_bowtie[1]/(true_positives_bowtie[1]+false_positives_bowtie[1])*100, 2),\
                        round(true_positives_bowtie[2]/(true_positives_bowtie[2]+false_positives_bowtie[2])*100, 2)]

    bwa_mem_fscore = [ round(2*precision*recall/(precision+recall), 2) for precision, recall in zip(bwa_mem_precision, bwa_mem_recall)]
    bowtie_fscore = [ round (2*precision*recall/(precision+recall), 2) for precision, recall in zip(bowtie_precision, bowtie_recall)]
    
    
    width = 0.08
    labels = ['errors set 1', 'errors set 2', 'errors set 3']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_fscore_rects = ax.bar(label_loc - width, bwa_mem_fscore, width=0.06, color="green", label="BWA MEM")
    bowtie_fscore_rects = ax.bar(label_loc + width, bowtie_fscore, width=0.06, color="orange", label="Bowtie")
   
    ax.set_ylabel("F_score (%)", fontsize=14)
    ax.set_title("F_score of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_fscore_rects, padding=3)
    ax.bar_label(bowtie_fscore_rects, padding=3)
    ax.legend()
    fig.tight_layout()
    plt.show()
    
         

def compare_multiple():
    genome_name="Escherichia phage"
    sim_sam_paths=["./testing/Staphylococcus borealis/1/GCF_003042555.1_ASM304255v1_genomic.sam","./testing/Staphylococcus borealis/2/GCF_003042555.1_ASM304255v1_genomic.sam","./testing/Staphylococcus borealis/3/GCF_003042555.1_ASM304255v1_genomic.sam"] 
    bwa_sam_paths=["./testing/Staphylococcus borealis/1/GCF_003042555.1_ASM304255v1_genomic_bwa_mem.sam","./testing/Staphylococcus borealis/2/GCF_003042555.1_ASM304255v1_genomic_bwa_mem.sam","./testing/Staphylococcus borealis/3/GCF_003042555.1_ASM304255v1_genomic_bwa_mem.sam"] 
    bowtie_sam_paths=["./testing/Staphylococcus borealis/1/GCF_003042555.1_ASM304255v1_genomic_bowtie.sam","./testing/Staphylococcus borealis/2/GCF_003042555.1_ASM304255v1_genomic_bowtie.sam", "./testing/Staphylococcus borealis/3/GCF_003042555.1_ASM304255v1_genomic_bowtie.sam"] 
    bwa_data = []
    bowtie_data = []  

    for i in range(3):
        bwa_data.append(compare(sim_sam_paths[i], bwa_sam_paths[i], 10))
        bowtie_data.append(compare(sim_sam_paths[i], bowtie_sam_paths[i], 7))
   
    get_accuracy_histogram(bwa_data, bowtie_data, genome_name)
    get_precision_histogram(bwa_data, bowtie_data, genome_name)
    get_recall_histogram(bwa_data, bowtie_data, genome_name)
    get_fcore_histogram(bwa_data, bowtie_data, genome_name)


compare_multiple()


