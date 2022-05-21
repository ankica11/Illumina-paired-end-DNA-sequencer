
from tkinter import font
import matplotlib.pyplot as plt
import numpy as np
from comparator import compare
import pandas as pd

from statistics import get_statistics



def get_accuracy_histogram(bwa_data, bowtie_data, genome_name):

    #bwa_data = [(0,1,2,3), (0,1,2,3), (0,1,2,3), (0,1,2,3), (0,1,2,3)]
    
    tools = ["BWA MEM", "Bowtie"]
   
    bwa_data_low = bwa_data[0]
    bwa_data_medium = bwa_data[1]
    bwa_data_high = bwa_data[2]
    bowtie_data_low = bowtie_data[0]
    bowtie_data_medium = bowtie_data[1]
    bowtie_data_high = bowtie_data[2]

    total = bwa_data[0][1] + bwa_data[0][2] + bwa_data[0][3]

    bwa_mem_aligned = [round(item[1]/total*100, 2) for item in bwa_data]
    bwa_mem_misaligned = [round(item[2]/total*100, 2) for item in bwa_data]
    bwa_mem_unmapped = [round(item[3]/total*100, 2) for item in bwa_data]

    bowtie_aligned = [round(item[1]/total*100, 2) for item in bowtie_data]
    bowtie_misaligned = [round(item[2]/total*100, 2) for item in bowtie_data]
    bowtie_unmapped = [round(item[3]/total*100, 2) for item in bowtie_data]

    
    width = 0.13
    labels = ['Set 1', 'Set 2', 'Set 3', 'Set 4', 'Set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_aligned_rects = ax.bar(label_loc - 3*width, bwa_mem_aligned, width=0.1, color="#F24BC0", label="BWA MEM aligned")
    bowtie_aligned_rects = ax.bar(label_loc - 2*width, bowtie_aligned, width=0.1, color="#3666AC", label="Bowtie aligned")
   
    bwa_mem_misaligned_rects = ax.bar(label_loc - width, bwa_mem_misaligned, width=0.1, color="#F24BC0", label="BWA MEM misaligned", alpha=0.6)
    bowtie_misaligned_rects = ax.bar(label_loc, bowtie_misaligned, width=0.1, color="#3666AC", label="Bowtie misaligned", alpha=0.6)
  
    bwa_mem_unmapped_rects = ax.bar(label_loc + width, bwa_mem_unmapped, width=0.1, color="#F24BC0", label="BWA MEM unmapped", alpha=0.2)
    bowtie_unmapped_rects = ax.bar(label_loc + 2*width, bowtie_unmapped, width=0.1, color="#3666AC", label="Bowtie unmapped", alpha=0.2)
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
    return (bwa_mem_aligned, bwa_mem_misaligned, bwa_mem_unmapped, bowtie_aligned, bowtie_misaligned, bowtie_unmapped)

def get_precision_histogram(bwa_data, bowtie_data, genome_name):
    

    true_positives_bwa_mem = [item[1] for item in bwa_data]
    true_positives_bowtie =  [item[1] for item in bowtie_data]
    false_positives_bwa_mem = [item[2] for item in bwa_data]
    false_positives_bowtie = [item[2] for item in bowtie_data]

    #precision = TP/(TP+FP)
    bwa_mem_precision = [round(true_positive/(true_positive + false_positive)*100, 2) for true_positive, false_positive in zip(true_positives_bwa_mem, false_positives_bwa_mem)]
    bowtie_precision = [round(true_positive/(true_positive + false_positive)*100, 2) for true_positive, false_positive in zip(true_positives_bowtie, false_positives_bowtie)]
  
    width = 0.15
    labels = ['Set 1', 'Set 2', 'Set 3', 'Set 4', 'Set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_precision_rects = ax.bar(label_loc - width, bwa_mem_precision, width=0.2, color="#F24BC0", label="BWA MEM")
    bowtie_precision_rects = ax.bar(label_loc + width, bowtie_precision, width=0.2, color="#3666AC", label="Bowtie")
   
    ax.set_ylabel("Precision (%)", fontsize=14)
    ax.set_title("Precision of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_precision_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_precision_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()
    return bwa_mem_precision, bowtie_precision



def get_recall_histogram(bwa_data, bowtie_data, genome_name):
    true_positives_bwa_mem = [item[1] for item in bwa_data]
    true_positives_bowtie = [item[1] for item in bowtie_data]
    false_negatives_bwa_mem = [item[3] for item in bwa_data]
    false_negatives_bowtie = [item[3] for item in bowtie_data]

    #recall = TP/(TP+FN)
    bwa_mem_recall = [round(true_positive/(true_positive + false_negative)*100, 2) for true_positive, false_negative in zip(true_positives_bwa_mem, false_negatives_bwa_mem)]
    bowtie_recall = [round(true_positive/(true_positive + false_negative)*100, 2) for true_positive, false_negative in zip(true_positives_bowtie, false_negatives_bowtie)]
   
    
    width = 0.15
    labels = ['Set 1', 'Set 2', 'Set 3', 'Set 4', 'Set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_recall_rects = ax.bar(label_loc - width, bwa_mem_recall, width=0.2, color="#F24BC0", label="BWA MEM")
    bowtie_recall_rects = ax.bar(label_loc + width, bowtie_recall, width=0.2, color="#3666AC", label="Bowtie")
   
    ax.set_ylabel("Recall (%)", fontsize=14)
    ax.set_title("Recall (sensitivity) of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_recall_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_recall_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()
    return bwa_mem_recall, bowtie_recall


def get_fcore_histogram(bwa_data, bowtie_data, genome_name):
    true_positives_bwa_mem = [item[1] for item in bwa_data]
    true_positives_bowtie = [item[1] for item in bowtie_data]
    false_negatives_bwa_mem = [item[3] for item in bwa_data]
    false_negatives_bowtie = [item[3] for item in bowtie_data]
    false_positives_bwa_mem = [item[2] for item in bwa_data]
    false_positives_bowtie = [item[2] for item in bowtie_data]

    bwa_mem_recall = [round(true_positive/(true_positive + false_negative)*100, 2) for true_positive, false_negative in zip(true_positives_bwa_mem, false_negatives_bwa_mem)]
    bowtie_recall = [round(true_positive/(true_positive + false_negative)*100, 2) for true_positive, false_negative in zip(true_positives_bowtie, false_negatives_bowtie)]
    bwa_mem_precision = [round(true_positive/(true_positive + false_positive)*100, 2) for true_positive, false_positive in zip(true_positives_bwa_mem, false_positives_bwa_mem)]
    bowtie_precision = [round(true_positive/(true_positive + false_positive)*100, 2) for true_positive, false_positive in zip(true_positives_bowtie, false_positives_bowtie)]
   

    bwa_mem_fscore = [ round(2*precision*recall/(precision+recall), 2) for precision, recall in zip(bwa_mem_precision, bwa_mem_recall)]
    bowtie_fscore = [ round (2*precision*recall/(precision+recall), 2) for precision, recall in zip(bowtie_precision, bowtie_recall)]
    
    
    width = 0.15
    labels = ['Set 1', 'Set 2', 'Set 3', 'Set 4', 'Set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_fscore_rects = ax.bar(label_loc - width, bwa_mem_fscore, width=0.2, color="#F24BC0", label="BWA MEM")
    bowtie_fscore_rects = ax.bar(label_loc + width, bowtie_fscore, width=0.2, color="#3666AC", label="Bowtie")
   
    ax.set_ylabel("F_score (%)", fontsize=14)
    ax.set_title("F_score of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_fscore_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_fscore_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()
    return bwa_mem_fscore, bowtie_fscore
    
         

def compare_multiple():
    directory_path = "./testing2/Clostridium/"
    filename = "Clostridium_tetani"
    genome_name = "Clostridium tetani"
    sim_sam_paths = [directory_path+str(i)+filename+".sam" for i in range(5)]
    bwa_sam_paths = [directory_path+str(i)+filename+"_bwa_mem.sam" for i in range(5)]
    bowtie_sam_paths = [directory_path+str(i)+filename+"_bowtie.sam" for i in range(5)]
    bwa_data = []
    bowtie_data = []  

    for i in range(5):
        bwa_data.append(compare(sim_sam_paths[i], bwa_sam_paths[i], 10))
        bowtie_data.append(compare(sim_sam_paths[i], bowtie_sam_paths[i], 7))

 
   
    bw_aln, bw_mis, bw_un, bow_aln, bow_mis, bow_un = get_accuracy_histogram(bwa_data, bowtie_data, genome_name)
    bw_pre, bow_pre = get_precision_histogram(bwa_data, bowtie_data, genome_name)
    bw_rec, bow_rec = get_recall_histogram(bwa_data, bowtie_data, genome_name)
    bw_fs, bow_fs = get_fcore_histogram(bwa_data, bowtie_data, genome_name)

    path = "./testing2/Clostridium/Clostridium_tetani.csv"
    get_statistics(bw_aln, bow_aln, bw_mis, bow_mis, bw_un, bow_un, bw_pre, bow_pre, bw_rec, bow_rec, bw_fs, bow_fs, path)



compare_multiple()




