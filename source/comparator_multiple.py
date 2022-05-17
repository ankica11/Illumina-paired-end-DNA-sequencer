
from tkinter import font
import matplotlib.pyplot as plt
import numpy as np
from comparator import compare


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

    
    width = 0.15
    labels = ['set 1', 'set 2', 'set 3', 'set 4', 'set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_aligned_rects = ax.bar(label_loc - 3*width, bwa_mem_aligned, width=0.1, color="green", label="BWA MEM aligned")
    bowtie_aligned_rects = ax.bar(label_loc - 2*width, bowtie_aligned, width=0.1, color="orange", label="Bowtie aligned")
   
    bwa_mem_misaligned_rects = ax.bar(label_loc - width, bwa_mem_misaligned, width=0.1, color="green", label="BWA MEM misaligned", alpha=0.6)
    bowtie_misaligned_rects = ax.bar(label_loc, bowtie_misaligned, width=0.1, color="orange", label="Bowtie misaligned", alpha=0.6)
  
    bwa_mem_unmapped_rects = ax.bar(label_loc + width, bwa_mem_unmapped, width=0.1, color="green", label="BWA MEM unmapped", alpha=0.2)
    bowtie_unmapped_rects = ax.bar(label_loc + 2*width, bowtie_unmapped, width=0.1, color="orange", label="Bowtie unmapped", alpha=0.2)
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
    

    true_positives_bwa_mem = [item[1] for item in bwa_data]
    true_positives_bowtie =  [item[1] for item in bowtie_data]
    false_positives_bwa_mem = [item[2] for item in bwa_data]
    false_positives_bowtie = [item[2] for item in bowtie_data]

    #precision = TP/(TP+FP)
    bwa_mem_precision = [round(true_positive/(true_positive + false_positive)*100, 2) for true_positive, false_positive in zip(true_positives_bwa_mem, false_positives_bwa_mem)]
    bowtie_precision = [round(true_positive/(true_positive + false_positive)*100, 2) for true_positive, false_positive in zip(true_positives_bowtie, false_positives_bowtie)]
  
    width = 0.15
    labels = ['set 1', 'set 2', 'set 3', 'set 4', 'set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_precision_rects = ax.bar(label_loc - width, bwa_mem_precision, width=0.2, color="green", label="BWA MEM")
    bowtie_precision_rects = ax.bar(label_loc + width, bowtie_precision, width=0.2, color="orange", label="Bowtie")
   
    ax.set_ylabel("Precision (%)", fontsize=14)
    ax.set_title("Precision of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_precision_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_precision_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()



def get_recall_histogram(bwa_data, bowtie_data, genome_name):
    true_positives_bwa_mem = [item[1] for item in bwa_data]
    true_positives_bowtie = [item[1] for item in bowtie_data]
    false_negatives_bwa_mem = [item[3] for item in bwa_data]
    false_negatives_bowtie = [item[3] for item in bowtie_data]

    #recall = TP/(TP+FN)
    bwa_mem_recall = [round(true_positive/(true_positive + false_negative)*100, 2) for true_positive, false_negative in zip(true_positives_bwa_mem, false_negatives_bwa_mem)]
    bowtie_recall = [round(true_positive/(true_positive + false_negative)*100, 2) for true_positive, false_negative in zip(true_positives_bowtie, false_negatives_bowtie)]
   
    
    width = 0.15
    labels = ['set 1', 'set 2', 'set 3', 'set 4', 'set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_recall_rects = ax.bar(label_loc - width, bwa_mem_recall, width=0.2, color="green", label="BWA MEM")
    bowtie_recall_rects = ax.bar(label_loc + width, bowtie_recall, width=0.2, color="orange", label="Bowtie")
   
    ax.set_ylabel("Recall (%)", fontsize=14)
    ax.set_title("Recall (sensitivity) of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_recall_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_recall_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()
    


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
    labels = ['set 1', 'set 2', 'set 3', 'set 4', 'set 5']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_fscore_rects = ax.bar(label_loc - width, bwa_mem_fscore, width=0.2, color="green", label="BWA MEM")
    bowtie_fscore_rects = ax.bar(label_loc + width, bowtie_fscore, width=0.2, color="orange", label="Bowtie")
   
    ax.set_ylabel("F_score (%)", fontsize=14)
    ax.set_title("F_score of aligners for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_fscore_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_fscore_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()
    
         

def compare_multiple():
    genome_name="Clostridium tetani"
    sim_sam_paths=["./testing2/Clostridium/4/0Clostridium_tetani.sam", "./testing2/Clostridium/4/1Clostridium_tetani.sam", "./testing2/Clostridium/4/2Clostridium_tetani.sam", "./testing2/Clostridium/4/3Clostridium_tetani.sam", "./testing2/Clostridium/4/4Clostridium_tetani.sam"] 
    bwa_sam_paths=["./testing2/Clostridium/4/0Clostridium_tetani_bwa_mem.sam", "./testing2/Clostridium/4/1Clostridium_tetani_bwa_mem.sam", "./testing2/Clostridium/4/2Clostridium_tetani_bwa_mem.sam", "./testing2/Clostridium/4/3Clostridium_tetani_bwa_mem.sam", "./testing2/Clostridium/4/4Clostridium_tetani_bwa_mem.sam"]
    bowtie_sam_paths=["./testing2/Clostridium/4/0Clostridium_tetani_bwa_mem.sam", "./testing2/Clostridium/4/1Clostridium_tetani_bowtie.sam", "./testing2/Clostridium/4/2Clostridium_tetani_bowtie.sam", "./testing2/Clostridium/4/3Clostridium_tetani_bowtie.sam", "./testing2/Clostridium/4/4Clostridium_tetani_bowtie.sam"]
    bwa_data = []
    bowtie_data = []  

    for i in range(5):
        bwa_data.append(compare(sim_sam_paths[i], bwa_sam_paths[i], 10))
        bowtie_data.append(compare(sim_sam_paths[i], bowtie_sam_paths[i], 7))
   
    get_accuracy_histogram(bwa_data, bowtie_data, genome_name)
    get_precision_histogram(bwa_data, bowtie_data, genome_name)
    get_recall_histogram(bwa_data, bowtie_data, genome_name)
    get_fcore_histogram(bwa_data, bowtie_data, genome_name)


compare_multiple()


