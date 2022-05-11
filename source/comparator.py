import matplotlib.pyplot as plt
import numpy as np
import re
from alignment import parse_SAM_header

def calculate_mapping_quality_distribution(mapq_dict, mapq, align_type):
    if mapq < 10:
        mapq_dict[align_type]['0-10'] += 1
    elif (mapq >=10) and (mapq < 20):
        mapq_dict[align_type]['10-20'] += 1
    elif (mapq >=20) and (mapq < 30):
        mapq_dict[align_type]['20-30'] += 1
    elif (mapq >=30) and (mapq < 40):
        mapq_dict[align_type]['30-40'] += 1
    elif (mapq >=40) and (mapq < 50):
        mapq_dict[align_type]['40-50'] += 1
    else:
        mapq_dict[align_type]['50-60'] += 1 

    
def get_mapq_distribution_bar(mapq_distribution_dict, aligned_reads, misaligned_reads, unmapped_reads):
    
    bwa_mem_aligned = [round(num/(aligned_reads+misaligned_reads)*100) for num in mapq_distribution_dict['aligned'].values()]
    bwa_mem_misaligned = [round(num/(misaligned_reads+aligned_reads)*100) for num in mapq_distribution_dict['misaligned'].values()]
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
        else:
            bottoms_bwa_aligned.append(0)
            bottoms_bwa_misaligned.append(0) 
    labels = ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60'] # labels on x axis
    label_loc = np.arange(len(labels)) # positions of labels on x axis, equidistance
    width = 0.2 # width of the bar
    fig, ax = plt.subplots() # creates figure, plot and axis
    bwa_mem_aligned_rects = ax.bar(label_loc - width/2, bwa_mem_aligned, width=0.1, color="green", label = "BWA MEM aligned", bottom=bottoms_bwa_aligned) # draws bars on plot
    bwa_mem_misaligned_rects = ax.bar(label_loc - width/2, bwa_mem_misaligned, width=0.1, color="green", label = "BWA MEM misaligned", alpha=0.3, bottom=bottoms_bwa_misaligned)
    ax.set_xlabel("Mapping quality")
    ax.set_ylabel("Number of reads")
    ax.set_title("Mapping quality distribution of aligned and misaligned reads")
    ax.set_xticklabels(labels)
    ax.set_xticks(label_loc)
    #ax.bar_label(bwa_mem_aligned_rects, padding=3)
    #ax.bar_label(bwa_mem_misaligned_rects, padding=3)
    ax.legend()
    fig.tight_layout()
    plt.show()

def get_percision_bar():
    return
def get_recall_bar():
    return
def get_AUC_func():
    return
def get_fscore_func():
    return        
           

def compare(simulator_sam_file, tool_sam_file):
    unmapped_reads = 0
    misaligned_reads = 0
    aligned_reads = 0
    mapq_distribution_dict = {'aligned' : {'0-10':0, '10-20':0, '20-30':0, '30-40':0, '40-50':0, '50-60':0},\
                 'misaligned' : {'0-10':0, '10-20':0, '20-30':0, '30-40':0, '40-50':0, '50-60':0}}
    with open(simulator_sam_file, "r") as simulator_sam_handle, open(tool_sam_file, "r") as tool_sam_handle:
        sim_sam_header = parse_SAM_header(simulator_sam_handle)
        tool_sam_header = parse_SAM_header(tool_sam_handle)
        for sim_sam_line, tool_sam_line in zip(simulator_sam_handle, tool_sam_handle):
            sim_sam_line_rs = sim_sam_line.rstrip()
            tool_sam_line_rs = tool_sam_line.rstrip()
            sim_sam_fields = sim_sam_line_rs.split("\t")
            tool_sam_fields = tool_sam_line_rs.split("\t")
            sim_sam_alignment_pos = sim_sam_fields[3]
            tool_sam_alignment_pos = tool_sam_fields[3]
            tool_sam_alignment_flag = int(tool_sam_fields[1])
            tool_sam_alignment_mapq = int(tool_sam_fields[4])
            if (tool_sam_alignment_flag & 4) != 0:
                unmapped_reads += 1
            elif sim_sam_alignment_pos != tool_sam_alignment_pos:
                misaligned_reads += 1
                calculate_mapping_quality_distribution(mapq_distribution_dict, tool_sam_alignment_mapq, 'misaligned')
            else:
                aligned_reads += 1
                calculate_mapping_quality_distribution(mapq_distribution_dict, tool_sam_alignment_mapq, 'aligned')
        
        return(mapq_distribution_dict, aligned_reads, misaligned_reads, unmapped_reads)
            
        #print("misaligned: ", misaligned_reads)
        #print("correctly aligned: ", aligned_reads)
        #print("Unmapped reads: ", unmapped_reads)
        #print(mapq_distribution_dict.items())

def compare_(sim_sam_path, bwa_sam_path, bowtie_sam_path):
    mapq_distribution_dict_bwa, aligned_bwa, misaligned_bwa, unmapped_bwa = compare(sim_sam_path, bwa_sam_path)
    #mapq_distribution_dict_bowtie, aligned_bowtie, misaligned_bowtie, unmapped_bowtie = compare(sim_sam_path, bowtie_sam_path)

    # Data visualization, graphical representation
    get_mapq_distribution_bar(mapq_distribution_dict_bwa, aligned_bwa, misaligned_bwa, unmapped_bwa)
    #get_mapq_distribution_bar(mapq_distribution_dict_bowtie, aligned_bowtie, misaligned_bowtie, unmapped_bowtie)
    #get_percision_bar()
    #get_recall_bar()
    #get_AUC_func()
    #get_fscore_func()

    

simulator_sam_file = "./out_sam/ls_orchid.sam"
bwa_sam_file = "./bwa_sam/ls_orchid_cgc.sam"
bowtie_sam_file = "./bwa_sam/ls_orchid_cgc.sam"

compare_(simulator_sam_file, bwa_sam_file, bowtie_sam_file)