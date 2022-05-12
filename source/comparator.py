import matplotlib.pyplot as plt
import numpy as np
import re
from alignment import parse_SAM_header

def calculate_mapping_quality_distribution(mapq_dict, mapq, align_type, factor): # BWA MEM mapping quality range [0-60]
    if mapq < factor:
        mapq_dict[align_type]['0-10'] += 1
    elif (mapq >= factor) and (mapq < (factor*2)):
        mapq_dict[align_type]['10-20'] += 1
    elif (mapq >= (factor*2)) and (mapq < (factor*3)):
        mapq_dict[align_type]['20-30'] += 1
    elif (mapq >= (factor*3)) and (mapq < (factor*4)):
        mapq_dict[align_type]['30-40'] += 1
    elif (mapq >= (factor*4)) and (mapq < (factor*5)):
        mapq_dict[align_type]['40-50'] += 1
    else:
        mapq_dict[align_type]['50-60'] += 1 



    
def get_mapq_distribution_bar(bwa_data, bowtie_data, genome_name):
    bwa_mapq_distribution = bwa_data[0]
    bwa_aligned_num = bwa_data[1]
    bwa_misaligned_num = bwa_data[2]
    bwa_unmapped_num = bwa_data[3]
    bowtie_mapq_distribution = bowtie_data[0]
    bowtie_aligned_num = bowtie_data[1]
    bowtie_misaligned_num = bowtie_data[2]
    bowtie_unmapped_num = bowtie_data[3]     
    bwa_aligned = [round(num/bwa_aligned_num*100) for num in bwa_mapq_distribution['aligned'].values()]
    bwa_misaligned = [round(num/bwa_misaligned_num*100) for num in bwa_mapq_distribution['misaligned'].values()]
    bowtie_aligned= [round(num/bowtie_aligned_num*100) for num in bowtie_mapq_distribution['aligned'].values()]
    bowtie_misaligned = [round(num/bowtie_misaligned_num*100) for num in bowtie_mapq_distribution['misaligned'].values()]
    
    
    bottoms_bwa_aligned = []
    bottoms_bwa_misaligned = []
    bottoms_bowtie_aligned = []
    bottoms_bowtie_misaligned = []
    for i in range(len(bwa_aligned)):
        item1=bwa_aligned[i]
        item2=bwa_misaligned[i]
        item3=bowtie_aligned[i]
        item4=bowtie_misaligned[i]
        if item1 < item2:
            bottoms_bwa_misaligned.append(item1)
            bwa_misaligned[i] -= item1
            bottoms_bwa_aligned.append(0)
        elif item1 > item2:
            bottoms_bwa_aligned.append(item2)
            bwa_aligned[i] -= item2
            bottoms_bwa_misaligned.append(0)
        else:
            bottoms_bwa_aligned.append(0)
            bottoms_bwa_misaligned.append(0) 
        if item3 < item4:
            bottoms_bowtie_misaligned.append(item3)
            bowtie_misaligned[i] -= item3
            bottoms_bowtie_aligned.append(0)
        elif item3 > item4:
            bottoms_bowtie_aligned.append(item4)
            bowtie_aligned[i] -= item4
            bottoms_bowtie_misaligned.append(0)
        else:
            bottoms_bowtie_aligned.append(0)
            bottoms_bowtie_misaligned.append(0) 
    
    labels = ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60'] # labels on x axis
    label_loc = np.arange(len(labels)) # positions of labels on x axis, equidistance
    width = 0.2 # width of the bar
    fig, ax = plt.subplots() # creates figure, plot and axis
    bwa_mem_aligned_rects = ax.bar(label_loc - width/2, bwa_aligned, width=0.1, color="green", label = "BWA MEM aligned", bottom=bottoms_bwa_aligned) # draws bars on plot
    bwa_mem_misaligned_rects = ax.bar(label_loc - width/2, bwa_misaligned, width=0.1, color="green", label = "BWA MEM misaligned", alpha=0.3, bottom=bottoms_bwa_misaligned)
    bowtie_aligned_rects = ax.bar(label_loc + width/2, bowtie_aligned, width=0.1, color="orange", label = "Bowtie aligned", bottom=bottoms_bowtie_aligned) # draws bars on plot
    bowtie_misaligned_rects = ax.bar(label_loc + width/2, bowtie_misaligned, width=0.1, color="orange", label = "Bowtie misaligned", alpha=0.3, bottom=bottoms_bowtie_misaligned)
   
    ax.set_xlabel("Mapping quality", fontsize=14)
    ax.set_ylabel("Percentage of reads", fontsize=14)
    ax.set_title("Mapping quality distribution of aligned and misaligned reads of {} genome".format(genome_name), fontsize=14)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    #ax.bar_label(bwa_mem_aligned_rects, padding=3)
    #ax.bar_label(bwa_mem_misaligned_rects, padding=3)
    ax.legend()
    fig.tight_layout()
    plt.show()

def get_percision_bar():
    return
def get_recall_bar(bwa_data, bowtie_data, genome_name):
    
    tools=["BWA MEM", "Bowtie"]
    label_loc = np.arange(1)
    bwa_aligned_num = bwa_data[1]
    bowtie_aligned_num = bowtie_data[1]
    bwa_total = bwa_data[1] + bwa_data[2] + bwa_data[3]
    bowtie_total = bowtie_data[1] + bowtie_data[2] + bowtie_data[3]
    fig, ax = plt.subplots()
    recall_aligned=[round(bwa_aligned_num/bwa_total*100, 2), round(bowtie_aligned_num/bowtie_total*100, 2)] 
    recall_rects = ax.bar(tools, recall_aligned, width=0.5, color="orange")
    ax.set_xlabel("Tools", fontsize=14)
    ax.set_ylabel("Recall of aligned reads (%)", fontsize=14)
    ax.set_title("Recall function of aligned reads, {} genome".format(genome_name), fontsize=14)
    ax.bar_label(recall_rects, padding=3, fontsize=13, color="green")
    fig.tight_layout()
    plt.show()

    return
def get_AUC_func():
    return
def get_fscore_func():
    return        
           

def compare(simulator_sam_file, tool_sam_file, factor):
    
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
                calculate_mapping_quality_distribution(mapq_distribution_dict, tool_sam_alignment_mapq, 'misaligned', factor)
                
            else:
                aligned_reads += 1
                calculate_mapping_quality_distribution(mapq_distribution_dict, tool_sam_alignment_mapq, 'aligned', factor)
                
        
        return(mapq_distribution_dict, aligned_reads, misaligned_reads, unmapped_reads)
            

def compare_(sim_sam_path, bwa_sam_path, bowtie_sam_path, genome_name):
    bwa_data = compare(sim_sam_path, bwa_sam_path, 10)
    bowtie_data = compare(sim_sam_path, bowtie_sam_path, 7)
    
    get_recall_bar(bwa_data, bowtie_data, genome_name)
    # Data visualization, graphical representation
    #get_mapq_distribution_bar(bwa_data, bowtie_data, genome_name)
   
    #get_percision_bar(bwa_data, bowtie_data)
    
    #get_AUC_func()
    #get_fscore_func()

    

simulator_sam_file = "./out_sam/ls_orchid.sam"
bwa_sam_file = "./bwa_sam/ls_orchid_bwa_mem.sam"
bowtie_sam_file = "./bowtie_sam/ls_orchid.sam"
genome_name="ls_orchid"

compare_(simulator_sam_file, bwa_sam_file, bowtie_sam_file, genome_name)