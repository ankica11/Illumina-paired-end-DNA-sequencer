from cProfile import label
from datetime import datetime
from io import BufferedRWPair
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
    if bwa_aligned_num == 0:
        bwa_aligned_num = 100
    if bowtie_aligned_num == 0:
        bowtie_aligned_num = 100
    if bwa_misaligned_num == 0:
        bwa_misaligned_num = 100
    if bowtie_misaligned_num == 0:
        bowtie_misaligned_num = 100    
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
    width = 0.4 # width of the bar
    fig, ax = plt.subplots() # creates figure, plot and axis
    bwa_mem_aligned_rects = ax.bar(label_loc - width/2, bwa_aligned, width=0.3, color="#F24BC0", label = "BWA MEM aligned", bottom=bottoms_bwa_aligned) # draws bars on plot
    bwa_mem_misaligned_rects = ax.bar(label_loc - width/2, bwa_misaligned, width=0.3, color="#F24BC0", label = "BWA MEM misaligned", alpha=0.3, bottom=bottoms_bwa_misaligned)
    bowtie_aligned_rects = ax.bar(label_loc + width/2, bowtie_aligned, width=0.3, color="#3666AC", label = "Bowtie aligned", bottom=bottoms_bowtie_aligned) # draws bars on plot
    bowtie_misaligned_rects = ax.bar(label_loc + width/2, bowtie_misaligned, width=0.3, color="#3666AC", label = "Bowtie misaligned", alpha=0.3, bottom=bottoms_bowtie_misaligned)
   
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


def get_accuracy_histogram_single_mode(bwa_data, bowtie_data, genome_name):
    
    tools = ["BWA MEM", "Bowtie"]
   

    bwa_aligned_num = bwa_data[1]
    bwa_misaligned_num = bwa_data[2]
    bwa_unmapped_num = bwa_data[3]
    bwa_total = bwa_data[1] + bwa_data[2] + bwa_data[3]
    bowtie_aligned_num = bowtie_data[1]
    bowtie_misaligned_num = bowtie_data[2]
    bowtie_unmapped_num = bowtie_data[3]
    bowtie_total = bowtie_data[1] + bowtie_data[2] + bowtie_data[3]
    bwa_mem_statistics = [round(bwa_aligned_num/bwa_total*100, 2), round(bwa_misaligned_num/bwa_total*100, 2), round(bwa_unmapped_num/bwa_total*100, 2) ]
    bowtie_statistics = [round(bowtie_aligned_num/bwa_total*100, 2), round(bowtie_misaligned_num/bwa_total*100, 2), round(bowtie_unmapped_num/bwa_total*100, 2) ] 
   
    width = 0.3
    labels = ['aligned', 'misaligned', 'unmapped']
    label_loc = np.arange(len(labels))
    fig, ax = plt.subplots()
    bwa_mem_rects = ax.bar(label_loc - width/2, bwa_mem_statistics, width=0.2, color="#F24BC0", label="BWA MEM")
    bowtie_rects = ax.bar(label_loc + width/2, bowtie_statistics, width=0.2, color="#3666AC", label="Bowtie")
    ax.set_ylabel("Accuracy (%)", fontsize=14)
    ax.set_title("Percentage of aligned, misaligned and unmapped reads for {} genome".format(genome_name), fontsize=15)
    ax.set_xticks(label_loc)
    ax.set_xticklabels(labels)
    ax.bar_label(bwa_mem_rects, padding=3, fontsize=13)
    ax.bar_label(bowtie_rects, padding=3, fontsize=13)
    ax.legend()
    fig.tight_layout()
    plt.show()




def get_precision_histogram_single_mode(bwa_data, bowtie_data, genome_name):
    bwa_true_positives = bwa_data[1]
    bwa_false_positives = bwa_data[2]
    
    bowtie_true_positives = bowtie_data[1]
    bowtie_false_positives = bowtie_data[2]

    precision_statistics = [round(bwa_true_positives/(bwa_true_positives+bwa_false_positives)*100, 2),\
                            round(bowtie_true_positives/(bowtie_true_positives+bowtie_false_positives)*100, 2)]
    
    width = 0.03
    labels = ['BWA MEM', 'Bowtie']
    fig, ax = plt.subplots()
    precision_rects = ax.bar(labels, precision_statistics, label="precision", width=0.3)
    ax.set_ylabel("Precision (%)", fontsize=14)
    ax.set_title("Precision of aligners for {} genome".format(genome_name), fontsize=15)

    ax.bar_label(precision_rects, padding=3, fontsize=13)
   
   
    fig.tight_layout()
    plt.show()



     
# recall = sensitivity = TP/(TP+FN)   
def get_recall_histogram_single(bwa_data, bowtie_data, genome_name):
    # True positives are correctly aligned reads, false negatives are unmapped reads
    bwa_true_positives = bwa_data[1]
    bwa_false_negatives = bwa_data[3] 
    
    bowtie_true_positives = bowtie_data[1]
    bowtie_false_negatives = bowtie_data[3]

    recall_statistics = [round(bwa_true_positives/(bwa_true_positives+bwa_false_negatives)*100, 2),\
                            round(bowtie_true_positives/(bowtie_true_positives+bowtie_false_negatives)*100, 2)]
    
    width = 0.03
    labels = ['BWA MEM', 'Bowtie']
    fig, ax = plt.subplots()
    precision_rects = ax.bar(labels, recall_statistics, label="recall", width=0.3)
    ax.set_ylabel("Recall (%)", fontsize=14)
    ax.set_title("Recall(sensitivity) of aligners for {} genome".format(genome_name), fontsize=15)

    ax.bar_label(precision_rects, padding=3, fontsize=13)
   
   
    fig.tight_layout()
    plt.show()
    



def get_fcore_histogram_single_mode(bwa_data, bowtie_data, genome_name):
    bwa_true_positives = bwa_data[1]
    bwa_false_negatives = bwa_data[3] 
    bwa_false_positives = bwa_data[2]

    precision_bwa = bwa_true_positives / (bwa_true_positives + bwa_false_positives)
    recall_bwa = bwa_true_positives / (bwa_true_positives + bwa_false_negatives)
    fscore_bwa = 2*precision_bwa*recall_bwa/(precision_bwa + recall_bwa)
    
    bowtie_true_positives = bowtie_data[1]
    bowtie_false_negatives = bowtie_data[3]
    bowtie_false_positives = bowtie_data[2]

    precision_bowtie = bowtie_true_positives / (bowtie_true_positives + bowtie_false_positives)
    recall_bowtie = bowtie_true_positives / (bowtie_true_positives + bowtie_false_negatives)
    fscore_bowtie = 2*precision_bowtie*recall_bowtie/(precision_bowtie + recall_bowtie)

    fscore_statistics = [round(fscore_bwa, 4), round(fscore_bowtie, 4)]

    width = 0.03
    labels = ['BWA MEM', 'Bowtie']
    fig, ax = plt.subplots()
    precision_rects = ax.bar(labels, fscore_statistics, label="fscore", width=0.3)
    ax.set_ylabel("F_score (%)", fontsize=14)
    ax.set_title("F_score of aligners for {} genome".format(genome_name), fontsize=15)

    ax.bar_label(precision_rects, padding=3, fontsize=13)
   
   
    fig.tight_layout()
    plt.show()

    
    
    return 



def compare(simulator_sam_file, tool_sam_file, factor):
    """ Takes file paths of simulator SAM file and tool SAM file as arguments,
        calculates mapping quality distribution of correctly aligned and misaligned reads,
        number of correctly aligned, misaligned and unmapped reads. Mapping qualities are
        divided into 6 categories: 0-10, 10-20, 20-30, 30-40, 40-50, 50-60; since BWA MEM mapq 
        range is [0,60] and Bowtie mapq range is [0, 42] so Bowtie mapq values are linearized 
        in order to fit the BWA MEM mapping quality scale. Read is considered as correctly aligned
        only if its mapping position in tool's SAM file is equal to its position in simulator's SAM file.
        If that's not the case and read is not unmapped, then read is considered as misaligned.
        Alignments with 0x0004 flag set are unmapped.
        In terms of classification correctly aligned reads are considered as TRUE POSITIVES as both simulator and
        tool have mapped them the same way. Wrongly mapped reads or misaligned reads are considered as FALSE POSITIVES because 
        tool has indeed mapped them but not at the right position (like simulator) and they can later lead to false positives variant calls.
        Unmapped reads are FALSE NEGATIVES since tool hasn't aligned them at all, nor correctly nor wrongly. Unmapped 
        reads are not as critical as misaligned since they are easily filtered out and don't bring confusion in later stages of
        DNA analysis.
    """
    
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
            

def compare_single(args):
    t_start = datetime.now()
    print("DnaSeqSimAna: Comparison of BWA MEM and Bowtie has started at {}...".format(t_start))

    sim_sam_path = args.sim
    bwa_sam_path = args.bwa_mem
    bowtie_sam_path = args.bowtie
    genome_name = args.gen_name
    bwa_data = compare(sim_sam_path, bwa_sam_path, 10)
    bowtie_data = compare(sim_sam_path, bowtie_sam_path, 7)


    # Data visualization, graphical representation single mode
    get_mapq_distribution_bar(bwa_data, bowtie_data, genome_name)
    get_accuracy_histogram_single_mode(bwa_data, bowtie_data, genome_name)
    get_precision_histogram_single_mode(bwa_data, bowtie_data, genome_name)
    get_recall_histogram_single(bwa_data, bowtie_data, genome_name)
    get_fcore_histogram_single_mode(bwa_data, bowtie_data, genome_name)

    t_end = datetime.now()
    print("DnaSeqSimAna: Comparison finished in {}.".format(t_end))
   
    


#simulator_sam_file = "./testing/Staphylococcus borealis/3/GCF_003042555.1_ASM304255v1_genomic.sam"
#bwa_sam_file = "./testing/Staphylococcus borealis/3/GCF_003042555.1_ASM304255v1_genomic_bwa_mem.sam"
#bowtie_sam_file = "./testing/Staphylococcus borealis/3/GCF_003042555.1_ASM304255v1_genomic_bowtie.sam"
#genome_name="Staphylococcus borealis"

#compare_single(simulator_sam_file, bwa_sam_file, bowtie_sam_file, genome_name)

