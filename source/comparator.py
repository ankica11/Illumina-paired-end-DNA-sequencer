
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
            
        print("misaligned: ", misaligned_reads)
        print("correctly aligned: ", aligned_reads)
        print("Unmapped reads: ", unmapped_reads)
        print(mapq_distribution_dict.items())
            
           

simulator_sam_file = "./out_sam/ls_orchid.sam"
tool_sam_file = "./bwa_sam/ls_orchid_cgc.sam"

compare(simulator_sam_file, tool_sam_file)