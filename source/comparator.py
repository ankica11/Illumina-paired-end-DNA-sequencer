
import re


def compare(simulator_sam_file, tool_sam_file):
    misaligned_reads = 0
    correctly_aligned_reads = 0
    with open(simulator_sam_file, "r") as simulator_sam_handle, open(tool_sam_file, "r") as tool_sam_handle:
        
        pos = 0
        header_line = tool_sam_handle.readline()
        while header_line.startswith("@"):
            pos = tool_sam_handle.tell()
            header_line = tool_sam_handle.readline()
        tool_sam_handle.seek(pos)


        for sim_sam_line, tool_sam_line in zip(simulator_sam_handle, tool_sam_handle):
            sim_sam_line_splited = re.split(r"\t+", sim_sam_line)
            tool_sam_line_splited = re.split(r"\t+", tool_sam_line)
            if sim_sam_line_splited[3] != tool_sam_line_splited[3]:
                misaligned_reads += 1
            else:
                correctly_aligned_reads += 1
           
            #break
        print("misaligned: ", misaligned_reads)
        print("correctly aligned: ", correctly_aligned_reads)
            
           

simulator_sam_file = "./out_sam/ls_orchid.sam"
tool_sam_file = "./bwa_sam/ls_orchid_cgc.sam"

compare(simulator_sam_file, tool_sam_file)