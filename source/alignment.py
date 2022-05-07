def generate_alignment_and_write_2_SAM(sam_file_handle, read_record, aln_pos, mate_aln_pos, mate_distance,\
                                       snv_pos, ins_pos, del_pos):
     sam_file_format = "{}\tFLAG\t{}\t{}\tMAPQ\t{}\t=\t{}\t{}\t{}\t{}\tTAG\tVTYPE\tVALUE\n"
     read_id = read_record.id
     ord_num = read_record.annotations["ord_num"]
     seq_num = read_record.annotations["seq_num"]
     alignment_id = "DNA-SEQ-SIM-ANA" + ":" + str(seq_num) + ":" + str(ord_num)
     seq_name = read_record.name
     read = read_record.seq
     if read_record.annotations["pairity"] == 2:
        quality_scores_string = read_record.letter_annotations["phred_quality"][::-1] # q-scores for read 2 in SAM file shpuld'nt be reversed, so as read sequence 
     else: 
        quality_scores_string = read_record.letter_annotations["phred_quality"]
     cigar = generate_cigar_string(read, snv_pos, ins_pos, del_pos)
     
     sam_file_handle.write(sam_file_format.format(alignment_id, seq_name, aln_pos, cigar, mate_aln_pos, mate_distance, read, quality_scores_string))
    

def generate_cigar_string(read, snv_positions, insertions_positions, deletions_positions):
    cigar_codes = {'match' : 'M', 'snv' : 'S', 'deletion' : 'D', 'insertion' : 'I'}
    cigar = []
    curr_del = 0
    curr_ins = 0
    curr_snv = 0
    curr_count = 0
    curr_code = cigar_codes['match']
    for read_pos in range(len(read)):
        if curr_snv < len(snv_positions) and read_pos == snv_positions[curr_snv]:
            if curr_code == cigar_codes['snv']:
                curr_count += 1
            else:
                cigar.append(str(curr_count)+curr_code)
                curr_count = 1
                curr_code = cigar_codes['snv']
            curr_snv += 1
        elif curr_ins < len(insertions_positions) and read_pos == insertions_positions[curr_ins]:
            if curr_code == cigar_codes['insertion']:
                curr_count += 1
            else:
                cigar.append(str(curr_count) + curr_code)
                curr_count = 1
                curr_code = cigar_codes['insertion']
            curr_ins += 1
        elif curr_del < len(deletions_positions) and read_pos == deletions_positions[curr_del]:
            if curr_code == cigar_codes['deletion']:
                curr_count += 1
            else:
                cigar.append(str(curr_count) + curr_code)
                curr_count = 1
                curr_code = cigar_codes['deletion']
            curr_del += 1
        else:
            if curr_code == cigar_codes['match']:
                curr_count += 1
            else:
                cigar.append(str(curr_count)+curr_code)
                curr_count=1
                curr_code = cigar_codes['match']
    cigar.append(str(curr_count)+curr_code)
    if cigar[0][0] == str(0):
        cigar=cigar[1:]
    return ''.join(cigar)   