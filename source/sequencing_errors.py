import numpy as np

nucleotides = ['A', 'T', 'C', 'G']

def simulate_snv(read, snv_rate):
    return

def simulate_insertion(read, insertion_rate):
    return

def simulate_deletion(read, deletion_rate, additional_nucleotide):
   
    return


def simulate_sequencing_errors(read, snv_rate, insertion_rate, deletion_rate, end_position, ref_sequence_name, genome):
    #print("Applying sequencing errors...")
    read_length=len(read)
    #ref_sequence = str(genome[ref_sequence_name].seq).upper() #fasta parsing with seqio
    #ref_sequence = str(genome[ref_sequence_name]).upper()     #pesadijsko parsiranje faste
    #ref_sequence = genome
    ref_sequence = genome[ref_sequence_name].seq
    del_positions=[]
    ins_positions=[]
    snv_positions=[]
    del_num=0
    ins_num=0
    for err_pos in range(read_length):
        if np.random.random() <= snv_rate:                                            # doesn't change read length
            nucleotides_snv=list(nucleotides)
            if read[err_pos] != 'N':
                nucleotides_snv.remove(read[err_pos])
                new_nucleotide=np.random.choice(nucleotides_snv)
            else: 
                new_nucleotide='N'
            old=read[err_pos]
            snv_positions.append(err_pos)
            #print("Read before: ", read)
            read=read[:err_pos] + new_nucleotide + read[err_pos+1:]
            #print("SNV pos: {}, old: {}, new: {}".format(err_pos,old,new_nucleotide))
            #print("Read:        ", read)

            
        elif np.random.random() <= insertion_rate:                                    # increasses read length but we want to keep it the same 
                                                                                      # so we cut off one base from the end of read 
            inserted_nucleotide=np.random.choice(nucleotides)
            ins_positions.append(err_pos)
            #print("Read before: ", read)
            read=read[:err_pos] + inserted_nucleotide + read[err_pos:-1]
            #print("INS pos: {}, inserted: {}".format(err_pos, inserted_nucleotide))
            #print("Read:        ", read)
            ins_num += 1
                                                                                      
        elif np.random.random() <= deletion_rate:                                     # decreasses read length so its needed to add new additional
            deleted_nucleotide=read[err_pos]
            del_positions.append(err_pos)                                             # base at the end of read to preserve read length the new base is taken from
                                                                                      # reference sequence or if there is no bases left in sequence the new base is 
                                                                                      # chosen randomly                                      
            next_position_ref_seq = end_position + del_num - ins_num
            if next_position_ref_seq < len(ref_sequence):
                additional_nucleotide = ref_sequence[next_position_ref_seq]
               
            else: 
                additional_nucleotide=np.random.choice(nucleotides)
            
            #print("Read before: ", read)
            read=read[:err_pos] + read[err_pos+1:] + additional_nucleotide
            #print("DEL pos: {}, deleted: {}, additional: {}".format(err_pos, deleted_nucleotide, additional_nucleotide))
            #print("Read:        ", read)
            del_num += 1
            

    return (read, snv_positions, ins_positions, del_positions)











