import numpy as np
from Bio.SeqRecord import SeqRecord

def get_num_of_fragments_2_sequence(coverage, seq_length, paired_read_length):
    """ Function that calculates how many time is required to sequence dna contig into reads
       to satisfy given coverage value. The function uses Lander_Waterman equation to estimate 
       the desired number of reads, and number of inserts (1 insert "produces" 1 pair of reads).
    """
    num = round((coverage*seq_length)/(2*paired_read_length))
    
    return num

def calculate_sizes(single_read_length, ins_size, dna_seq_size): 
    """ Function fixes sizes of single read and insert. The size fix is necessary in 
        a sitation where insert or read sizes are bigger than dna sequence length, 
        because this simulator sequence every dna contig (sequence in genome, chromozomes, genes...)
        separately, it doesn't concatenate all sequences in genome and then generate reads like bwa mem does.

    """
    insert_size = ins_size
    read_length = single_read_length
    if insert_size > dna_seq_size:
        insert_size = dna_seq_size
        if single_read_length > insert_size:
            read_length = insert_size
    return (read_length, insert_size)


def encode_q_score(q_score):
    """ Function receives as an argument an integer value which represents quality score of base 
        and returns its ASCII value accoriding to rule Phred+33 (ASCII representation of quality scores
        is in interval (!, ~)) 
    """
    qs = round(q_score)
    if qs < 0:
        qs = 0
    elif qs > 41:
        qs = 41
    return chr(qs + 33) #Pread+33 encoding (0, 41) scores for illumina 1.8+

def generate_quality_scores(average_quality, read_length, standard_deviation=3):
    if average_quality > 41:
        average_quality = 41
    elif average_quality < 0:
        average_quality = 0

    q_scores_distribution = np.random.normal(average_quality, standard_deviation, read_length)
    q_scores=list(map(encode_q_score, q_scores_distribution))
    
    return ''.join(q_scores)


def generate_q_score_quality_string(q_scores):
    return ''.join(q_scores)

#generate_quality_scores(40,3,40)
#print(generate_q_score_quality_string(generate_quality_scores(40,3,40)))




def generate_read(fastq_file_handle, readRecord, standard_deviation=3):
    pairity = readRecord.annotations["pairity"]
    ord_num = readRecord.annotations["ord_num"]
    seq_num = readRecord.annotations["seq_num"]
    read_id = readRecord.name + ":" + str(ord_num)
    if pairity == 2:
        read_seq = readRecord.seq.reverse_complement()
    else:
        read_seq = readRecord.seq
    readRecord.id = read_id
    quality_scores = generate_quality_scores(readRecord.annotations["average_quality"], len(read_seq), standard_deviation)
    readRecord.letter_annotations["phred_quality"] = quality_scores
    read_fastq_format = "@{}/{}\n{}\n+\n{}\n"
    fastq_file_handle.write(read_fastq_format.format(read_id, pairity, read_seq, quality_scores))
    