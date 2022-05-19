from datetime import datetime
from email.errors import FirstHeaderLineIsContinuationDefect
import timeit
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
import os

from fasta import fasta_parsing
from sequencing_errors import simulate_sequencing_errors
from read import calculate_sizes, generate_read
from read import get_num_of_fragments_2_sequence
from alignment import generate_alignment_and_write_2_SAM
from alignment import write_SAM_header



genome_sequences_dict = {}          # for smaller fasta files, the dictionary is kept in RAM so access is faster
genome_sequences_index = {}         # for larger fasta files, the file is kept on disk and is loaded into RAM when the data is needed, slower access 



        
def multi2single_sequence_fasta():
    for seq_name, seq in genome_sequences_dict.items():
        fasta_single = open(seq_name+".fa","w")
        fasta_single.write(str(seq.seq))



#------------------------------------FASTQ and SAM generating----------------------------------------------------------------------------------------------
def generate_fastq_files(coverage, single_read_length, average_quality, ins_size, fasta_filename, snv_rate, ins_rate, del_rate,\
                         genome_sequences_dict, fastq_path="out_fastq", sam_path="out_sam"):
    """ Function that loops through the reference genome which is presented as a dictionary {sequence_name, SeqRecord object}
        and does the sequencing process on each sequence/contig in genome separately as a result of this process pairs of reads
        are generated and written into _1.fastq and _2.fastq files. The function also generates the SAM file with aligned reads
        to reference genome and other info about alignments itself. 
        Arguments: 
            coverage - the average number of reads that need to "cover" a single nucleotide in dna sequence during the aligning process
            paired_read_length - required read length
            average quality - required average quality of nucleotides, the greater the quality the lower the probability for the base to be called wrongly
            ins_size - the length of dna fragment to be sequenced into pair of reads in opposite directions 
            fasta_filename - name of fasta file with stored reference genome
            snv_rate - single nucleotide variations rate, the probability that during the sequencing process one base is substituted with another
            ins_rate - insertion rate, probability that read contain additional bases that don't exist in reference genome
            del_rate - deletion rate, probability that some bases are missing in reads compared to reference genome
    
    """
    t_start=datetime.now()
    
    out_dir_fastq_path = fastq_path      # directory for output fastq files of simulator, if not defined by user default is out_fastq
    
    if not os.path.exists(out_dir_fastq_path):
        os.system('mkdir ' + out_dir_fastq_path)
    #creating empty fastq files in out directory
    fastq_filename_1=fasta_filename + "_1.fastq"
    fastq_filename_2=fasta_filename + "_2.fastq"
    fastq_files = [os.path.join(out_dir_fastq_path, fastq_filename_1),\
                   os.path.join(out_dir_fastq_path, fastq_filename_2)]

    #removing existing files
    for file in fastq_files:
        if os.path.exists(file):
            os.remove(file)
    (fastq1_handle, fastq2_handle) = [open(filename, "w") for filename in fastq_files]

    #creating empty SAM file
    out_sam_dir_path = sam_path      # directory for output sam files of simulator, if not defined by user default is out_sam
    if not os.path.exists(out_sam_dir_path):
        os.system('mkdir ' + out_sam_dir_path)
    sam_filename = fasta_filename + ".sam"
    sam_file_path = os.path.join(out_sam_dir_path, sam_filename)
    if os.path.exists(sam_file_path):
        os.remove(sam_file_path)
    sam_handle = open(sam_file_path, "w")
    write_SAM_header(genome_sequences_dict, sam_handle, fasta_filename, fastq_filename_1, fastq_filename_2)
    # looping through the genome sequences and for each sequence choosing fragments positions (insert positions) to sequence into 2-paired reads
    # genome sequence dict is a dictionary (sequence name, SeqRecord object)
    # SeqRecord is a class from Biopython library and containes fields like SEQUENCE NAME, SEQUENCE ID, SEQUENCE(itself), plus additonal fields
    for seq_num, rec_seq in enumerate(genome_sequences_dict.values(), start=1):
        seq_name = rec_seq.id
        dna_seq = rec_seq.seq.upper()
        #if len(dna_seq) < ins_size:     # maybe we shouldn't even process sequences smaller than insert size????
        #    continue
        (paired_read_length, insert_size) = calculate_sizes(single_read_length, ins_size, len(dna_seq))
        num_of_fragments = get_num_of_fragments_2_sequence(coverage, len(dna_seq), paired_read_length) # calculating number of fragments that will satisfy the required coverage (how many reads are needed to cover each base in dna sequence) 
       
        generated_read_pairs_num = 0 # how many pairs of reads are generated

        for i in range(num_of_fragments):
            if(len(dna_seq)>insert_size): 
                pos = np.random.randint(0, len(dna_seq) - insert_size)
            else:
                pos = 0
            fragment=dna_seq[pos:(pos+insert_size)] # extraxting dna fragment from dna sequence, each fragment is sequenced into pair of reads
            generated_read_pairs_num += 1
            
            
            # -----------------------Creating Read 1 with added sequencing errors-----------------------------------------------------------------------
            read1_clean = fragment[0:paired_read_length] # Read 1 of pair is sequenced in regular direction (3'->5')
            #print("Read1 before: {}".format(read1_clean))
            (read1_mutated, snv_pos1, ins_pos1, del_pos1) = simulate_sequencing_errors(read1_clean, snv_rate, ins_rate, del_rate,\
                                                                                       pos+paired_read_length, seq_name, genome_sequences_dict)
            read1 = SeqRecord(read1_mutated, id="", name=seq_name,\
                              annotations={'pairity' : 1, 'average_quality' : average_quality, 'ord_num' : generated_read_pairs_num, 'seq_num':seq_num})
            generate_read(fastq1_handle, read1)
            #print("Read1 after:  {}".format(read1.seq))


            # ------------------------Creating Read 2 with added sequencing errors----------------------------------------------------------------------
            read2_clean = fragment[(insert_size-paired_read_length):insert_size] # Read 2 of pair is sequenced in reverse direction (5'->3')
            #print("Read2 before: {}".format(read2_clean))
            (read2_mutated, snv_pos2, ins_pos2, del_pos2) = simulate_sequencing_errors(read2_clean, snv_rate, ins_rate, del_rate,\
                                                                                       pos+insert_size, seq_name, genome_sequences_dict)
            read2 = SeqRecord(read2_mutated, id="", name=seq_name,\
                              annotations={'pairity' : 2, 'average_quality' : average_quality, 'ord_num' : generated_read_pairs_num, 'seq_num':seq_num})
            generate_read(fastq2_handle, read2)
            #print("Read2 after:  {}".format(read2.seq))
            
            
            
            #--------------------------Generating alignments for both reads in the same pair and writing them to SAM file----------------------
            generate_alignment_and_write_2_SAM(sam_handle, read1, pos + 1, pos+insert_size-paired_read_length+1, insert_size,\
                                               snv_pos1, ins_pos1, del_pos1)
            generate_alignment_and_write_2_SAM(sam_handle, read2, pos+insert_size-paired_read_length+1, pos+1, -insert_size,\
                                               snv_pos2, ins_pos2, del_pos2)
            
            

    fastq1_handle.close()
    fastq2_handle.close()
    sam_handle.close()

    t_end = datetime.now()
    print("Generated 2x{} Bytes of FASTQ files at {} directory.\nGenerated {} Bytes of SAM file at {} directory.".format(os.path.getsize(fastq_files[0]), out_dir_fastq_path,\
                                                  os.path.getsize(sam_file_path), out_sam_dir_path))

nucleotides = ['A', 'T', 'C', 'G']




def simulate_sequencing(args):

    t_start = datetime.now()
    print("DnaSeqSimAna: Illumina paired-end sequencing simulation has started at {}...".format(t_start))

    global genome_sequences_index
    global genome_sequences_dict
    fasta_file_path = args.fasta
    coverage = args.coverage
    average_quality = args.avg_quality
    insert_size = args.insert_size
    paired_read_length = args.read_length
    snv_rate = args.snv_rate
    ins_rate = args.ins_rate
    del_rate = args.del_rate
    fastq_path = args.out_fastq
    sam_path = args.out_sam
    
    
    genome_sequences_dict = fasta_parsing(fasta_file_path) # time: 1.123226 0.956
    fasta_filename = os.path.split(fasta_file_path)[1]
    generate_fastq_files(coverage, paired_read_length, average_quality, insert_size, fasta_filename.replace(".fasta",""), snv_rate,\
        ins_rate, del_rate, genome_sequences_dict, fastq_path, sam_path)

    t_end = datetime.now()
    print("DnaSeqSimAna: Finished sequencing in {}.".format(t_end-t_start))
    
    
#t1 = timeit.timeit(lambda: simulate(), number=1)
#print(t1)    

#simulate_sequencing()   

