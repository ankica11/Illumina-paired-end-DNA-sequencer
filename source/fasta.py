from datetime import datetime
import os
from Bio import SeqIO
import timeit



SMALL_FILE_SIZE_TRESHHOLD = 80*2**20  # files smaller that 3MB are considered small


#---------------------------------------parsiranje ulaznog FASTA fajla-------------------------------------------------------------------
def fasta_parsing(fasta_filename):

    t_start=datetime.now()
    #src_dir="./genomes" #directory for genome samples
    #if not os.path.exists(src_dir):
    #    os.system('mkdir ' + src_dir)
    #genome_sample_file_path = os.path.join(src_dir, fasta_filename)
    genome_sample_file_path=fasta_filename
    print("FASTA parsing of file {}, size: {} Bytes has started...".format(fasta_filename, os.path.getsize(genome_sample_file_path)))
    if os.path.getsize(genome_sample_file_path)<=SMALL_FILE_SIZE_TRESHHOLD:
        with open(genome_sample_file_path) as handle:
            genome_sequences_dict = SeqIO.to_dict(SeqIO.parse(handle,"fasta")) #faster but memory dependent cant be used for large files
    else:
        genome_sequences_dict = SeqIO.index(genome_sample_file_path,"fasta")   #slower but can be used for larger files

    t_end=datetime.now()
    print("Finished FASTA file parsing in {}.".format(t_end-t_start))
        
    return genome_sequences_dict



#t1 = timeit.timeit(lambda: fasta_parsing("GCF_000766835.1_Aquila_chrysaetos-1.0.2_cds_from_genomic.fasta"), number=1)
#print("time: {}".format(t1))



def fasta_parsing_extra_large(fasta_filename):
    return


#t1=timeit.timeit(lambda: fasta_parsing("GCF_000766835.1_Aquila_chrysaetos-1.0.2_cds_from_genomic.fasta"), number=1)
#t2=timeit.timeit(lambda: read_genome("GCF_000766835.1_Aquila_chrysaetos-1.0.2_cds_from_genomic.fasta"), number=1)
#print("fasta parsing func time: {}\n read genome func time: {}".format(t1, t2))