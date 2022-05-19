import sys
import traceback

import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from matplotlib import lines
from matplotlib.pyplot import close
import numpy
from regex import D
from fasta import fasta_parsing
from read import calculate_sizes
from read import generate_read
from alignment import generate_cigar_string
from alignment import generate_alignment_and_write_2_SAM
from sequencing_errors import simulate_sequencing_errors
from comparator import compare

# Testing simulator's functionalities

# fasta.py tests

def fasta_parsing_test():
    test_genome = "./genomes/sample.fasta"

    genome_dict = fasta_parsing(test_genome)
    genome_dict_test = {'NC_001612.1' : 'ATTAAACAGCCTGTGGGTTGTACCCACCCACAGGGCCCACTGGGCGCTAGCACTCTGATTCTACGGAATCCTTGTGCGCCATTACGAATGCTATCATAACTATTGCAGACTCTATGCCGATGACTCCACGACGAGTCTTTCTAGTTTTCAGTTCGTGTTT',\
                        'NC_001613.1' : 'TTAAAACAGCCTGTGGGTTGTACCCACCCACAGGGCCCACTGGGCGCTAGCACTCTGATTCTACGGAATCCTTGTGCGCC'}
    try:
        for (sequence_name, sequence_record), (sequence_name_test, sequence_test) in zip(genome_dict.items(), genome_dict_test.items()):
            assert sequence_name == sequence_name_test and sequence_record.seq == sequence_test
        print('-------Assertion test for function "fasta_parsing" successfully finished!-------')
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1)



# read.py tests

def calculate_sizes_test():
    try:
        sequence_length = 100
        insert_size = 150
        read_length = 120
        assert calculate_sizes(read_length, insert_size, sequence_length) == (100, 100)
        sequence_length = 100
        insert_size = 150
        read_length = 80
        assert calculate_sizes(read_length, insert_size, sequence_length) == (80, 100)
        print('-------Assertion test for function "calculate_sizes" successfully finished!-------')
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1)

def generate_read_test():
    try: 
        read_1 = SeqRecord(Seq('AAATCGTA'), id="", name="sequence_name",\
                            annotations={'pairity' : 1, 'average_quality' : 10, 'ord_num' : 1, 'seq_num':1})
        with open("./out_fastq/sample_test_1.fastq", "w") as test_file_handle:
            generate_read(test_file_handle, read_1)
        with open("./out_fastq/sample_test_1.fastq", "r") as test_file_handle:
            lines = test_file_handle.read().split()
            assert lines[0] == "@sequence_name:1/1"
            assert lines[1] == "AAATCGTA"
            assert lines[2] == "+"
            assert len(lines[3]) == 8
            for c in lines[3]:
                assert (c >= '!') and (c <= 'J')
        read_1.annotations['pairity'] = 2
        with open("./out_fastq/sample_test_1.fastq", "w") as test_file_handle:
            generate_read(test_file_handle, read_1)
        with open("./out_fastq/sample_test_1.fastq", "r") as test_file_handle:
            lines = test_file_handle.read().split()
            assert lines[1] == "TACGATTT"
        print('-------Assertion test for function "generate_read" successfully finished!-------')
        
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1)


         
# alignment.py tests

def generate_cigar_string_test():
    try:
        snv_pos = [0, 3, 9]
        ins_pos = [2, 6, 13]
        del_pos = [11, 15, 16]
        cigar_test = "1S1M1I1S2M1I2M1S1M1D1M1I1M2D3M"
        read = "AAAATTTTCCCCGGGGATCG"
        assert generate_cigar_string(read, snv_pos, ins_pos, del_pos) == cigar_test
        print('-------Assertion test for function "generate_read" successfully finished!-------')
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1)


def generate_alignement_and_write_2_SAM_test():
    snv_pos = [0, 3]
    ins_pos = [2, 6]
    del_pos = [1]
    cigar_test = "1S1D1I1S2M1I3M"
    read_1 = SeqRecord(Seq('AAATCGTAGG'), id="", name="sequence_name",\
                            annotations={'pairity' : 2, 'average_quality' : 40, 'ord_num' : 1, 'seq_num':1})
    try:
        with open("./out_fastq/sample_test_1.fastq", "w") as fastq, open("./out_sam/sample_test.sam", "w") as sam:
            generate_read(fastq, read_1)
            generate_alignment_and_write_2_SAM(sam, read_1, 10, 50, 70, snv_pos, ins_pos, del_pos)
        with open("./out_fastq/sample_test_1.fastq", "r") as fastq, open("./out_sam/sample_test.sam", "r") as sam:
            fields = sam.read().split()
            lines = fastq.read().split()
            assert fields[10] == lines[3][::-1]
            assert fields[5] == cigar_test
        print('-------Assertion test for function "generate_alignement_and_write_2_SAM" successfully finished!-------')
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1) 



# sequencing_errors.py tests

def simulate_sequencing_errors_test():
    try:
        ref_seq = SeqRecord(Seq("CCCCTTATGGGTGTTTTACGTAGAAATACTATTTTTGCTTATTTTGACGATCCCCGATACAGAAAAGATAAAAAGGGTTCAGGAATTGTTAAATTTAGATATAGGACCCTAGAAGAAGAATATAGGACTCGAGCGGAAGACTCAGAAGAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), "seq_name", "seq_name")
        read1 = Seq("CCCCTTATGGGTGTTTTACGTAGAAATACTATTTTTGCTTATTTTGACGATCCCCGATACAGAAAAGATAAAAAGGGTTCAGGAATTGTTAAATTTAGATATAGGACCCTAGAAGAAGAATATAGGACTCGAGCGGAAGACTCAGAAGAG")
        genome = {'seq_name' : ref_seq}
        snv_rate = 0.007
        ins_rate = 0.005
        del_rate = 0.003
        snv_num = 0
        ins_num = 0
        del_num = 0
        num_of_reads = 1000
        for i in range(0, num_of_reads):
                _, snv_pos, ins_pos, del_pos = simulate_sequencing_errors(read1, snv_rate, ins_rate, del_rate, 16, "seq_name", genome)
                snv_num += len(snv_pos)
                ins_num += len(ins_pos)
                del_num += len(del_pos)
        #print(snv_num/(500*len(read1)))
        #print(ins_num/(500*len(read1)))
        #print(del_num/(500*len(read1)))
        assert  ((snv_rate - 1) <= snv_num/(num_of_reads*len(read1))) and (snv_num/(num_of_reads*len(read1)) <= (snv_rate + 1))
        assert  ((ins_rate - 1) <= ins_num/(num_of_reads*len(read1))) and (ins_num/(num_of_reads*len(read1)) <= (ins_rate + 1))
        assert  ((del_rate - 1) <= del_num/(num_of_reads*len(read1))) and (del_num/(num_of_reads*len(read1)) <= (del_rate + 1))
        print('-------Assertion test for function "simulate_sequencing_errors" successfully finished!-------')
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1) 


# comparator.py tests
def compare_test():
    try:
        sim_path = "./out_sam/test_sim.sam"
        tool_path = "./out_sam/test_tool.sam"
        mapq_dict, aln_num, mis_num, unmapped_num = compare(sim_path, tool_path, 10)
        assert aln_num == 9
        assert mis_num == 6
        assert unmapped_num == 3
        assert mapq_dict['aligned']['40-50'] == 1
        assert mapq_dict['aligned']['50-60'] == 8
        assert mapq_dict['misaligned']['10-20'] == 6
        assert sum(mapq_dict["aligned"].values()) == aln_num
        assert sum(mapq_dict["misaligned"].values()) == mis_num
        print('-------Assertion test for function "compare" successfully finished!-------')
    except AssertionError:
        type, value, tb = sys.exc_info()
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]

        print('Err: Assertion test failed for function "{}" in {}\nLine: {}, statement: {}.'.format(func, filename, line, text))
        exit(1)  




def start_tests():
    print("DNASeqSimAna testing mode has started...")
    print("Testing fasta.py module...")
    fasta_parsing_test()
    print("Testing read.py module...")
    calculate_sizes_test()
    generate_read_test()
    print("Testing alignment.py module...")
    generate_cigar_string_test()
    generate_alignement_and_write_2_SAM_test()
    print("Testing sequencing_errors.py module...")
    simulate_sequencing_errors_test()
    print("Testing comparator.py module...")
    compare_test()
    print("Finished testing!")


start_tests()

    
 


           







