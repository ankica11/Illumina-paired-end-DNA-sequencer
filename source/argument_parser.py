from argparse import ArgumentParser

from simulator import simulate_sequencing


def parse_arguments():

    parser = ArgumentParser(prog="dnaSeqSimAna", description="Illumina pair end DNA sequencing simulator and comparator")
    parser.add_argument("--version", "-v", action="version", version="dnaSeqSimAna v1.2022")
    parser.add_argument("--fasta", "-fa", help="FASTA file path with reference genome", type=str)
    parser.add_argument('--avg_quality', help="Average nucletiode quality Illumina 1.8 scale [0,41]", type=int, default=40)
    parser.add_argument('--coverage', help="Coverage", type=float, default=3)
    parser.add_argument('--read_length', help="Read length", type=int, default=150)
    parser.add_argument('--insert_size', help="Insert size", type=int, default=400)
    parser.add_argument('--snv_rate', help="SNV error rate", type=float, default=0.005)
    parser.add_argument('--ins_rate', help="Insertion error rate", type=float, default=0.005)
    parser.add_argument('--del_rate', help="Deleteion error rate", type=float, default=0.005)
    parser.add_argument('--out_sam', help="Directory path for output SAM files", type=str, default="./out_sam")
    parser.add_argument('--out_fastq', help="Directory path for output FASTQ files", type=str, default="./out_fastq")
    parser.add_argument('--bwa_mem', help="BWA MEM SAM file path", type=str)
    parser.add_argument('--bowtie', help="Bowtie SAM file path", type=str)
    parser.add_argument('--sim', help="Simulator SAM file path", type=str)
    parser.add_argument('--mode', help="--mode compare for comparison mode, --mode sim for simulation mode, --mode multi for multiple comparison. It's necssary to start simulation first and then comparison.", type=str, default="sim")
    parser.add_argument('--gen_name', help="Genome name", type=str)
    parser.add_argument('--sam_list', help="List of paths with SAM files to compare", nargs="+")
    parser.add_argument('--num_of_sets', help="Number of sets for each genome", type=int)
    # for multiple comparison its neccessary that bowtie SAM files have _bowtie in filename, and bwa mem SAM files _bwa_mem
    # also SAM files generated from same fastq files, from same set need to have same id
    # every SAM group of files needs to have id as ord number starting from 0 and ending with number of sets - 1
    # example of command: ./source/dna_seq_sim_ana.py --mode multi --sam list path/to/sam/files/genome1_name 
    # path/to/sam/files/genome2_name path/to/sam/files/genome3_name
    return parser.parse_args()

