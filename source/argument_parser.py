from argparse import ArgumentParser

from simulator import simulate_sequencing

#dokle bre vise alo 
def parse_arguments():

    parser = ArgumentParser(prog="dnaSeqSimAna", description="Illumina pair end DNA sequencing simulator and comparator")
    parser.add_argument("--version", "-v", action="version", version="dnaSeqSimAna v1.2022")
    parser.add_argument("--fasta", "-fa", help="FASTA file path with reference genome", type=str)
    parser.add_argument('--avg_quality', help="Average nucletiode quality Illumina 1.8 scale [0,40]", type=int, default=40)
    parser.add_argument('--coverage', help="Coverage", type=int, default=2)
    parser.add_argument('--read_length', help="Read length", type=int, default=150)
    parser.add_argument('--insert_size', help="Insert size", type=int, default=300)
    parser.add_argument('--snv_rate', help="SNV error rate", type=float, default=0.005)
    parser.add_argument('--ins_rate', help="Insertion error rate", type=float, default=0.005)
    parser.add_argument('--del_rate', help="Deleteion error rate", type=float, default=0.005)
    parser.add_argument('--out_sam', help="Directory path for output SAM files", type=str, default="./out_sam")
    parser.add_argument('--out_fastq', help="Directory path for output FASTQ files", type=str, default="./out_fastq")
    parser.add_argument('--bwa_mem', help="BWA MEM SAM file path", type=str)
    parser.add_argument('--bowtie', help="Bowtie SAM file path", type=str)
    parser.add_argument('--sim', help="Simulator SAM file path", type=str)
    parser.add_argument('--mode', help="--mode compare for comparison mode, --mode sim for simulation mode. It's necssary to start simulation first and then comparison.", type=str, default="sim")
    parser.add_argument('--gen_name', help="Genome name", type=str)

    return parser.parse_args()

