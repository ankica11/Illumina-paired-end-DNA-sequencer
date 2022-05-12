from argparse import ArgumentParser

from simulator import simulate_sequencing

#dokle bre vise alo 
def parse_arguments():

    parser = ArgumentParser(prog="dnaSeqSimAna", description="Illumina pair end DNA sequencing simulator and comparator")
    parser.add_argument("--version", "-v", action="version", version="dnaSeqSimAna v1.2022")
    parser.add_argument("--fasta", "-fa", help="Path FASTA file with reference genome", type=str)
    parser.add_argument('--avgquality', help="Average nucletiode quality", type=int, default=40)
    parser.add_argument('--coverage', help="Coverage", type=int, default=2)
    parser.add_argument('--readlength', help="Read length", type=int, default=150)
    parser.add_argument('--insertsize', help="Insert size", type=int, default=350)
    parser.add_argument('--snv_rate', help="Probability snv mutation", type=float, default=0)
    parser.add_argument('--ins_rate', help="Probability insert mutation", type=float, default=0)
    parser.add_argument('--del_rate', help="Probability delete mutation", type=float, default=0)
    

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    simulate_sequencing(args)