from argparse import ArgumentParser

from source.simulator import simulate_sequencing


def parse_arguments():

    parser = ArgumentParser(prog="dnaSeqSimAna", description="Illumina pair end DNA sequencing simulator and comparator")
    parser.add_argument("--version", "-v", action="version", version="dnaSeqSimAna v1.2022")
    parser.add_argument("--fasta", "-fa", help="Path FASTA file with reference genome", type=str)
    parser.add_argument('--avgquality', help="Average nucletiode quality", type=int, default=2)
    parser.add_argument('--coverage', help="Coverage", type=int, default=2)
    parser.add_argument('--readlength', help="Read length", type=int, default=30)
    parser.add_argument('--insertsize', help="Insert size", type=int, default=80)
    parser.add_argument('--probsnv', help="Probability snv mutation", type=float, default=0)
    parser.add_argument('--probins', help="Probability insert mutation", type=float, default=0)
    parser.add_argument('--probdel', help="Probability delete mutation", type=float, default=0)
    

    return parser.parse_args()


args = parse_arguments()

simulate_sequencing(args.fasta, args.coverage, args.avgquality, args.insertsize, args.readlength, args.probsnv, args.probins, args.probdel)