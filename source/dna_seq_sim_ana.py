



from argument_parser import parse_arguments
from simulator import simulate_sequencing
from comparator import compare_single


def main():
    args = parse_arguments()

    if(args.mode == "sim"):
        #print("starting simulation")
        if args.fasta == None:
            print("DnaSeqSimAna: Simulating finished with errors. You need to specify path to the FASTA file.")
        elif args.insert_size < args.read_length:
            print("DnaSeqSimAna: Simulating finished with errors. Insert size can't be smaller than read length!")
        else: simulate_sequencing(args)
    elif(args.mode == "compare"):
        if None in (args.sim, args.bwa_mem, args.bowtie):
            print("DnaSeqSimAna: Comparison finished with errors. Check if you specified required input arguments --sim, --bwa_mem and --bowtie.\nType --help for more info.")
        else:
            #print("starting comparison")
            compare_single(args)


if __name__ == "__main__":
    main()
   
