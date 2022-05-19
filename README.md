# Illumina-paired-end-DNA-sequencer
A Python program that simulates DNA sequencing of paired-end Illumina reads and their aligning to the reference genome. It's useful for comparison of different aligners like BWA-MEM or Bowtie. 

# Usage
## Simulating mode(default)
Example of usage: ./source/dna_seq_sim_ana.py --mode sim --fasta "C:/Users/user/path/to/FASTA/file.fasta" --insert_size 400 --read_lenght 150 --coverage 3 --avg_quality 40 --snv_rate 0.009 
--ins_rate 0.009 --del_rate 0.009 --out_sam "D:/dir/for/SAM/files" --out_fastq "D:/dir/for/FASTQ/files"

Arguments:
--fasta - path to FASTA file with reference genome (.fasta, .fa)
--insert_size - size of a fragment to be sequenced, can't be smaller than read length
--read_length - length of a single read in pair
--coverage - number of mapped reads to cover single nucleotide in genome, the bigger to more accurate variant calling
--avg_quality - average quality of nucleotides (bases) of reads in FASTQ files
--snv_rate - single nucleotide variation rate
--ins_rate - insertion rate
--del_rate - deletion rate
--out_sam - path to directory for output SAM files of simulator, default ./out_sam
--out_fastq - path to directory for output FASTQ files of simulator, default ./out_fastq

## Comparison mode
Example of usage: --mode compare --sim "D:/path/to/simulator/SAM/file" --bwa_mem "D:/path/to/bwa mem/SAM/file" --bowtie "D:/path/to/bowtie/SAM/file" --gen_name "Mycobacterium tuberculosis"

Arguments:
--sim - path to SAM file that simulator generated
--bwa_mem - path to SAM file that BWA MEM aligner generated
--bowtie - path to SAM file that Bowtie2 aligner generated
