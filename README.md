# Illumina-paired-end-DNA-sequencer
A Python program that simulates DNA sequencing of paired-end Illumina reads and their aligning to the reference genome. It's useful for comparison of different aligners like BWA-MEM or Bowtie. 

# Usage
## Simulating mode(default)
<b>Example of usage:</b> ./source/dna_seq_sim_ana.py --mode sim --fasta <i>"C:/Users/user/path/to/FASTA/file.fasta"</i> --insert_size 400 --read_lenght 150 --coverage 3 --avg_quality 40 --snv_rate 0.009 
--ins_rate 0.009 --del_rate 0.009 --out_sam <i>"D:/dir/for/SAM/files"</i> --out_fastq <i>"D:/dir/for/FASTQ/files"</i> </br>

<b>Arguments:</b> </br>
<b>--fasta </b> - path to FASTA file with reference genome (.fasta, .fa) </br>
<b>--insert_size </b> - size of a fragment to be sequenced, can't be smaller than read length </br>
<b>--read_length<?b> - length of a single read in pair <br/>
<b>--coverage</b> - number of mapped reads to cover single nucleotide in genome, the bigger to more accurate variant calling </br>
<b>--avg_quality</b> - average quality of nucleotides (bases) of reads in FASTQ files </br>
<b>--snv_rate</b> - single nucleotide variation rate </br>
<b>--ins_rate</b> - insertion rate </br>
<b>--del_rate</b> - deletion rate <br>
<b>--out_sam</b> - path to directory for output SAM files of simulator, default ./out_sam </br>
<b>--out_fastq</b> - path to directory for output FASTQ files of simulator, default ./out_fastq <br>

## Comparison mode
<b>Example of usage:</b> --mode compare --sim <i>"D:/path/to/simulator/SAM/file"</i> --bwa_mem <i>"D:/path/to/bwa mem/SAM/file"</i> --bowtie <i>"D:/path/to/bowtie/SAM/file"</i> --gen_name <i>"Mycobacterium tuberculosis"</i> </br>

<b>Arguments:</b> </br>
<b>--sim</b> - path to SAM file that simulator generated </br>
<b>--bwa_mem</b> - path to SAM file that BWA MEM aligner generated </br>
<b>--bowtie</b> - path to SAM file that Bowtie2 aligner generated </br>
