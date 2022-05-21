# Illumina paired-end DNA sequencer simulator
A Python program that simulates DNA sequencing of paired-end Illumina reads and their aligning to the reference genome. It's useful for comparison of different aligners like BWA-MEM or Bowtie. </br>

### Illumina paired-end sequencing 
__1.__ Sampling DNA and cutting it to smaller fragments </br>
__2.__ Attaching adapters</br>
__3.__ PCR amplification</br>
__4.__ Attaching template to surface/flowcel</br>
__5.__ PCR/bridge amplification (cluster creation)</br>
__6.__ Adding fluorescent bases and taking a 
picture after each cycle (repeat this many 
times)</br>
__7.__ Stack up images and read the sequence</br>
![image](https://user-images.githubusercontent.com/76231958/169331402-92aa6bac-2388-4eb8-801e-a9acd37f7549.png)

### Simulation of Illumina paired-end sequencing
__1. FASTA file parsing__ - fasta file containes reference genome. Function from Bio SeqIO library was used for parsing FASTA files into dictionary object {sequence_name : sequence} </br>
__2. Generating reads__ - 1 pair of reads is generated from 1 DNA fragment. Positions of fragments were chosen randomly. Number of reads is calculated using Lander-Waterman equation _N = sequence_length * coverage / (2 * read_length)_. Reads are sequenced from each genome sequence separately, not from whole genome, that way every sequence is covered equaly. </br>
__3. Simulating errors__ - Reads are traveresed according to given error rates. Real-life sequencers make mistakes while generating reads so we need to simulate that. Also 
in real-life sequencing reference genome and sampled genome are not identical due to variations and mutations. </br>
__4. Writing reads to FASTQ files__ - First read in pair is written to _read_1.fastq and_ and second read in pair is written to _read_2.fastq_ like reveresed complement. </br>
__5. Generating SAM file__ - When simulating sequencing process read's mapping positions are well known, so paralel to writing reads to fastq files we can write alignments of paired reads to SAM file.

![image](https://user-images.githubusercontent.com/76231958/169332015-1312016d-8280-4ee9-86b0-64893cf72df6.png)


### Comparison
Both Bowtie SAM file and BWA MEM SAM file are first compared with simulator's SAM file and number of correctly aligned, misaligned and unmapped reads are calculated alongside mapping qaulity distribution of correctly aligned and misaligned reads. Mapping qualities are divided into 6 categories: 0-10, 10-20, 20-30, 30-40, 40-50, 50-60; since BWA MEM mapq range is [0,60] and Bowtie mapq range is [0, 42] so Bowtie mapq values are linearized in order to fit the BWA MEM mapping quality scale. Read is considered as correctly aligned only if its mapping position in tool's SAM file is equal to its position in simulator's SAM file. If that's not the case and read is not unmapped, then read is considered as misaligned. Alignments with 0x0004 flag set are unmapped. In terms of classification correctly aligned reads are considered as TRUE POSITIVES as both simulator and tool have mapped them the same way. Wrongly mapped reads or misaligned reads are considered as FALSE POSITIVES because tool has indeed mapped them but not at the right position (like simulator) and they can later lead to false positives variant calls. Unmapped reads are FALSE NEGATIVES since tool hasn't aligned them at all, nor correctly nor wrongly. Unmapped reads are not as critical as misaligned since they are easily filtered out and don't bring confusion in later stages of DNA analysis.

![image](https://user-images.githubusercontent.com/76231958/169332386-013c93f1-cb60-488c-869d-02ff3b0c5b6e.png)



# Usage
## Simulating mode(default)
<b>Example of usage:</b> ./source/dna_seq_sim_ana.py --mode sim --fasta <i>"C:/Users/user/path/to/FASTA/file.fasta"</i> --insert_size 400 --read_lenght 150 --coverage 3 --avg_quality 40 --snv_rate 0.009 
--ins_rate 0.009 --del_rate 0.009 --out_sam <i>"D:/dir/for/SAM/files"</i> --out_fastq <i>"D:/dir/for/FASTQ/files"</i> </br>

<b>Arguments:</b> </br>
__--mode__ - choosing mode, if mode == sim, then sequencing simulation is started </br>
<b>--fasta </b> - path to FASTA file with reference genome (.fasta, .fa) </br>
<b>--insert_size </b> - size of a fragment to be sequenced, can't be smaller than read length </br>
<b>--read_length</b> - length of a single read in pair <br/>
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
__--mode__ - choosing mode, if mode == compare then comparison is started
<b>--sim</b> - path to SAM file that simulator generated </br>
<b>--bwa_mem</b> - path to SAM file that BWA MEM aligner generated </br>
<b>--bowtie</b> - path to SAM file that Bowtie2 aligner generated </br>
<b>--gen_name</b> - name of reference genome
