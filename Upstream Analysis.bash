#install the following topics 
Sudo apt install Sra-tools
Sudo apt install fastqc 
Sudo apt install Trimmomatic 
Sudo apt install Kallisto

#First Step: Download the Data from SRA using & Split Paird-end Files
# Download SRRSRR948835 Reads
Prefetch SRRSRR948835
Fastq-dump --split-files SRRSRR948835

# Download SRRSRR948836 Reads
Prefetch SRRSRR948836
Fastq-dump --split-files SRRSRR948836

#Second Step: Make Folder for Quality Check of the Reads
mkdir Qc_output

##Third Step: Quality Check of the Reads
Fastqc SRRSRR948836_R1 -o Qc_output 
Fastqc SRRSRR948836_R2 -o Qc_output 
Fastqc SRRSRR948835_R1 -o Qc_output 
Fastqc SRRSRR948835_R2 -o Qc_output

#Fourth Step: Create New Folder for Trimmed Sequence
mkdir Trimmed_output


#Fifth Step: Trimming the low quality sequence & adapters content using Trimmomatic
java -jar trimmomatic-0.39.jar PE -phred33  SRR948835_1.fastq SRR948835_2.fastq output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Sixth Step: Download the Reference Genome using Wget (Darmor Brassica Napus Transcriptome )
wget https://yanglab.hzau.edu.cn/static/bnir/assets/genomic_sequence/BnIRData/AACC.Brassica_napus/Darmor/v10/Brassica_napus.Darmor.v10.cds.fa.gz


#Seveth Step: Indexing the reference Transcriptome using Kallisto 
 kallisto index -i transcriptome.index Brassica_napus.ZS11.v0.cds.fa.gz

#Eighteth Step: Alignment & Quantification using Kalisto 
kallisto quant -i transcriptome.index -o output_dir -b 100 reads_1.fastq reads_2.fastq
