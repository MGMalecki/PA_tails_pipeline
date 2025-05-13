This pipeline should work with ungzipped fastq files located in directories fastq_R1/ fastq_R2/
This pipeline should still need genome indexes for human genome/transcriptome made with bowtie or STAR to be put in designated folders
All genome sequences and gtf file to generate STAR intex were from https://www.gencodegenes.org/human/
fasta/index_bowtie_transcriptome/prot_coding_
fasta/index_bowtie_genome/GRCh38_primary_genome
fasta/
For transcriptome file I used "gencode.v47.pc_transcripts.fa" and prefix "prot_coding_"
example code:

$ bowtie2-build --threads 16 gencode.v47.pc_transcripts.fa prot_coding_

The genome file I used GRCh38.primary_assembly.genome.fa and prefix "GRCh38_primary_genome"
example code:

$ bowtie2-build --threads 16 GRCh38.primary_assembly.genome.fa GRCh38_primary_genome

The bowtie generated indexes should have prefixes exactly as here or pipeline will not work 

Finally STAR generated index of masked genome (same genome file GRCh38.primary_assembly.genome.fa)

$ STAR --runThreadN 16 --runMode genomeGenerate --genomeDir index_masked_genome_star/ --genomeFastaFiles genome/GRCh38.primary_masked_AT_new.fa --sjdbGTFfile genome/gtf/gencode.v47.primary_assembly.basic.annotation.gtf --sjdbOverhang 79 --genomeSAindexNbases 14

To create masked genome I uploaded in OSF bed file that has locations of AT streatches that I identified using custom python script. Link to bed file:
https://osf.io/x2zeu/files/osfstorage

Once you have it mask the genome using bedtools:
Example code:

$ bedtools maskfasta -fi GRCh38.primary_assembly.genome.fa -bed AT_stretches_merged.bed -fo GRCh38.primary_masked_AT_new.fa

Next check if you cannot find streatches, for example like this:

grep -n 'TTTTTTTTTTTTTTTTTTTT' GRCh38.primary_masked_AT_new.fa | wc -l
grep -n 'AAAAAAAAAAAAAAAAAAAA' GRCh38.primary_masked_AT_new.fa | wc -l
grep -n 'AAAAAAAAAAAAAAAACAAA' GRCh38.primary_masked_AT_new.fa | wc -l
grep -n 'TTTTGTTTTTTTTTTTTTTT' GRCh38.primary_masked_AT_new.fa | wc -l

It should give no hits in masked genome but hits in original genome file.

ENVIORNMENT: I did not include conda bits in the snakemake, but this repo contains yaml file with software from my enviornment where this pipeline was working :)
before running use .yaml file to create appropriate enviornment and activate it.

FINALLY Once everything is in its place you should place fastq files in their folders and update config.yaml file to point input file names and sample name.
Run pipeline with command (-j should be adjusted to the system)
For new set of files it should be enough to add new fastq files to the folder and modify config.yaml and execute snakemake line again (please try;)

snakemake R_output/t47_final_small.rds --configfile config.yaml -j 16
