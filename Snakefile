sample = config["sample"]

rule align_read1:
    input:
        fastq=config["R1"]
    output:
        f"sam/{sample}_R1_aligned.sam"
    shell:
        """
        bowtie2 -x fasta/index_bowtie_transcriptome/prot_coding_ -U {input} -S {output} --no-unal
        """

rule extract_read1_names:
    input:
        f"sam/{sample}_R1_aligned.sam"
    output:
        f"sam/{sample}_R1_names.txt"
    shell:
        """
        awk '{{if ($1 !~ /^@/) print $1}}' {input} | sort | uniq > {output}
        """

rule filter_read2_mates:
    input:
        fastqR2 = config["R2"],
        fastqR1 = config["R1"],
        names = f"sam/{sample}_R1_names.txt"
    output:
        read2_aligned = f"fastq_read2/{sample}_R2_R1aligned.fq",
        read1_aligned = f"fastq_read1/{sample}_R1_R1aligned.fq"
    shell:
        """
        seqtk subseq {input.fastqR2} {input.names} > {output.read2_aligned}
        seqtk subseq {input.fastqR1} {input.names} > {output.read1_aligned}
        """

rule filter_read2_with_adapter_sequence:
    input:
        f"fastq_read2/{sample}_R2_R1aligned.fq"
    output:
        f"fastq_read2/{sample}_R2_adapter.fq"
    shell:
        """
        grep -E -A 2 -B 1 --no-group-separator "^[ATCGN]{{15}}GTCAG" {input} > {output}
        """

rule fastp_filtering:
    input:
        read1 = f"fastq_read1/{sample}_R1_R1aligned.fq",
        read2 = f"fastq_read2/{sample}_R2_adapter.fq"
    output:
        read1_fp = f"fastq_read1/{sample}_R1_fp.fq",
        read2_fp = f"fastq_read2/{sample}_R2_fp.fq",
        report = f"fastq_read2/{sample}_fp_report.html"
    shell:
        """
        fastp \
            -i {input.read1} \
            -I {input.read2} \
            -o {output.read1_fp} \
            -O {output.read2_fp} \
			-Q \
            --umi \
            --umi_loc=read2 \
            --umi_len=15 \
            --trim_front2=5 \
            -l 80 \
            -h {output.report}
        """		

rule align_R2_bowtie_genome_end_to_end:
	input: 
		f"fastq_read2/{sample}_R2_fp.fq"
	output: 
		f"R_input/R2_genome_mapped.tsv"
	shell:
		"""
		bowtie2 --no-unal --threads 16 -x fasta/index_bowtie_genome/GRCh38_primary_genome -U {input} -S sam/R2_aligned_to_genome.sam
		python R_scripts/extract_read_tags.py -i sam/R2_aligned_to_genome.sam -o {output}
		"""
	
rule align_R2_to_masked_genome:
	input:
		f"fastq_read2/{sample}_R2_fp.fq"
	output:
		f"R2_aligned_to_masked_genome/Aligned.out.sam"
	shell:
		"""
		STAR \
  --genomeDir fasta/index_masked_genome_star/ \
  --readFilesIn {input} \
  --runThreadN 12 \
  --outFileNamePrefix R2_aligned_to_masked_genome/ \
  --outSAMtype SAM \
  --outFilterMismatchNoverReadLmax 0.05 \
  --alignIntronMax 50000 \
  --outFilterMultimapNmax 1 \
  --outSAMunmapped Within
		"""
	
rule extract_columns_from_R1_R2:
    input:
        read1_sam = f"sam/{sample}_R1_aligned.sam",
        read2_sam = f"R2_aligned_to_masked_genome/Aligned.out.sam"
    output:
        read1_cols = f"R_input/{sample}_R1_aligned.tsv",
        read2_cols = f"R_input/{sample}_R2_aligned_to_masked_genome.tsv"
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} !/^@/ {{split($3, a, "|"); print $1, a[2], a[1], $2, $4, $5, $6}}' {input.read1_sam} > {output.read1_cols}
        awk 'BEGIN {{OFS="\t"}} $1 !~ /^@/ {{print $1, $2, $3, $4, $6, $10}}' {input.read2_sam} > {output.read2_cols}
        """

rule tail_analysis:
    input:
        tabR1 = "R_input/{sample}_R1_aligned.tsv",
        tabR2 = "R_input/{sample}_R2_aligned_to_masked_genome.tsv",
        tabR2_2 = "R_input/R2_genome_mapped.tsv"
    output:
        big_table = "R_output/{sample}_final_big.rds",
        small_table = "R_output/{sample}_final_small.rds"
    script:
        "R_scripts/run_tail_analysis.R"

