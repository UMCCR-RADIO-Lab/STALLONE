# By Andrew Pattison at the UMCCR

configfile: "config.yaml"
# These require command line input
overhang="49"
# Can be zcat or '-' for uncompressed files
gzip="zcat"

# These won't change but the reference files will be required
ref_genome_dir="reference"
ref_genome_fa="reference/GRCh37.fa"
dbsnp="reference/dbsnp-151.vcf.gz"
exons="reference/exons_main_chrs_refseq.bed"
coding="reference/RefSeq_exons_coding_sorted_no_overlap_main_chrs.bed"
exac="reference/exac.vcf.gz"
vep_cache="reference/VEP_GRCh37_cache"
# Could maybe use 'which vep' to get this
vep_path="/data/cephfs/punim0010/extras/Pattison/miniconda2/envs/STALLONE/bin/"

# Define all the outputs I want
rule all:
    input:
        final="output/final/TMB_all_samples.csv"

# Run atropos to trim low quality reads
#The 'X{200}' args remove homopolymeric seqs at the 3' end of reads 
# rule atropos_trim:
#     input:
#         #lambda wildcards: config["samples"][wildcards.sample],
#         fq1="data/samples/{sample}.R1_001.fastq.gz",
#         fq2="data/samples/{sample}.R1_001.fastq.gz"
#     output:
#         fq1_trim=temp("data/samples/{sample}_trimmed_R1.fastq.gz"),
#         fq2_trim=temp("data/samples/{sample}_trimmed_R2.fastq.gz"),
#         report="output/trim_reports/{sample}.trim_report.json"
#     params:
#         samp="{sample}"
#     threads: 5
#     shell:
#         "atropos trim --threads {threads} --quality-base 33 --format fastq "
#         "-a 'A{{200}}' -a 'C{{200}}' -a 'G{{200}}$' -a 'T{{200}}' --overlap 8 "
#         "--no-default-adapters --no-cache-adapters -A 'A{{200}}' -A 'C{{200}}' "
#         "-A 'G{{200}}' -A 'T{{200}}' -pe1 {input.fq1}  -pe2 {input.fq2} "
#         "-o {output.fq1_trim} -p {output.fq2_trim} --report-file {output.report} "
#         "--report-formats json --sample-id {params.samp} --quality-cutoff=5 --minimum-length=25"

# STAR first pass to get splice junctions
rule STAR_pass1:
    input:
        "data/samples/{sample}.R1_001.fastq.gz",
        "data/samples/{sample}.R2_001.fastq.gz"
    output:
        outbam=temp("output/bams/{sample}.bam"),
        sjout=temp("output/bams/{sample}SJ.out.tab")
    params:
        tempdir="output/bams/{sample}_temp",
        out_file_prefix="output/bams/{sample}",
        genomedir=ref_genome_dir,
        gz=gzip
    threads: 5
    shell:
        "STAR --genomeDir {params.genomedir} --readFilesCommand {params.gz} --outTmpDir {params.tempdir} "
        "--runThreadN {threads} --readFilesIn {input} --outFileNamePrefix {params.out_file_prefix} "
        "--outStd BAM_Unsorted > {output.outbam}"

# Generate a new genome for STAR with added splice junctions
rule STAR_genome_generate:
    input:
        "output/bams/{sample}SJ.out.tab",
    output:
        "output/genomes/{sample}/mockfile.txt"        
    params:
        ref_fa=ref_genome_fa,
        sjdb_overhang=overhang,
        tempdir="output/genomes/{sample}_temp2",
        gdir="output/genomes/{sample}"
    threads: 5
    shell:
        "STAR --runMode genomeGenerate --genomeDir {params.gdir} --genomeFastaFiles {params.ref_fa} --outTmpDir {params.tempdir} "
        "--sjdbFileChrStartEnd {input} --sjdbOverhang {params.sjdb_overhang} "
        "--runThreadN {threads} --limitOutSJcollapsed 500000 && touch {output}"

# Second pass including updated splice junctions
rule STAR_pass2:
    input:
        fq1="data/samples/{sample}.R1_001.fastq.gz",
        fq2="data/samples/{sample}.R2_001.fastq.gz",
        mockfile="output/genomes/{sample}/mockfile.txt"
    output:
        temp("output/bams/{sample}_second_pass.bam")
    params:
        tempdir="output/genomes/{sample}_temp3",
        genomedir="output/genomes/{sample}",
        gz=gzip,
        out_file_prefix="output/bams/{sample}_second_pass_"
    threads: 5
    shell:
        "STAR --outSAMtype BAM Unsorted --genomeDir {params.genomedir} --readFilesCommand {params.gz} --outTmpDir {params.tempdir} "
        "--readFilesIn {input.fq1} {input.fq2} --runThreadN {threads} --outFileNamePrefix {params.out_file_prefix} --outStd BAM_Unsorted > {output}"

# Add read groups that GATK likes with picard
rule Picard_RG:
    input:
        "output/bams/{sample}_second_pass.bam"
    output:
        temp("output/bams/{sample}_rg_added_sorted.bam")
    shell:
        "picard AddOrReplaceReadGroups I={input} O={output} "
        "SO=coordinate RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM={wildcards.sample}"
         
# Mark duplicates with Picard
# SORTING_COLLECTION_SIZE_RATIO reduced becuase I was running out of RAM
rule Picard_MD:
    input:
        "output/bams/{sample}_rg_added_sorted.bam"
    output:
        temp("output/bams/{sample}_dups_marked.bam")
    params:
        metrics="output/bams/{sample}_marked_dup_metrics.txt"
    shell:
        "picard MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.1 I={input} O={output} "
        "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={params.metrics}"

# Reorder the bam with picard
rule Picard_RB:
    input:
        "output/bams/{sample}_dups_marked.bam"
    output:
        temp("output/bams/{sample}_ordered.bam")
    params:
        ref_fa=ref_genome_fa
    shell:
        "picard ReorderSam REFERENCE={params.ref_fa} I={input} O={output} " 
        "VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"    

# Clip reads from splice junctions overhanging introns 
rule GATK_split:
    input:
        "output/bams/{sample}_ordered.bam"
    output:
        temp("output/bams/{sample}_split.bam")
    params:
        ref_fa=ref_genome_fa
    shell:
        "gatk SplitNCigarReads -R {params.ref_fa} -I {input} -O {output}"

# Reclibrate the bases of the bam file to correct for errors across a flow cell
rule GATK_BQSR_calc:
    input:
        "output/bams/{sample}_split.bam"
    output:
        temp("output/bams/{sample}_recal_data.table")
    params:
        dbsnp=dbsnp,
        ref_fa=ref_genome_fa
    shell:
        "gatk BaseRecalibrator -R {params.ref_fa} -I {input} --known-sites {params.dbsnp} -O {output}"

# Apply the base recalibration
rule GATK_apply_BQSR:
    input:
        split_bam="output/bams/{sample}_split.bam",
        bqsr_table="output/bams/{sample}_recal_data.table"
    output:
        "output/bams/{sample}_recal.bam"
    params:
        ref_fa=ref_genome_fa,
    shell:
        "gatk ApplyBQSR -R {params.ref_fa} -I {input.split_bam} -bqsr {input.bqsr_table} -O {output}"

# Get the coverage for the sample
# This function gets coverage from the coding portion of the genome only
# To more accurately calculate TMB
# We need at least 4 bases coverage since we need at least 4 altered bases to call a variant
rule Get_BAM_coverage:
    input:
        "output/bams/{sample}_recal.bam"
    output:
        cov_file="output/coverage/{sample}.quantized.bed.gz"
    params:
        coding=coding,
        prefix="output/coverage/{sample}"
    threads: 5
    shell:
        "mosdepth -t {threads} --by {params.coding} --quantize 4: {params.prefix} {input}"

# Get the total amount of the sample sufficiently covered to call variants
rule Sum_coverage:
    input:
        "output/coverage/{sample}.quantized.bed.gz"
    output:
        "output/total_coverage/{sample}_bases_covered.txt"        
    params:
        coding=coding,
    threads: 5
    shell:
        "gunzip {input} | awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}' >> {output}"
    
# GATK haplotpye caller
rule GATK_HC:
    input:
        "output/bams/{sample}_recal.bam"
    output:
        "output/vcfs/{sample}_GATK_HC_raw.vcf"
    params:
        ref_fa=ref_genome_fa,
        dbsnp=dbsnp,
        exons=exons
    threads: 5
    shell:
        "gatk HaplotypeCaller -R {params.ref_fa} -I {input} --intervals {params.exons} --dont-use-soft-clipped-bases true "
        "--native-pair-hmm-threads {threads} --standard-min-confidence-threshold-for-calling 20.0 "
        "-O {output} --dbsnp {params.dbsnp}"

# Take the resuls of the GATK haplotype caller and run vcf2maf
# I can then do the obvious filtering and then read them into R for some custom filtering
rule vcf2maf:
    input:
        "output/vcfs/{sample}_GATK_HC_raw.vcf"
    output:
        "output/raw_mafs/{sample}.maf"
    params:
        ref_fa=ref_genome_fa,
        exac=exac,
        vep_cache=vep_cache,
        vep_path=vep_path
    threads: 5
    shell:
        "vcf2maf.pl --input-vcf {input} --output-maf {output} --ref-fasta {params.ref_fa} --filter-vcf {params.exac} "
        "--vep-path {params.vep_path} --vep-data {params.vep_cache} --vep-forks 8 --tumor-id {wildcards.sample} "
        "--species homo_sapiens --retain-info DB,MQ,QD,FS,AD,DP"

# Make a rule to filter the raw mafs for quicker/custom R processing
rule Filter_raw_maf:
    input:
        "output/raw_mafs/{sample}.maf"
    output:
        "output/filtered_mafs/{sample}_loose_filter.maf"
    shell:
        "Rscript scripts/Filter_raw_maf.R {input} {output}"

# Run the R script to generate somatic variants and TMB
rule Get_somatic_variants:
    input:
        coverage=expand("output/total_coverage/{sample}_bases_covered.txt", sample=config["samples"]),
        mafs=expand("output/filtered_mafs/{sample}_loose_filter.maf", sample=config["samples"])
    output:
        maf="output/final/Somatic_varaints_all_samples.maf",
        TMB="output/final/TMB_all_samples.csv"
    shell:
        "Rscript scripts/Knockout_unwanted_variants.R {output.TMB} {output.maf} {input.mafs} {input.coverage}"









