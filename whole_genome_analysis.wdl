version 1.0

# Task to perform quality control using FastQC on paired-end reads
task FastQC {
    input {
        File input_fastq1   # Input FASTQ file (R1)
        File input_fastq2   # Input FASTQ file (R2)
        String output_dir   # Directory to store FastQC results
        Int memory          # Memory allocation
        Int cpu             # CPU allocation
        Int disk            # Disk space required
    }

    command {
        mkdir -p ~{output_dir}
        fastqc ~{input_fastq1} ~{input_fastq2} --outdir ~{output_dir}
    }

    output {
        File fastqc_zip1 = "~{output_dir}/$(basename ~{input_fastq1} .fastq).zip"
        File fastqc_zip2 = "~{output_dir}/$(basename ~{input_fastq2} .fastq).zip"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "biocontainers/fastqc:v0.11.9-2"
    }
}

# Task to trim paired-end reads using Trimmomatic
task Trimmomatic {
    input {
        File input_fastq1    # Input FASTQ file (R1)
        File input_fastq2    # Input FASTQ file (R2)
        String output_dir    # Directory to store trimmed FASTQ results
        String adapter_file  # File containing adapter sequences
        Int memory           # Memory allocation
        Int cpu              # CPU allocation
        Int disk             # Disk space required
    }

    command {
        mkdir -p ~{output_dir}
        trimmomatic PE -phred33 ~{input_fastq1} ~{input_fastq2} \
            ~{output_dir}/trimmed_paired_1.fastq ~{output_dir}/trimmed_unpaired_1.fastq \
            ~{output_dir}/trimmed_paired_2.fastq ~{output_dir}/trimmed_unpaired_2.fastq \
            ILLUMINACLIP:~{adapter_file}:2:30:10 \
            SLIDINGWINDOW:4:20 \
            MINLEN:36
    }

    output {
        File trimmed_paired_1 = "~{output_dir}/trimmed_paired_1.fastq"
        File trimmed_paired_2 = "~{output_dir}/trimmed_paired_2.fastq"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "biocontainers/trimmomatic:v0.39-1"
    }
}

# Task to align paired-end reads using BWA-MEM
task BWAAlign {
    input {
        File trimmed_paired_1   # Input trimmed FASTQ file (R1)
        File trimmed_paired_2   # Input trimmed FASTQ file (R2)
        String reference        # Reference genome file
        String output_prefix    # Prefix for output BAM files
        Int memory              # Memory allocation
        Int cpu                 # CPU allocation
        Int disk                # Disk space required
    }

    command {
        bwa mem -t ~{cpu} ~{reference} ~{trimmed_paired_1} ~{trimmed_paired_2} | \
        samtools view -Sb - > ${output_prefix}.bam

        samtools sort -o ${output_prefix}_sorted.bam ${output_prefix}.bam
        samtools index ${output_prefix}_sorted.bam
    }

    output {
        File bam_file = "${output_prefix}_sorted.bam"
        File bai_file = "${output_prefix}_sorted.bam.bai"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "biocontainers/bwa:v0.7.17-2"
    }
}

# Task to preprocess the raw sequencing data using GATK's BQSR (Base Quality Score Recalibration)
task BQSR {
    input {
        File bam_file             # Input BAM file
        File known_sites          # Known sites file for BQSR (e.g., known SNPs, indels)
        String reference          # Reference genome file
        Int memory                # Memory allocation
        Int cpu                   # CPU allocation
        Int disk                  # Disk space required
    }

    command {
        gatk BaseRecalibrator \
            -I ~{bam_file} \
            -R ~{reference} \
            --known-sites ~{known_sites} \
            -O recalibration_report.grp

        gatk ApplyBQSR \
            -R ~{reference} \
            -I ~{bam_file} \
            --bqsr-recal-file recalibration_report.grp \
            -O recalibrated.bam
    }

    output {
        File recalibrated_bam = "recalibrated.bam"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "broadinstitute/gatk:latest"
    }
}

# Task to call variants using GATK HaplotypeCaller
task HaplotypeCaller {
    input {
        File recalibrated_bam     # Input BAM file after BQSR
        String reference          # Reference genome file
        String output_prefix      # Prefix for output VCF files
        Int memory                # Memory allocation
        Int cpu                   # CPU allocation
        Int disk                  # Disk space required
    }

    command {
        gatk HaplotypeCaller \
            -R ~{reference} \
            -I ~{recalibrated_bam} \
            -O ${output_prefix}.vcf.gz \
            -ERC GVCF
    }

    output {
        File vcf = "${output_prefix}.vcf.gz"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "broadinstitute/gatk:latest"
    }
}

# Task to perform variant joint calling using GATK GenotypeGVCFs
task GenotypeGVCFs {
    input {
        Array[File] gvcfs          # Array of GVCF files from HaplotypeCaller
        String reference           # Reference genome file
        String output_prefix       # Prefix for output VCF files
        Int memory                 # Memory allocation
        Int cpu                    # CPU allocation
        Int disk                   # Disk space required
    }

    command {
        gatk GenotypeGVCFs \
            -R ~{reference} \
            -V gendb://genomicDB \
            -O ${output_prefix}.vcf.gz
    }

    output {
        File vcf = "${output_prefix}.vcf.gz"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "broadinstitute/gatk:latest"
    }
}

# Task to perform variant filtration using GATK's VariantFiltration
task VariantFiltration {
    input {
        File vcf_file              # Input VCF file
        String reference          # Reference genome file
        String output_prefix      # Prefix for output VCF files
        Int memory                # Memory allocation
        Int cpu                   # CPU allocation
        Int disk                  # Disk space required
    }

    command {
        gatk VariantFiltration \
            -R ~{reference} \
            -V ~{vcf_file} \
            -O ${output_prefix}_filtered.vcf.gz \
            --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
            --filter-name "filter"
    }

    output {
        File vcf_filtered = "${output_prefix}_filtered.vcf.gz"
    }

    runtime {
        memory: memory + "G"
        cpu: cpu
        disk: disk + "GB"
        docker: "broadinstitute/gatk:latest"
    }
}

# Workflow for whole genome analysis including trimming, alignment, and GATK processing
workflow WholeGenomeAnalysis {
    input {
        File fastq_file1          # Input FASTQ file (R1)
        File fastq_file2          # Input FASTQ file (R2)
        File known_sites          # Known sites file for BQSR
        String reference         # Reference genome file
        String output_prefix     # Prefix for output files
        String adapter_file      # Adapter sequences file
        Int memory_trim          # Memory for trimming
        Int cpu_trim             # CPU for trimming
        Int disk_trim            # Disk for trimming
        Int memory_align         # Memory for alignment
        Int cpu_align            # CPU for alignment
        Int disk_align           # Disk for alignment
        Int memory_bqsr          # Memory for BQSR
        Int cpu_bqsr             # CPU for BQSR
        Int disk_bqsr            # Disk for BQSR
        Int memory_hc            # Memory for HaplotypeCaller
        Int cpu_hc               # CPU for HaplotypeCaller
        Int disk_hc              # Disk for HaplotypeCaller
        Int memory_gg            # Memory for GenotypeGVCFs
        Int cpu_gg               # CPU for GenotypeGVCFs
        Int disk_gg              # Disk for GenotypeGVCFs
        Int memory_vf            # Memory for VariantFiltration
        Int cpu_vf               # CPU for VariantFiltration
        Int disk_vf              # Disk for VariantFiltration
    }

    call FastQC {
        input:
            input_fastq1=fastq_file1,
            input_fastq2=fastq_file2,
            output_dir="fastqc_results",
            memory=memory_trim,
            cpu=cpu_trim,
            disk=disk_trim
    }

    call Trimmomatic {
        input:
            input_fastq1=fastq_file1,
            input_fastq2=fastq_file2,
            output_dir="trimmed_results",
            adapter_file=adapter_file,
            memory=memory_trim,
            cpu=cpu_trim,
            disk=disk_trim
    }

    call BWAAlign {
        input:
            trimmed_paired_1=Trimmomatic.trimmed_paired_1,
            trimmed_paired_2=Trimmomatic.trimmed_paired_2,
            reference=reference,
            output_prefix="aligned",
            memory=memory_align,
            cpu=cpu_align,
            disk=disk_align
    }

    call BQSR {
        input:
            bam_file=BWAAlign.bam_file,
            known_sites=known_sites,
            reference=reference,
            memory=memory_bqsr,
            cpu=cpu_bqsr,
            disk=disk_bqsr
    }

    call HaplotypeCaller {
        input:
            recalibrated_bam=BQSR.recalibrated_bam,
            reference=reference,
            output_prefix="raw_variants",
            memory=memory_hc,
            cpu=cpu_hc,
            disk=disk_hc
    }

    call GenotypeGVCFs {
        input:
            gvcfs=[HaplotypeCaller.vcf],
            reference=reference,
            output_prefix="genotyped_variants",
            memory=memory_gg,
            cpu=cpu_gg,
            disk=disk_gg
    }

    call VariantFiltration {
        input:
            vcf_file=GenotypeGVCFs.vcf,
            reference=reference,
            output_prefix="filtered_variants",
            memory=memory_vf,
            cpu=cpu_vf,
            disk=disk_vf
    }

    output {
        File final_vcf = VariantFiltration.vcf_filtered
    }
}

