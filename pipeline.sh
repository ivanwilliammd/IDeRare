#!/bin/bash
set -euo pipefail

echo "STEP 0 : Setup Variable, Directory, and Ensure NVIDIA Docker Supported"
### WARNING : Don't use tilde (~) in the path

# -------------------
## Step 0a: Set variable necessary with sample name for RG
# -------------------
echo "STEP 0a: Set variable necessary"

BIN_VERSION=1.5.0
memory=60G
DEEPTRIO_MODEL=WES

glnexus_version=1.4.1
tiddit_version=3.6.1--py38h24c8ff8_0

# Father
father_name=V350145665_L04_B5EHOMdmhwXAACA-519

# Mother
mother_name=V350145665_L04_B5EHOMdmhwXAABA-517

# Proband
proband_name=V350145665_L04_B5EHOMdmhwXAAAA-515

proband_SM=Proband
proband_PU=AgilentV6
proband_PL=DNBSeq
proband_LB=WXS

exomiser_solo=V350145665_L04_B5EHOMdmhwXAAAA-515_exomiser_solo.yml
exomiser_solo_sv=V350145665_L04_B5EHOMdmhwXAAAA-515_tiddit_exomiser_solo.yml
exomiser_trio=V350145665_L04_B5EHOMdmhwXAAAA-515_exomiser_trio.yml

mother_SM=Mother
mother_PU=AgilentV6
mother_PL=DNBSeq
mother_LB=WXS

father_SM=Father
father_PU=AgilentV6
father_PL=DNBSeq
father_LB=WXS

trimming=false
solo_analysis=true
trio_analysis=true

# Check if BIN_VERSION is null
if [ -z "$BIN_VERSION" ]; then
    echo "BIN_VERSION is null, stopping the pipeline."
    exit 1
fi

# Check if memory is null
if [ -z "$memory" ]; then
    echo "memory is null, stopping the pipeline."
    exit 1
fi

# Check if DEEPTRIO_MODEL is null
if [ -z "$DEEPTRIO_MODEL" ]; then
    echo "DEEPTRIO_MODEL is null, stopping the pipeline."
    exit 1
fi

# Check if glnexus_version is null
if [ -z "$glnexus_version" ]; then
    echo "glnexus_version is null, stopping the pipeline."
    exit 1
fi

# Check if tiddit_version is null
if [ -z "$tiddit_version" ]; then
    echo "tiddit_version is null, stopping the pipeline."
    exit 1
fi

# Check if proband_name is null
if [ -z "$proband_name" ]; then
    echo "Input the filename for proband."
    exit 1
fi

# If trio_analysis is not null, check if mother_name and father_name is null
if [ "$trio_analysis" = true ]; then
    if [ -z "$mother_name" ]; then
        echo "Input the filename for mother."
        exit 1
    fi

    if [ -z "$father_name" ]; then
        echo "Input the filename for father."
        exit 1
    fi
fi

# -------------------
## Step 0b: Prepare directory and pedigree file
# -------------------
echo "STEP 0b: Prepare directory and pedigree file"

## Pedigree file on INPUT_DIR (Just file name only)
INPUT_DIR="/home/ivanwilliamharsono/Downloads/trio_patient/input"
PEDIGREE="trio.ped"

## Subfolder of INPUT_DIR
FASTQ_DIR=${INPUT_DIR}/A_FASTQ
SAM_DIR=${INPUT_DIR}/B_RAW_SAM_BAM

## Output Dir
OUTPUT_DIR="/home/ivanwilliamharsono/Downloads/trio_patient/output"
ANNOTATED_DIR="/home/ivanwilliamharsono/Downloads/trio_patient/annotated"
SV_DIR="/home/ivanwilliamharsono/Downloads/trio_patient/sv_tiddit"
EXOMISER_DIR="/home/ivanwilliamharsono/Downloads/trio_patient/exomiser"

## FASTA file on REFERENCE_DIR (Just file name only)
REFERENCE_DIR=/home/ivanwilliamharsono/Downloads/Database
REFSEQ_FASTA=hg38.p14.fa

snpEff_dir=/home/ivanwilliamharsono/Downloads/Sandbox/snpEff/data
snpEff_ver=GRCh38.p14
dbnsfp=/home/ivanwilliamharsono/Downloads/Database/dbNSFPv4.5a_custombuild.gz
dbSNP=/home/ivanwilliamharsono/Downloads/Database/Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz
ClinVar=/home/ivanwilliamharsono/Downloads/Database/clinvar_20240215.vcf.gz
chr_rename=/home/ivanwilliamharsono/Downloads/IDeRare/convert/vcf_rename.chrnames

cd /home/ivanwilliamharsono/Downloads/trio_patient

## Make Directory
mkdir -p ${INPUT_DIR} ${FASTQ_DIR} ${SAM_DIR} ${OUTPUT_DIR} ${OUTPUT_DIR}/intermediate_results_dir_proband ${OUTPUT_DIR}/intermediate_results_dir_trio ${ANNOTATED_DIR} ${SV_DIR} ${EXOMISER_DIR}

# -------------------
# STEP 1: QC - Run fastqp 
# -------------------
echo "STEP 1: QC - Run fastqp - Optional (if the input is raw untrimmed fastq)"

# Check if proband_name is null and trimming is true
if [ ! -z "$proband_name" ] && [ "$trimming" = true ]; then
    echo "STEP 1a : Proband"

    if [ -f "${FASTQ_DIR}/${proband_name}_1.fastq" ] && [ -f "${FASTQ_DIR}/${proband_name}_2.fastq" ]; then
        echo "Raw fastq exist, continue to Trimming"
    else
        echo "File not exist, please check the input directory"
        exit 1
    fi

    if [ -f "${FASTQ_DIR}/${proband_name}_1.fq.gz" ] && [ -f "${FASTQ_DIR}/${proband_name}_2.fq.gz" ]; then
        echo "Output fq.qz exist, skipping Trimming"
    else
        echo "File not exist, continue to Trimming"
        fastp -g -x -w $(nproc) \
            -D --dup_calc_accuracy 6 \
            --in1 ${FASTQ_DIR}/${proband_name}_1.fastq \
            --in2 ${FASTQ_DIR}/${proband_name}_2.fastq \
            --out1 ${FASTQ_DIR}/${proband_name}_1.fq.gz \
            --out2 ${FASTQ_DIR}/${proband_name}_2.fq.gz \
            -h ${FASTQ_DIR}/${proband_name}.html \
            -j ${FASTQ_DIR}/${proband_name}.json \
            -R ${proband_name}-${proband_SM}
    fi
fi


# Check if mother_name is null and trimming is true
if [ ! -z "$mother_name" ] && [ "$trimming" = true ] && [ "$trio_analysis" = true ]; then
    echo "STEP 1b : Mother"
    
    if [ -f "${FASTQ_DIR}/${mother_name}_1.fastq" ] && [ -f "${FASTQ_DIR}/${mother_name}_2.fastq" ]; then
        echo "Raw fastq exist, continue to Trimming"
    else
        echo "File not exist, please check the input directory"
        exit 1
    fi

    if [ -f "${FASTQ_DIR}/${mother_name}_1.fq.gz" ] && [ -f "${FASTQ_DIR}/${mother_name}_2.fq.gz" ]; then
        echo "Output fq.qz exist, skipping Trimming"
    else
        echo "File not exist, continue to Trimming"

        fastp -g -x -w $(nproc) \
            -D --dup_calc_accuracy 6 \
            --in1 ${FASTQ_DIR}/${mother_name}_1.fastq \
            --in2 ${FASTQ_DIR}/${mother_name}_2.fastq \
            --out1 ${FASTQ_DIR}/${mother_name}_1.fq.gz \
            --out2 ${FASTQ_DIR}/${mother_name}_2.fq.gz \
            -h ${FASTQ_DIR}/${mother_name}.html \
            -j ${FASTQ_DIR}/${mother_name}.json \
            -R ${mother_name}-${mother_SM}
    fi
fi

# Check if father_name is null and trimming is true
if [ ! -z "$father_name" ] && [ "$trimming" = true ] && [ "$trio_analysis" = true ]; then
    echo "STEP 1C : Father"
    
    if [ -f "${FASTQ_DIR}/${father_name}_1.fastq" ] && [ -f "${FASTQ_DIR}/${father_name}_2.fastq" ]; then
        echo "Raw fastq exist, continue to Trimming"
    else
        echo "File not exist, please check the input directory"
        exit 1
    fi

    if [ -f "${FASTQ_DIR}/${father_name}_1.fq.gz" ] && [ -f "${FASTQ_DIR}/${father_name}_2.fq.gz" ]; then
        echo "Output fq.qz exist, skipping Trimming"
    else
        echo "File not exist, continue to Trimming"
        fastp -g -x -w $(nproc) \
            -D --dup_calc_accuracy 6 \
            --in1 ${FASTQ_DIR}/${father_name}_1.fastq \
            --in2 ${FASTQ_DIR}/${father_name}_2.fastq \
            --out1 ${FASTQ_DIR}/${father_name}_1.fq.gz \
            --out2 ${FASTQ_DIR}/${father_name}_2.fq.gz \
            -h ${FASTQ_DIR}/${father_name}.html \
            -j ${FASTQ_DIR}/${father_name}.json \
            -R ${father_name}-${father_SM}
    fi
fi

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM2
# --------------------------------------
echo "STEP 2: Map to reference using BWA-MEM2"


# Check if proband_name is not null
if [ ! -z "$proband_name" ]; then
    echo "STEP 2a : Proband"

    if [ -f "${FASTQ_DIR}/${proband_name}_1.fq.gz" ] && [ -f "${FASTQ_DIR}/${proband_name}_2.fq.gz" ]; then
        echo "Paired end fq.qz exist, continue to BWA-MEM2"
    else
        echo "File not exist, please check the input directory"
        exit 1
    fi

    if [-f "${SAM_DIR}/${proband_name}_raw.sam" ] || [ -f "${INPUT_DIR}/${proband_name}.bam" ]; then
        echo "SAM file or final BAM file exist, skipping BWA-MEM2"
    else
        # BWA MEM2 Alignment
        bwa-mem2 mem -R "@RG\tID:${proband_name}\tSM:${proband_SM}\tPU:${proband_PU}\tPL:${proband_PL}\tLB:${proband_LB}" \
            -t $(nproc) ${REFERENCE_DIR}/${REFSEQ_FASTA} \
            ${FASTQ_DIR}/${proband_name}_1.fq.gz \
            ${FASTQ_DIR}/${proband_name}_2.fq.gz \
            > ${SAM_DIR}/${proband_name}_raw.sam
    fi
    
fi


# Check if mother_name is not null and trio_analysis is true
if [ ! -z "$mother_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 2b : Mother"

    if [ -f "${FASTQ_DIR}/${mother_name}_1.fq.gz" ] && [ -f "${FASTQ_DIR}/${mother_name}_2.fq.gz" ]; then
        echo "Paired end fq.qz exist, continue to BWA-MEM2"
    else
        echo "File not exist, please check the input directory"
        exit 1
    fi
    
    if [-f "${SAM_DIR}/${mother_name}_raw.sam" ] || [ -f "${INPUT_DIR}/${mother_name}.bam" ]; then
        echo "SAM file or final BAM file exist, skipping BWA-MEM2"
    else
        # BWA MEM2 Alignment
        bwa-mem2 mem -R "@RG\tID:${mother_name}\tSM:${mother_SM}\tPU:${mother_PU}\tPL:${mother_PL}\tLB:${mother_LB}" \
            -t $(nproc) ${REFERENCE_DIR}/${REFSEQ_FASTA} \
            ${FASTQ_DIR}/${mother_name}_1.fq.gz \
            ${FASTQ_DIR}/${mother_name}_2.fq.gz \
            > ${SAM_DIR}/${mother_name}_raw.sam
    fi
fi



# Check if father_name is not null and trio_analysis is true
if [ ! -z "$father_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 2c : Father"
    
    if [ -f "${FASTQ_DIR}/${father_name}_1.fq.gz" ] && [ -f "${FASTQ_DIR}/${father_name}_2.fq.gz" ]; then
        echo "Paired end fq.qz exist, continue to BWA-MEM2"
    else
        echo "File not exist, please check the input directory"
        exit 1
    fi

    if [-f "${SAM_DIR}/${father_name}_raw.sam" ] || [ -f "${INPUT_DIR}/${father_name}.bam" ]; then
        echo "SAM file or final BAM file exist, skipping BWA-MEM2"
    else
        # BWA MEM2 Alignment
        bwa-mem2 mem -R "@RG\tID:${father_name}\tSM:${father_SM}\tPU:${father_PU}\tPL:${father_PL}\tLB:${father_LB}" \
            -t $(nproc) ${REFERENCE_DIR}/${REFSEQ_FASTA} \
            ${FASTQ_DIR}/${father_name}_1.fq.gz \
            ${FASTQ_DIR}/${father_name}_2.fq.gz \
            > ${SAM_DIR}/${father_name}_raw.sam
    fi
fi

# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - sambamba
# -----------------------------------------
echo "STEP 3: Mark Duplicates and Sort - Sambamba"

# Conversion of SAM to BAM & Markdup

# Check if proband_name is not null
# 2nd pass : check if ${INPUT_DIR}/${proband_name}.bam not exist yet --> then run the function, else skipped

if [ ! -z "$proband_name" ]; then
    echo "STEP 3a : Proband"

    if [ -f "${INPUT_DIR}/${proband_name}.bam" ]; then
        echo "Final BAM file exist, skipping Mark Duplicates and Sort"
    else
        sambamba view -p -t=$(nproc) -l=9 \
            -S ${SAM_DIR}/${proband_name}_raw.sam \
            -f=bam -o=${SAM_DIR}/${proband_name}_raw.bam

        sambamba markdup -r -p -t=$(nproc) -l=9 \
            ${SAM_DIR}/${proband_name}_raw.bam \
            ${SAM_DIR}/${proband_name}_dedup.bam

        sambamba sort -m=${memory} -p -t=$(nproc) -l=9 \
            ${SAM_DIR}/${proband_name}_dedup.bam \
            -o=${INPUT_DIR}/${proband_name}.bam
    fi
fi

# Check if mother_name is not null and trio_analysis is true
if [ ! -z "$mother_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 3b : Mother"

    if [ -f "${INPUT_DIR}/${mother_name}.bam" ]; then
        echo "Final BAM file exist, skipping Mark Duplicates and Sort"
    else
        sambamba view -p -t=$(nproc) -l=9 \
            -S ${SAM_DIR}/${mother_name}_raw.sam \
            -f=bam -o=${SAM_DIR}/${mother_name}_raw.bam

        sambamba markdup -r -p -t=$(nproc) -l=9 \
            ${SAM_DIR}/${mother_name}_raw.bam \
            ${SAM_DIR}/${mother_name}_dedup.bam

        sambamba sort -m=${memory} -p -t=$(nproc) -l=9 \
            ${SAM_DIR}/${mother_name}_dedup.bam \
            -o=${INPUT_DIR}/${mother_name}.bam
    fi
fi

# Check if father_name is not null and trio_analysis is true
if [ ! -z "$father_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 3c : Father"
    
    if [ -f "${INPUT_DIR}/${father_name}.bam" ]; then
        echo "Final BAM file exist, skipping Mark Duplicates and Sort"
    else
        sambamba view -p -t=$(nproc) -l=9 \
            -S ${SAM_DIR}/${father_name}_raw.sam \
            -f=bam -o=${SAM_DIR}/${father_name}_raw.bam

        sambamba markdup -r -p -t=$(nproc) -l=9 \
            ${SAM_DIR}/${father_name}_raw.bam \
            ${SAM_DIR}/${father_name}_dedup.bam

        sambamba sort -m=${memory} -p -t=$(nproc) -l=9 \
            ${SAM_DIR}/${father_name}_dedup.bam \
            -o=${INPUT_DIR}/${father_name}.bam
    fi
fi

### Remove Intermediate SAM BAM file as it consumes too much spaces

# Check if the directory exists
if [ -d "$SAM_DIR" ]; then
    echo "Remove all data inside SAM directory to save harddisk space"
    rm -r "$SAM_DIR"
else
    echo "Directory $SAM_DIR does not exist."
fi

# ----------------------------------------------
# STEP 4: Variant Calling
# ----------------------------------------------

if [ ! -z "$proband_name" ] && [ "$solo_analysis" = true ]; then
    echo "STEP 4a: Variant Calling Proband DeepVariant"

    # Further check if ${INPUT_DIR}/${proband_name}.bam all exists
    if [ -f "${INPUT_DIR}/${proband_name}.bam" ]; then
        if [ -f "${OUTPUT_DIR}/${proband_name}_proband.vcf.gz" ]; then
            echo "VCF file exist, skipping DeepVariant"
        else
            echo "BAM file exist, continue to DeepVariant"
            docker run --gpus 1 \
                -v "${INPUT_DIR}":"/input"   \
                -v "${OUTPUT_DIR}":"/output"  \
                -v "${REFERENCE_DIR}":"/reference" \
                google/deepvariant:"${BIN_VERSION}-gpu" \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type ${DEEPTRIO_MODEL} \
                --ref /reference/${REFSEQ_FASTA} \
                --reads /input/${proband_name}.bam \
                --num_shards $(nproc)  \
                --intermediate_results_dir /output/intermediate_results_dir_proband \
                --output_gvcf /output/${proband_name}_proband.g.vcf.gz \
                --output_vcf /output/${proband_name}_proband.vcf.gz

            ## Remove Intermediate_results_dir to save spaces
            rm -r -f ${OUTPUT_DIR}/intermediate_results_dir_proband

            ## Stop Docker to save spaces and memory
            docker stop $(docker ps -aq)
            docker rm $(docker ps -aq)
        fi
    else
        echo "BAM file not exist, please check the input directory 3a"
        exit 1
    fi
fi

if [ ! -z "$proband_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 4b: Variant Calling DeepTrio"

    # Further check if ${INPUT_DIR}/${father_name}.bam and ${INPUT_DIR}/${mother_name}.bam and ${INPUT_DIR}/${proband_name}.bam all exists
    if [ -f "${INPUT_DIR}/${father_name}.bam" ] && [ -f "${INPUT_DIR}/${mother_name}.bam" ] && [ -f "${INPUT_DIR}/${proband_name}.bam" ]; then

        if [ -f "${OUTPUT_DIR}/${proband_name}.g.vcf.gz" ] && [ -f "${OUTPUT_DIR}/${father_name}.g.vcf.gz" ] && [ -f "${OUTPUT_DIR}/${mother_name}.g.vcf.gz" ]; then
            echo "gVCF file exist, skipping DeepTrio"
        else
            echo "All BAM files exist, continue to DeepTrio"
            docker run --gpus 1 \
                -v "${INPUT_DIR}":"/input"   \
                -v "${OUTPUT_DIR}":"/output"  \
                -v "${REFERENCE_DIR}":"/reference" \
                google/deepvariant:deeptrio-"${BIN_VERSION}-gpu" \
                /opt/deepvariant/bin/deeptrio/run_deeptrio \
                --model_type ${DEEPTRIO_MODEL} \
                --ref /reference/${REFSEQ_FASTA} \
                --reads_child /input/${proband_name}.bam \
                --reads_parent1 /input/${father_name}.bam \
                --reads_parent2 /input/${mother_name}.bam \
                --output_vcf_child /output/${proband_name}.output.vcf.gz \
                --output_vcf_parent1 /output/${father_name}.output.vcf.gz \
                --output_vcf_parent2 /output/${mother_name}.output.vcf.gz \
                --sample_name_child 'child' \
                --sample_name_parent1 'father' \
                --sample_name_parent2 'mother' \
                --num_shards $(nproc)  \
                --intermediate_results_dir /output/intermediate_results_dir_trio \
                --output_gvcf_child /output/${proband_name}.g.vcf.gz \
                --output_gvcf_parent1 /output/${father_name}.g.vcf.gz \
                --output_gvcf_parent2 /output/${mother_name}.g.vcf.gz

            ## Remove Intermediate_results_dir to save spaces
            rm -r -f ${OUTPUT_DIR}/intermediate_results_dir_trio

            ## Stop Docker to save spaces and memory
            docker stop $(docker ps -aq)
            docker rm $(docker ps -aq)
        fi
    else
        echo "One of the BAM files not exist, please check the input proband.bam, father.bam, mother.bam or redo the step 3a, 3b, 3c"
        exit 1
    fi
fi

# ----------------------------------------------
# STEP 5: Merge gVCF files with GLnexus
# ----------------------------------------------
# Check if solo_analysis=true then skip, if trio_analysis=true then continue
if [ ! -z "$proband_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 5: Merge gVCF files with GLnexus"

    if [ -f "${OUTPUT_DIR}/${proband_name}_trio_merged.vcf.gz" ]; then
        echo "Merged VCF file exist, skipping GLnexus"
    else
        docker run \
            -v "${OUTPUT_DIR}":"/output" \
            ghcr.io/dnanexus-rnd/glnexus:v${glnexus_version} \
            /usr/local/bin/glnexus_cli \
            --config DeepVariant_unfiltered \
            /output/${proband_name}.g.vcf.gz \
            /output/${father_name}.g.vcf.gz \
            /output/${mother_name}.g.vcf.gz \
            | bcftools view -Oz -o ${OUTPUT_DIR}/${proband_name}_trio_merged.vcf.gz 

        ## Stop Docker to save spaces and memory
        docker stop $(docker ps -aq)
        docker rm $(docker ps -aq)
    fi

    # ----------------------------------------------
    # STEP 6: Calculate Mendelian Violation Rate using RTG Tools
    # ----------------------------------------------
    echo "STEP 6: Calculate Mendelian Violation Rate using RTG Tools"
    #### Reference : https://www.animalgenome.org/bioinfo/resources/manuals/RTGOperationsManual.pdf


    if [ -f "${OUTPUT_DIR}/${proband_name}_trio_merged.vcf.gz" ]; then
        echo "Merged VCF file exist, continue to RTG Tools"

        if [ -d ${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf ];
        then
            echo "${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf exists. Didn't create sdf."
        else
            echo "${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf directory does not exist. Creating it"
            docker run \
                -v "${INPUT_DIR}":"/input" \
                -v "${REFERENCE_DIR}":"/reference" \
                realtimegenomics/rtg-tools format \
                -o /reference/${REFSEQ_FASTA}.sdf "/reference/${REFSEQ_FASTA}"
        fi


        docker run \
            -v "${INPUT_DIR}":"/input" \
            -v "${REFERENCE_DIR}":"/reference" \
            -v "${OUTPUT_DIR}":"/output" \
            realtimegenomics/rtg-tools vcfstats \
            "/output/${proband_name}_trio_merged.vcf.gz" \
            | tee output/deepvariant.${proband_name}_trio.vcfstats.txt

        docker run \
            -v "${INPUT_DIR}":"/input" \
            -v "${REFERENCE_DIR}":"/reference" \
            -v "${OUTPUT_DIR}":"/output" \
            realtimegenomics/rtg-tools mendelian \
            -i "/output/${proband_name}_trio_merged.vcf.gz" \
            -o "/output/${proband_name}_trio_annotated.output.vcf.gz" \
            --pedigree=/input/trio.ped \
            -t /reference/${REFSEQ_FASTA}.sdf \
            | tee output/deepvariant.${proband_name}_trio.mendelian.txt

        # -------------------
        # STEP 6b: Docker Stop and Kill All docker Container Process
        # -------------------
        docker stop $(docker ps -aq)
        docker rm $(docker ps -aq)
    else
        echo "Merged VCF file not exist, please check the input directory"
        exit 1
    fi
fi

if [ ! -z "$proband_name" ] && [ "$solo_analysis" = true ]; then
    echo "SKIPPED STEP 5: Merge gVCF files with GLnexus, due to solo analysis (Trio)"
    echo "SKIPPED STEP 6: Calculate Mendelian Violation Rate using RTG Tools (Trio)"
fi

# ## Jump to direct to annotation session

# -------------------
# STEP 7: Annotate Variants - SnpEff SnpSift
# -------------------

if [ ! -z "$proband_name" ] && [ "$solo_analysis" = true ]; then
    echo "STEP 7a: Annotate Variants Proband - SNPEff with latest database, SnpSift ClinVar, SnpSift dbNSFP"
    
    if [ ! -f "${OUTPUT_DIR}/${proband_name}_proband.vcf.gz" ]; then
        echo "VCF file not exist, please check the input directory"
        exit 1
    fi

    # If output file not exist, run annotate
    if [ ! -f "${proband_name}-converted-deepVariant.vcf" ] && [ ! -f "${proband_name}-converted-deepVariant.vcf.gz" ]; then
        # Create GATK > dbnsfp Chromosome
        bcftools annotate \
            --rename-chrs "${chr_rename}" \
            --threads "$(nproc)" -Oz \
            -o "${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf" \
            "${OUTPUT_DIR}/${proband_name}_proband.vcf.gz"
    else
        echo "Deepvariant converted VCF file exists, skipping conversion"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf.gz" ] && [ -f "${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf" ]; then
        # SnpEff with recent GRCh38.p14 database
        SnpEff -v "${snpEff_ver}" -dataDir "${snpEff_dir}" \
            -s "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.html" \
            "${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf"
    else
        echo "SnpEff VCF file exists, skipping SnpEff"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf.gz" ] && [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf" ]; then
        # SnpSift Annotate ClinVar
        SnpSift annotate -v "${dbSNP}" \
            "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf"
    else
        echo "Annotated dbSNP file exists, skipping annotation"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf.gz" ] && [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf" ]; then
        SnpSift annotate -v "${ClinVar}" \
            "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf"
    else
        echo "Annotated ClinVar file exists, skipping annotation"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf.gz" ] && [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf" ]; then
        # Annotate using dbNSFP for SNP Only (Indel Give 0 Annotation Result)
        SnpSift dbnsfp -v -db "${dbnsfp}" \
            "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf"
    else
        echo "Annotated dbNSFP file exists, skipping annotation"
    fi

    # Compress file size using bgzip if the previous .vcf file exist
    if [ -f "${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf" ] && [ ! -f "${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf"
    fi
fi

# Similar to above, but for trio_analysis

if [ ! -z "$proband_name" ] && [ "$trio_analysis" = true ]; then
    echo "STEP 7b: Annotate Variants Trio - SNP Sift"

    # If output file not exist yet, and previous file exist
    if [ ! -f "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf" ] && [ ! -f "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf.gz" ]; then
        # Create GATK > dbnsfp Chromosome
        bcftools annotate \
            --rename-chrs "${chr_rename}" \
            --threads "$(nproc)" -Oz \
            -o "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf" \
            "${OUTPUT_DIR}/${proband_name}_trio_annotated.output.vcf.gz"
    else
        echo "DeepTrio converted VCF file exists, skipping conversion"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf.gz" ] && [ -f "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf" ]; then
        # SnpEff with recent GRCh38.p14 database
        SnpEff -v "${snpEff_ver}" -dataDir "${snpEff_dir}" \
            -s "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.html" \
            "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf"
    else
        echo "SnpEff VCF file exists, skipping SnpEff"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf.gz" ] && [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf" ]; then
        # SnpSift
        SnpSift annotate -v "${dbSNP}" \
            "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf"
    else
        echo "Annotated dbSNP file exists, skipping annotation"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf.gz" ] && [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf" ]; then
        # SnpSift Annotate ClinVar
        SnpSift annotate -v "${ClinVar}" \
            "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf"
    else
        echo "Annotated ClinVar file exists, skipping annotation"
    fi

    # If output file not exist yet, and previous file exist
    if [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf.gz" ] && [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf" ]; then
        # Annotate using dbNSFP for SNP Only (Indel Give 0 Annotation Result)
        SnpSift dbnsfp -v -db "${dbnsfp}" \
            "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf" \
            > "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf"
    else
        echo "Annotated dbNSFP file exists, skipping annotation"
    fi

    # Compress file size using bgzip if the previous .vcf file exist
    if [ -f "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf" ] && [ ! -f "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf"
    fi

    if [ -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf" ] && [ ! -f "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf.gz" ]; then
        bgzip --threads "$(nproc)" "${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf"
    fi
fi

# Additional SV Step using Tiddit @ SV_DIR

# Check if proband_name is not null and annotated tiddit output is not yet exist in folder
if [ ! -z "$proband_name" ] && [ ! -f "${SV_DIR}/output.filtered.dbnsfp.vcf" ]; then
    echo "STEP 8: SV calling using Tiddit, SnpEff annotation, SnpSift dbNSFP annotation"

    if [ -f "${INPUT_DIR}/${proband_name}.bam" ]; then
        echo "Final BAM file exist, continue to Tiddit"
        docker run \
            -v "${INPUT_DIR}":"/input"   \
            -v "${SV_DIR}":"/output"  \
            -v "${REFERENCE_DIR}":"/reference" \
            quay.io/biocontainers/tiddit:${tiddit_version} tiddit \
            --sv \
            --ref /reference/${REFSEQ_FASTA} \
            --bam /input/${proband_name}.bam \
            --skip_assembly \
            --threads $(nproc) \
            -o /output/output

        # Stop Docker
        docker stop $(docker ps -aq)
        docker rm $(docker ps -aq)

        grep -E "#|PASS" ${SV_DIR}/output.vcf > ${SV_DIR}/output.filtered.vcf

        SnpEff -v ${snpEff_ver} -dataDir ${snpEff_dir} \
            ${SV_DIR}/output.filtered.vcf > ${SV_DIR}/output.filtered.snpeff.vcf

        SnpSift dbnsfp -v -db ${dbnsfp} \
            ${SV_DIR}/output.filtered.snpeff.vcf > ${SV_DIR}/output.filtered.dbnsfp.vcf
    else
        echo "BAM file not exist, please check the input directory"
        exit 1
    fi
else
    echo "SV output folder exist, and there are already final data, skipping Tiddit"  
fi

# Step using Exomiser @ Exomiser root folder

# Check if solo_analysis is true
if [ "$solo_analysis" = true ]; then
    echo "STEP 9: Exomiser SNP Solo Analysis"
    exomiser --analysis ${exomiser_solo}
fi

# Check if trio_analysis is true
if [ "$trio_analysis" = true ]; then
    echo "STEP 9: Exomiser SNP Trio Analysis"
    exomiser --analysis ${exomiser_trio}
fi

# Check if SV VCF is true
if [ -f "${SV_DIR}/output.filtered.dbnsfp.vcf" ]; then
    echo "STEP 9: Exomiser SV Solo Analysis"
    exomiser --analysis ${exomiser_solo_sv}
fi