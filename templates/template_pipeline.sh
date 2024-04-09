#!/bin/bash
set -euo pipefail

echo "STEP 0 : Setup Variable, Directory, and Ensure NVIDIA Docker Supported"
### WARNING : Don't use tilde (~) in the path

# -------------------
## Step 0a: Set variable necessary with sample name for RG
# -------------------
echo "STEP 0a: Set variable necessary"

BIN_VERSION={{dv_version}}
memory={{max_mem}}
DEEPTRIO_MODEL={{dv_model}}

glnexus_version={{glnexus_version}}
tiddit_version={{tiddit_version}}

# Father
father_name={{father}}

# Mother
mother_name={{mother}}

# Proband
proband_name={{proband}}

proband_SM=Proband
proband_PU={{library}}
proband_PL={{method}}
proband_LB=WXS

exomiser_solo={{proband}}_exomiser_solo.yml
exomiser_solo_sv={{proband}}_tiddit_exomiser_solo.yml
exomiser_trio={{proband}}_exomiser_trio.yml

mother_SM=Mother
mother_PU={{library}}
mother_PL={{method}}
mother_LB=WXS

father_SM=Father
father_PU={{library}}
father_PL={{method}}
father_LB=WXS

trimming={{trimming}}
solo_analysis={{solo_analysis}}
trio_analysis={{trio_analysis}}

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

# -------------------
## Step 0b: Prepare directory and pedigree file
# -------------------
echo "STEP 0b: Prepare directory and pedigree file"

## Pedigree file on INPUT_DIR (Just file name only)
INPUT_DIR="{{data_dir}}/input"
PEDIGREE="trio.ped"

## Subfolder of INPUT_DIR
FASTQ_DIR=${INPUT_DIR}/A_FASTQ
SAM_DIR=${INPUT_DIR}/B_RAW_SAM_BAM

## Output Dir
OUTPUT_DIR="{{data_dir}}/output"
ANNOTATED_DIR="{{data_dir}}/annotated"
SV_DIR="{{data_dir}}/sv_tiddit"
EXOMISER_DIR="{{data_dir}}/exomiser"

## FASTA file on REFERENCE_DIR (Just file name only)
REFERENCE_DIR={{ref_dir}}
REFSEQ_FASTA={{ref_fasta}}

snpEff_dir={{snpEff_dir}}
snpEff_ver={{snpEff_ver}}
dbnsfp={{dbNSFP_file}}
dbSNP={{dbSNP_file}}
ClinVar={{ClinVar_file}}
chr_rename={{chr_rename}}

cd {{data_dir}}

## Make Directory
mkdir -p ${INPUT_DIR} ${FASTQ_DIR} ${SAM_DIR} ${OUTPUT_DIR} ${OUTPUT_DIR}/intermediate_results_dir_proband ${OUTPUT_DIR}/intermediate_results_dir_trio ${ANNOTATED_DIR} ${SV_DIR} ${EXOMISER_DIR}

# -------------------
# STEP 1: QC - Run fastqp 
# -------------------
echo "STEP 1: QC - Run fastqp - Optional (if the input is raw untrimmed fastq)"
echo "STEP 1a : Proband"

# Check if proband_name is null
if [ -z "$proband_name" ]; then
    echo "proband_name is null, stopping the pipeline."
    exit 1
else
    echo "proband_name is not null, continuing alignment."
    
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

echo "STEP 1b : Mother"

# Check if mother_name is null
if [ -z "$mother_name" ]; then
    echo "mother_name is null, stopping the pipeline."
    exit 1
else
    echo "mother_name is not null, continuing alignment."
    
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

echo "STEP 1C : Father"
if [ -z "$father_name" ]; then
    echo "father_name is null, stopping the pipeline."
    exit 1
else
    echo "father_name is not null, continuing alignment."
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

# --------------------------------------
# STEP 2: Map to reference using BWA-MEM2
# --------------------------------------
echo "STEP 2: Map to reference using BWA-MEM2"

echo "STEP 2a : Proband"

# BWA MEM2 Alignment
bwa-mem2 mem -R "@RG\tID:${proband_name}\tSM:${proband_SM}\tPU:${proband_PU}\tPL:${proband_PL}\tLB:${proband_LB}" \
    -t $(nproc) ${REFERENCE_DIR}/${REFSEQ_FASTA} \
    ${FASTQ_DIR}/${proband_name}_1.fq.gz \
    ${FASTQ_DIR}/${proband_name}_2.fq.gz \
    > ${SAM_DIR}/${proband_name}_raw.sam

echo "STEP 2b : Mother"

bwa-mem2 mem -R "@RG\tID:${mother_name}\tSM:${mother_SM}\tPU:${mother_PU}\tPL:${mother_PL}\tLB:${mother_LB}" \
    -t $(nproc) ${REFERENCE_DIR}/${REFSEQ_FASTA} \
    ${FASTQ_DIR}/${mother_name}_1.fq.gz \
    ${FASTQ_DIR}/${mother_name}_2.fq.gz \
    > ${SAM_DIR}/${mother_name}_raw.sam

echo "STEP 2c : Father"

bwa-mem2 mem -R "@RG\tID:${father_name}\tSM:${father_SM}\tPU:${father_PU}\tPL:${father_PL}\tLB:${father_LB}" \
    -t $(nproc) ${REFERENCE_DIR}/${REFSEQ_FASTA} \
    ${FASTQ_DIR}/${father_name}_1.fq.gz \
    ${FASTQ_DIR}/${father_name}_2.fq.gz \
    > ${SAM_DIR}/${father_name}_raw.sam

# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - sambamba
# -----------------------------------------
echo "STEP 3: Mark Duplicates and Sort - Sambamba"

# Conversion of SAM to BAM & Markdup
echo "STEP 3a : Proband"

sambamba view -p -t=$(nproc) -l=9 \
    -S ${SAM_DIR}/${proband_name}_raw.sam \
    -f=bam -o=${SAM_DIR}/${proband_name}_raw.bam

sambamba markdup -r -p -t=$(nproc) -l=9 \
    ${SAM_DIR}/${proband_name}_raw.bam \
    ${SAM_DIR}/${proband_name}_dedup.bam

sambamba sort -m=${memory} -p -t=$(nproc) -l=9 \
    ${SAM_DIR}/${proband_name}_dedup.bam \
    -o=${INPUT_DIR}/${proband_name}.bam

echo "STEP 3b : Mother"

sambamba view -p -t=$(nproc) -l=9 \
    -S ${SAM_DIR}/${mother_name}_raw.sam \
    -f=bam -o=${SAM_DIR}/${mother_name}_raw.bam

sambamba markdup -r -p -t=$(nproc) -l=9 \
    ${SAM_DIR}/${mother_name}_raw.bam \
    ${SAM_DIR}/${mother_name}_dedup.bam

sambamba sort -m=${memory} -p -t=$(nproc) -l=9 \
    ${SAM_DIR}/${mother_name}_dedup.bam \
    -o=${INPUT_DIR}/${mother_name}.bam

echo "STEP 3c : Father"

sambamba view -p -t=$(nproc) -l=9 \
    -S ${SAM_DIR}/${father_name}_raw.sam \
    -f=bam -o=${SAM_DIR}/${father_name}_raw.bam

sambamba markdup -r -p -t=$(nproc) -l=9 \
    ${SAM_DIR}/${father_name}_raw.bam \
    ${SAM_DIR}/${father_name}_dedup.bam

sambamba sort -m=${memory} -p -t=$(nproc) -l=9 \
    ${SAM_DIR}/${father_name}_dedup.bam \
    -o=${INPUT_DIR}/${father_name}.bam

### Remove Intermediate SAM BAM file as it consumes too much spaces

#### Remove all in one folder 
rm ${SAM_DIR}/*.sam ${SAM_DIR}/*.bam

#### Remove one by one
rm ${SAM_DIR}/${proband_name}_raw.sam ${SAM_DIR}/${mother_name}_raw.sam ${SAM_DIR}/${father_name}_raw.sam
rm ${SAM_DIR}/${proband_name}_raw.bam ${SAM_DIR}/${mother_name}_raw.bam ${SAM_DIR}/${father_name}_raw.bam
rm ${SAM_DIR}/${proband_name}_dedup.bam ${SAM_DIR}/${mother_name}_dedup.bam ${SAM_DIR}/${father_name}_dedup.bam

# ----------------------------------------------
# STEP 4: Variant Calling
# ----------------------------------------------
echo "STEP 4a: Variant Calling Proband DeepVariant"

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

echo "STEP 4b: Variant Calling DeepTrio"

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

# ----------------------------------------------
# STEP 5: Merge gVCF files with GLnexus
# ----------------------------------------------
echo "STEP 5: Merge gVCF files with GLnexus"

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

# ----------------------------------------------
# STEP 6: Calculate Mendelian Violation Rate using RTG Tools
# ----------------------------------------------
echo "STEP 6: Calculate Mendelian Violation Rate using RTG Tools"
#### Reference : https://www.animalgenome.org/bioinfo/resources/manuals/RTGOperationsManual.pdf

if [ -d ${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf ];
then
    echo "${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf exists. Removing directory."
    rm -r -f ${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf
else
	echo "${REFERENCE_DIR}/${REFSEQ_FASTA}.sdf directory does not exist. Continue"
fi

docker run \
    -v "${INPUT_DIR}":"/input" \
    -v "${REFERENCE_DIR}":"/reference" \
    realtimegenomics/rtg-tools format \
    -o /reference/${REFSEQ_FASTA}.sdf "/reference/${REFSEQ_FASTA}"

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

# ## Jump to direct to annotation session

# -------------------
# STEP 7: Annotate Variants - SnpEff SnpSift
# -------------------
echo "STEP 7a: Annotate Variants Proband - SNPEff with latest database, SnpSift ClinVar, SnpSift dbNSFP"

# Create GATK > dbnsfp Chromosome
bcftools annotate \
  --rename-chrs ${chr_rename} \
  --threads nproc -Oz \
  -o ${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf \
  ${OUTPUT_DIR}/${proband_name}_proband.vcf.gz

# SnpEff with recent GRCh38.p14 database
SnpEff -v ${snpEff_ver} -dataDir ${snpEff_dir} \
    -s ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.html \
    ${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf

# SnpSift Annotate ClinVar
SnpSift annotate -v ${dbSNP} \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf

SnpSift annotate -v ${ClinVar} \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf

# Annotate using dbNSFP for SNP Only (Indel Give 0 Annotation Result)
SnpSift dbnsfp -v -db ${dbnsfp} \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf

# Final SnpEff for known vs unknown dbSNP
SnpEff -v ${snpEff_ver} -dataDir ${snpEff_dir} \
    -s ${ANNOTATED_DIR}/${proband_name}-2ndSnpEff-dbSNP-ClinVar-deepVariant.html \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-2ndSnpEff-dbSNP-ClinVar-deepVariant.vcf

# Compress file size using bgzip
bgzip --threads nproc ${OUTPUT_DIR}/${proband_name}-converted-deepVariant.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepVariant.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepVariant.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepVariant.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf


echo "STEP 7b: Annotate Variants Trio - SNP Sift"

# Create GATK > dbnsfp Chromosome
bcftools annotate \
  --rename-chrs ${chr_rename} \
  --threads nproc -Oz \
  -o ${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf \
  ${OUTPUT_DIR}/${proband_name}_trio_annotated.output.vcf.gz

# SnpEff with recent GRCh38.p14 database
SnpEff -v ${snpEff_ver} -dataDir ${snpEff_dir} \
    -s ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.html \
    ${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf

# SnpSift Annotate dbSNP & ClinVar
SnpSift annotate -v ${dbSNP} \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf

SnpSift annotate -v ${ClinVar} \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf

# Annotate using dbNSFP for SNP Only (Indel Give 0 Annotation Result)
SnpSift dbnsfp -v -db ${dbnsfp} \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf

# Final SnpEff for known vs unknown dbSNP
SnpEff -v ${snpEff_ver} -dataDir ${snpEff_dir} \
    -s ${ANNOTATED_DIR}/${proband_name}-2ndSnpEff-dbSNP-ClinVar-deepTrio.html \
    ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf \
    > ${ANNOTATED_DIR}/${proband_name}-2ndSnpEff-dbSNP-ClinVar-deepTrio.vcf

# Compress file size using bgzip
bgzip --threads nproc ${OUTPUT_DIR}/${proband_name}-converted-deepTrio.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-deepTrio.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-deepTrio.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-deepTrio.vcf
bgzip --threads nproc ${ANNOTATED_DIR}/${proband_name}-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf


# echo "STEP 8: SV calling - Tiddit"
# # Additional SV Step using Tiddit @ SV_DIR

# docker run \
#     -v "${INPUT_DIR}":"/input"   \
#     -v "${SV_DIR}":"/output"  \
#     -v "${REFERENCE_DIR}":"/reference" \
#     quay.io/biocontainers/tiddit:${tiddit_version} tiddit \
#     --sv \
#     --ref /reference/${REFSEQ_FASTA} \
#     --bam /input/${proband_name}.bam \
#     --skip_assembly \
#     --threads $(nproc) \
#     -o /output/output

# # Stop Docker
# docker stop $(docker ps -aq)
# docker rm $(docker ps -aq)

# grep -E "#|PASS" ${SV_DIR}/output.vcf > ${SV_DIR}/output.filtered.vcf

# SnpEff -v ${snpEff_ver} -dataDir ${snpEff_dir} \
#     ${SV_DIR}/output.filtered.vcf > ${SV_DIR}/output.filtered.snpeff.vcf

# SnpSift dbnsfp -v -db ${dbnsfp} \
#     ${SV_DIR}/output.filtered.snpeff.vcf > ${SV_DIR}/output.filtered.dbnsfp.vcf

echo "STEP 9: Exomiser Analysis"
# Step using Exomiser @ Exomiser root folder
exomiser --analysis ${exomiser_solo}
exomiser --analysis ${exomiser_trio}
exomiser --analysis ${exomiser_solo_sv}