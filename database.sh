#!/bin/bash
set -euo pipefail

# This script is used to download all the application dependencies used in the paper to the selected Path, make sure you have already create the path
CURRENT_PATH=$(pwd)
TARGET_PATH="$HOME/Downloads/Database"

# Ensure the directory exists
mkdir -p "$TARGET_PATH"

# Change to the target directory
cd "$TARGET_PATH"

echo "Starting Download All Neccessary Database, ensure stable connection"

echo "01. Downloading Reference Sequence GRCh38.p14 from UCSC"
## Reference Sequence GRCh38.p13
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p14/hg38.p14.fa.gz -O hg38.p14.fa.gz
gzip -d hg38.p14.fa.gz

echo "Run indexing with bwa-mem2, ensure your Memory + Swap > 76GB for smoother experiences"
bwa-mem2 index hg38.p14.fa

echo "Create index using samtools"
samtools faidx hg38.p14.fa

echo "02. Download and process dbNSFP from https://sites.google.com/site/jpopgen/dbNSFP?pli=1"
## Paper used dbNSFP4.4a, This link download dbNSFP4.5a, see full release update on https://sites.google.com/site/jpopgen/dbNSFP/changelog?authuser=0
dbNSFP_version=4.5a

echo "02a. Downloading dbNSFP${dbNSFP_version}"
wget https://usf.box.com/shared/static/2hzcx5s6p1xui7oen16xqzndfrkt8l9l -O dbNSFP${dbNSFP_version}.zip
unzip dbNSFP${dbNSFP_version}.zip

echo "02c. Processing the file into a single file version"

### Create Header
zcat dbNSFP${dbNSFP_version}_variant.chr1.gz | head -n 1 | bgzip > header.gz

### NOTE: bgzip parameter --threads X represents number of threads
cat dbNSFP${dbNSFP_version}_variant.chr{1..22}.gz dbNSFP${dbNSFP_version}_variant.chrX.gz dbNSFP${dbNSFP_version}_variant.chrY.gz dbNSFP${dbNSFP_version}_variant.chrM.gz | zgrep -v '#chr' | bgzip --threads nproc > dbNSFPv${dbNSFP_version}_custom.gz

### Add header back into file
cat header.gz dbNSFPv${dbNSFP_version}_custom.gz > dbNSFPv${dbNSFP_version}_custombuild.gz

zcat dbNSFPv${dbNSFP_version}_custombuild.gz | awk '!/^#/{print $0}' | cut -f 1 | sort | uniq | wc -l

### Create tabix index
tabix -s 1 -b 2 -e 2 dbNSFPv${dbNSFP_version}_custombuild.gz


echo "03. Download and Processing dbSNP from NCBI FTP"

echo "03a. Downloading dbSNP b156"
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz -O GCF_000001405.40.gz
wget https://ftp.ncbi.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz.tbi -O GCF_000001405.40.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt -O GCF_000001405.40_GRCh38.p14_assembly_report.txt

echo "03b. Processing the file into SnpSift ready file"
bcftools annotate \
  --rename-chrs ${CURRENT_PATH}/convert/GCF_000001405.40_GRCh38.p14_assembly_report_revised_snpsift.chrnames \
  --threads nproc -Oz \
  -o Homo_sapiens_assembly38.dbsnp156_snpsift.vcf.gz \
  GCF_000001405.40.gz


echo "03c. Check and filter out forbidden charset annotation"
## Reference : https://pcingola.github.io/SnpEff/ss_faq/
SnpSift vcfCheck Homo_sapiens_assembly38.dbsnp156_snpsift.vcf.gz 2>&1 | grep "INFO field" | cut -f 2 -d "'" | sort | uniq -c

zcat Homo_sapiens_assembly38.dbsnp156_snpsift.vcf.gz \
    | sed 's/\&base_change=/\&base_change%3D/g' \
    | sed 's/A=;/A%3D;/' \
    | sed 's/C=;/C%3D;/' \
    | sed 's/G=;/G%3D;/' \
    | sed 's/T=;/T%3D;/' \
    | sed 's/=,;/%3D,/' \
    > Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf

# Let's also compress and index the new file so we can use it as a database
bgzip Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf 
tabix Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz 

SnpSift vcfCheck Homo_sapiens_assembly38.dbsnp156_snpsift.fixed.vcf.gz 2>&1 | grep "INFO field" | cut -f 2 -d "'" | sort | uniq -c


echo "04. Download and Processing ClinVar from NCBI FTP"

echo "04a. Downloading ClinVar data"
# Paper used 20230514, while the most recent one is 20240215
# File version of 20230514 : 
## https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2023/clinvar_20230514.vcf.gz
## https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2023/clinvar_20230514_papu.vcf.gz.tbi
clinvar_version=20240215

wget https://ftp.ncbi.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20240215.vcf.gz -O clinvar_${clinvar_version}.vcf.gz
tabix clinvar_${clinvar_version}.vcf.gz

echo "04b. SnpSift check"
SnpSift vcfCheck clinvar_${clinvar_version}.vcf.gz 2>&1 | grep "INFO field" | cut -f 2 -d "'" | sort | uniq -c