#!/bin/bash
set -euo pipefail

# This script is used to download all the application dependencies used in the paper to the selected Path, make sure you have already create the path
TARGET_PATH="$HOME/Downloads/Sandbox"

# Parse arguments
SKIP_EXECUTABLE=false

for arg in "$@"; do
    case $arg in
        --skip-executable)
            SKIP_EXECUTABLE=true
            ;;
        *)
            ;;
    esac
done

# Download Exomiser data from https://data.monarchinitiative.org/exomiser/data/index.html & extract
exomiser_data_version=2309

# Ensure the directory exists
mkdir -p "$TARGET_PATH"

# Change to the target directory
cd "$TARGET_PATH"

# Skip steps 1-6 if requested
if [ "$SKIP_EXECUTABLE" = true ]; then
    echo "Skipping steps 1-6 as per request."
else
    echo "Download Neccesary Application"

    echo "01. Downloading fastp"
    # Download fastp from https://github.com/OpenGene/fastp/archive/refs/tags/v0.23.4.zip & extract
    wget http://opengene.org/fastp/fastp.0.23.4
    mv fastp.0.23.4 fastp
    chmod a+x ./fastp

    echo "02. Downloading bwa-mem2"
    # Download bwa-mem2 from https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 & extract
    wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 -O bwa-mem2-2.2.1_x64-linux.tar.bz2
    tar -xvf bwa-mem2-2.2.1_x64-linux.tar.bz2
    mv bwa-mem2-2.2.1_x64-linux bwa-mem2
    chmod a+x ./bwa-mem2/bwa-mem2

    echo "03. Downloading sambamba"
    # Download executable sambamba from https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz & extract
    # Used in paper is v1.0.0 conda distribution, but there is no different with the 1.0.1 version
    wget https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz -O sambamba-1.0.1-linux-amd64-static.gz
    gunzip sambamba-1.0.1-linux-amd64-static.gz
    mv sambamba-1.0.1-linux-amd64-static sambamba
    chmod a+x ./sambamba

    echo "04. Downloading bcftools"
    # Download bcftools and samtools from https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 , https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 & extract
    # Used in paper is v1.17, but there is no different with the core feature of --rename-chrs is still the same
    wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 -O bcftools-1.19.tar.bz2
    tar -xvf bcftools-1.19.tar.bz2
    # Make the bcftools executable
    cd bcftools-1.19
    make
    chmod a+x ./bcftools
    cd ../

    wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 -O samtools-1.19.2.tar.bz2
    tar -xvf samtools-1.19.2.tar.bz2
    # Make the samtools executable
    cd samtools-1.19.2
    make
    chmod a+x ./samtools
    cd ../

    echo "05. Downloading SnpEff and SnpSift"
    # Download SnpEff from https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip & extract
    wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
    chmod a+x ./snpEff/snpEff.jar ./snpEff/SnpSift.jar

    echo "05a. Downloading SnpEff Database"
    java -jar ./snpEff/snpEff.jar download GRCh38.105

    echo "06. Downloading Exomiser"
    # Download Exomiser from https://github.com/exomiser/Exomiser/releases/download/13.3.0/exomiser-cli-13.3.0-distribution.zip & extract
    wget https://github.com/exomiser/Exomiser/releases/download/13.3.0/exomiser-cli-13.3.0-distribution.zip -O exomiser-cli-13.3.0-distribution.zip
    unzip exomiser-cli-13.3.0-distribution.zip
    mv exomiser-cli-13.3.0-distribution exomiser
    chmod a+x ./exomiser/exomiser-cli-13.3.0.jar

    echo "06b. Downloading Exomiser Data"
    # Paper use 2302 database, but the latest version is 2309
    mkdir -p exomiser/data/exomiser-data
    wget https://data.monarchinitiative.org/exomiser/data/${exomiser_data_version}_hg38.zip -O exomiser/data/exomiser-data/${exomiser_data_version}_hg38.zip
    wget https://data.monarchinitiative.org/exomiser/data/${exomiser_data_version}_phenotype.zip -O exomiser/data/exomiser-data/${exomiser_data_version}_phenotype.zip
    unzip exomiser/data/exomiser-data/${exomiser_data_version}_hg38.zip -d exomiser/data/exomiser-data
    unzip exomiser/data/exomiser-data/${exomiser_data_version}_phenotype.zip -d exomiser/data/exomiser-data

echo "06c. Create application properties"
cat <<EOF > exomiser/application.properties
exomiser.data-directory=$TARGET_PATH/exomiser/data/exomiser-data
exomiser.hg38.data-version=$exomiser_data_version
exomiser.phenotype.data-version=$exomiser_data_version
EOF

fi


echo "07. Downloading Tiddit"
# Paper used Tiddit v3.6.0, but the is no different in sv feature compared to v3.6.1 (minor bug fix)
docker run quay.io/biocontainers/tiddit:3.6.1--py38h24c8ff8_0

echo "STEP 8: Installing Dependency of NVIDIA Docker 2"

# (1) Install nvidia driver:
# https://linuxhint.com/install-cuda-ubuntu/
sudo apt-get "${APT_ARGS[@]}" update
sudo apt-get "${APT_ARGS[@]}" install \
  build-essential \
  curl \
  "linux-headers-$(uname -r)" \
  nvidia-cuda-toolkit

# See https://www.tensorflow.org/install/source#gpu for versions required.
if ! dpkg-query -W cuda-11-3; then
  echo "Installing CUDA..."
  UBUNTU_VERSION="2004"
  curl -O https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
  sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
  # From https://forums.developer.nvidia.com/t/notice-cuda-linux-repository-key-rotation/212772
  sudo -H apt-key adv --fetch-keys "http://developer.download.nvidia.com/compute/cuda/repos/ubuntu${UBUNTU_VERSION}/x86_64/3bf863cc.pub"
  sudo add-apt-repository -y "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/ /"
  sudo -H apt-get update "${APT_ARGS[@]}" > /dev/null
  sudo -H apt-get full-upgrade "${APT_ARGS[@]}" > /dev/null
  sudo -H apt-get install "${APT_ARGS[@]}" cuda-11-3
fi

echo "Checking for CUDNN..."
if [[ ! -e /usr/local/cuda-11/include/cudnn.h ]]; then
  echo "Installing CUDNN..."
  CUDNN_TAR_FILE="cudnn-11.3-linux-x64-v8.2.0.53.tgz"
  wget -q https://developer.download.nvidia.com/compute/redist/cudnn/v8.2.0/${CUDNN_TAR_FILE}
  tar -xzvf ${CUDNN_TAR_FILE}
  sudo cp -P cuda/include/cudnn.h /usr/local/cuda-11/include
  sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-11/lib64/
  sudo cp -P cuda/lib64/libcudnn* /usr/local/cuda-11/lib64/
  sudo chmod a+r /usr/local/cuda-11/lib64/libcudnn*
  sudo ldconfig
fi
# (3) Install nvidia docker:
# https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html#installing-on-ubuntu-and-debian
# Add the package repositories
distribution=$(. /etc/os-release;echo "$ID$VERSION_ID")
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L "https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list" | sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update && sudo apt-get install "${APT_ARGS[@]}" nvidia-docker2
sudo systemctl restart docker

fi 

echo "STEP 8: Dry Run NVIDIA Docker 2 using GPU"
#### Test nvidia-smi with the latest official CUDA image
docker run --gpus 1 nvidia/cuda:11.3.0-base-ubuntu20.04 nvidia-smi

echo "09. Downloading DeepVariant Docker Image"
## Paper used DeepVariant v1.5.0 & GLNexus 1.2.7 
## but the latest version is DeepVariant v1.6.0 & GLNexus v1.4.1, 
## Main difference is that DeepVariant 1.6.0 occupy around 16GB of memory (Keras), 
## While DeepVariant 1.5.0 used around 6GB of memory and trained on Slim
## See this link for details : https://github.com/google/deepvariant/blob/r1.6/docs/FAQ.md 
## Uncomment the 1.6.0 version if you have more than 16GB of GPU memory
# docker pull google/deepvariant:1.6.0-gpu
# docker pull google/deepvariant:deeptrio-1.6.0-gpu
docker pull google/deepvariant:1.5.0-gpu
docker pull google/deepvariant:deeptrio-1.5.0-gpu
docker pull realtimegenomics/rtg-tools
docker pull ghcr.io/dnanexus-rnd/glnexus:v1.4.1

echo "Final : Add all the application to the .bashrc"

# Add aliases to ~/.bashrc\
echo "" >> ~/.bashrc
echo "alias fastp=\"$TARGET_PATH/fastp\";" >> ~/.bashrc
echo "alias bwa-mem2=\"$TARGET_PATH/bwa-mem2/bwa-mem2\";" >> ~/.bashrc
echo "alias sambamba=\"$TARGET_PATH/sambamba\";" >> ~/.bashrc
echo "alias bcftools=\"$TARGET_PATH/bcftools-1.19/bcftools\";" >> ~/.bashrc
echo "alias samtools=\"$TARGET_PATH/samtools-1.19.2/samtools\";" >> ~/.bashrc
echo "alias SnpEff=\"java -jar $TARGET_PATH/snpEff/snpEff.jar\";" >> ~/.bashrc
echo "alias SnpSift=\"java -jar $TARGET_PATH/snpEff/SnpSift.jar\";" >> ~/.bashrc
echo "alias exomiser=\"java -jar $TARGET_PATH/exomiser/exomiser-cli-13.3.0.jar\";" >> ~/.bashrc

# Source ~/.bashrc to apply changes
source ~/.bashrc
echo "Aliases added to ~/.bashrc"

# Install Python Dependencies Required
echo "Removing old iderare environment"
conda remove --name iderare --all -y
conda create --name iderare python=3.8 -y
conda activate iderare
echo "Install Python Dependencies"
conda install pandas pyyaml numpy scipy matplotlib -y
pip install iderare-pheno
conda deactivate iderare
