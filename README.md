# IDeRare

IDeRare or *"Indonesia Exome Rare Disease Variant Discovery Pipeline"* is simple and ready to use variant discovery pipeline to discover rare disease variants from exome sequencing data.

## Authored by 
Ivan William Harsono<sup>a</sup>, Yulia Ariani<sup>b</sup>, Beben Benyamin<sup>c,d,e</sup>, Fadilah Fadilah<sup>f,g</sup>, Dwi Ari Pujianto<sup>b</sup>, Cut Nurul Hafifah<sup>h</sup>

<sup>a</sup>Doctoral Program in Biomedical Sciences, Faculty of Medicine, Universitas Indonesia, Jakarta, Indonesia.<br> 
<sup>b</sup>Department of Medical Biology, Faculty of Medicine, Universitas Indonesia, Jakarta, Indonesia.<br> 
<sup>c</sup>Australian Centre for Precision Health, University of South Australia, Adelaide, SA, 5000, Australia. <br>
<sup>d</sup>UniSA Allied Health and Human Performance, University of South Australia, Adelaide, SA, 5000, Australia. <br>
<sup>e</sup>South Australian Health and Medical Research Institute (SAHMRI), University of South Australia, Adelaide, SA, 5000, Australia. <br>
<sup>f</sup>Department of Medical Chemistry, Faculty of Medicine, Universitas Indonesia, Jalan Salemba Raya number 4, Jakarta, 10430, Indonesia.<br>
<sup>g</sup>Bioinformatics Core Facilities - IMERI, Faculty of Medicine, Universitas Indonesia, Jalan Salemba Raya number 6, Jakarta, 10430, Indonesia .<br>
<sup>h</sup>Department of Child Health, Dr. Cipto Mangunkusumo Hospital, Faculty of Medicine, University of Indonesia, Jakarta, Indonesia. <br>


**Note:** Currently IDeRare paper is being considered journal submission. The citation will be updated once the paper is published.

## Description
- This pipeline is designed to be used in Linux environment
- Original paper may used different version of tools, and the prerequisite used in this pipeline is the latest version of the tools
- This pipeline is designed and tested with Indonesia rare disease trio patient, but it should be also usable for general cases of rare disease variant discovery from Exome Sequences data given paired end .fq.gz file and HPO data(s)
- Ensure you have at least 250GB free for database and application setup, and 100GB free for each Trio family exome set
- The .yaml file path are assuming all the folder are stored in ```Downloads``` folder with subfolder of ```Database``` (for RefSeq, dbNSFP, dbSNP, ClinVar), ```Sandbox``` (for application and its database), ```IDeRare``` (git cloned folder)

## Quick Install
1. Clone this repository
```bash
git clone https://github.com/ivanwilliammd/IDeRare
```
2. Have a Linux environment (Ubuntu or Ubuntu-like 22.04 LTS distro is recommended)
3. Install [Docker](https://docs.docker.com/engine/install/ubuntu/) and [Anaconda - optional](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)- see [Prerequisite.md](installation/Prerequisite.md) for more details
4. Run dependency installation script and database script
```bash
cd installation
source install_dependencies.sh
source download_database.sh
cd ../
```
5. **Optional** : translate the coded clinical information (mixture of SNOMED-CT, Orphanet, HPO) separated by enter at ```clinical_data.txt``` file. <br><br>
This will update the ```iderare.yml``` *hpo_ids* section with the HPO code, and show the non-matched SNOMED-CT code in terminal.
```bash
conda activate iderare
python iderare_phentl.py
```
<br> Example of ```clinical_data.txt``` file
```txt
SNOMEDCT:65959000 // Differential Diagnosis 1 : Beta thalassemia (disorder)
SNOMEDCT:11179002 // Differential Diagnosis 2 : Glycogen storage disease, type IV (disorder)
SNOMEDCT:258211005 // Inheritance pattern : Autosomal recessive inheritance (qualifier value)
SNOMEDCT:165397008 // Laboratory finding : Hemoglobin low (finding)
SNOMEDCT:80515008 // Physical examination finding : Hepatomegaly
OMIM:230800 // Gaucher Disease
ORPHA:848 // Glycogen storage disease, type IV
HP:0001396 // Cholestasis
HP:0002240 //  Hepatomegaly
LOINC:2862-1|L // Hypoalbuminemia
```
1. Set the data, directory file reference and trio information on ```iderare.yml```.<br><br> 
**Note** : all exome files should be located in the ```input/A_FASTQ``` folder of absolute path setup by ```data_dir``` at ```iderare.yml```
<br><br>
![File Structure](picture/01.png)
![Example File](picture/02.png)
1. Run the bash script
```bash
source iderare.sh
```