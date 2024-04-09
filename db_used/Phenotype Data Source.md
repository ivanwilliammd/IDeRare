# Phenotype and Orphanet Data Source

### Requirements
1. HPO to SNOMED map : ```hp.obo``` from [HPO Ontology download](https://hpo.jax.org/app/data/ontology)
2. HPO complete ontology data from [HPO annotation download](https://hpo.jax.org/app/data/annotations)
3. SNOMED CT to ORPHA map : See further instruction below
4. ORPHA to OMIM map : ```ORDO_en_4.4.owl``` from [ORDO download link at Bioportal](https://bioportal.bioontology.org/ontologies/ORDO)
5. LOINC to HPO map : ```loinc2hpo-annotations.tsv``` from [LOINC 2 HPO annotation](https://github.com/TheJacksonLaboratory/loinc2hpoAnnotation.git)


### Build HPO to SNOMED map & ORPHANET to HPO Map Subset
1. Snapshot downloaded from HPO website at 2024-03-10 attached at ```rawdl_20240310```.
2. Run ```Generate HPO Subset.ipynb``` to generate the subset.
3. Resulting file saved at ```subset``` folder named ```hpo2snomed_subset.tsv``` and ```orpha2hpo_subset.tsv```. *Note : 2 is to*

### Build SNOMED CT to ORPHANET Subset
1. Download 2 (two) zip : SNOMED-CT RF2 & SNOMED-CT Orphanet Map from [MLDS](https://mlds.ihtsdotools.org/#/landing) (accessible from each countries SNOMED CT National Release Center)
2. Copy the following file Navigate though the following file : 
   - Concept : ```Snapshot/Terminology/sct2_Concept_Snapshot_INT_*.txt``` file.
   - Description : ```sct2_Description_Snapshot-en_INT_*.txt``` file.
   - Language Reference Set : ```sct2_sRefset_LanguageSnapshot-en_INT_*.txt``` file.
   - Replace * with the release version you downloaded, replace INT with the language code of the release (if you didn't use reference).
3. Run ```Generate SNOMED ORPHA Subset.ipynb``` to generate the subset.
4. Resulting file between ```20240301 International SNOMED-CT RF2``` and ```20230731 Orphanet Map``` saved at ```subset``` folder named ```snomed2orpha_subset.tsv```.