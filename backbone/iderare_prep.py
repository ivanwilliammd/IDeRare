#!/usr/bin/python

import os
import yaml

with open("iderare.yml") as i:
    y = yaml.safe_load(i)

    # Load data for trio analysis
    data_dir = y['analysis']['data_dir']
    proband = y['analysis']['proband']
    mother = y['analysis']['mother']
    father = y['analysis']['father']

    proband_gender = y['analysis']['proband_gender']
    proband_phen = y['analysis']['proband_phen']
    mother_phen = y['analysis']['mother_phen']
    father_phen = y['analysis']['father_phen']
    hpo_ids = y['analysis']['hpo_ids']

    library = y['analysis']['library']
    method = y['analysis']['method']

    # Load data for setup database path / source
    dv_version = y['setup']['dv_version']
    dv_model = y['setup']['dv_model']

    glnexus_version = y['setup']['glnexus_version']
    tiddit_version = y['setup']['tiddit_version']
    max_mem = y['setup']['max_mem']

    ref_dir = y['setup']['ref_dir']
    ref_fasta = y['setup']['ref_fasta']

    snpEff_dir = y['setup']['snpEff_dir']
    snpEff_ver = y['setup']['snpEff_ver']

    exomiser_dir = y['setup']['exomiser_dir']
    exomiser_data_ver = y['setup']['exomiser_data_ver']
    
    dbNSFP_file = y['setup']['dbNSFP_file']
    dbSNP_file = y['setup']['dbSNP_file']
    ClinVar_file = y['setup']['ClinVar_file']

    chr_rename = y['setup']['chr_rename']


# Open template_pipeline.yml and replace {{variable}} with the value from iderare.yml
with open("templates/template_pipeline.sh") as t:
    template = t.read()
    template = template.replace("{{data_dir}}", data_dir)
    template = template.replace("{{proband}}", proband)
    template = template.replace("{{mother}}", mother)
    template = template.replace("{{father}}", father)

    template = template.replace("{{library}}", library)
    template = template.replace("{{method}}", method)

    template = template.replace("{{dv_version}}", dv_version)
    template = template.replace("{{dv_model}}", dv_model)
    template = template.replace("{{max_mem}}", max_mem)

    template = template.replace("{{glnexus_version}}", glnexus_version)
    template = template.replace("{{tiddit_version}}", tiddit_version)

    template = template.replace("{{ref_dir}}", ref_dir)
    template = template.replace("{{ref_fasta}}", ref_fasta)

    template = template.replace("{{snpEff_dir}}", snpEff_dir)
    template = template.replace("{{snpEff_ver}}", snpEff_ver)

    template = template.replace("{{dbNSFP_file}}", dbNSFP_file)
    template = template.replace("{{dbSNP_file}}", dbSNP_file)
    template = template.replace("{{ClinVar_file}}", ClinVar_file)

    template = template.replace("{{chr_rename}}", chr_rename)

    # Export the new pipeline.yml
    with open("pipeline.sh", "w+") as p:
        p.write(template)
        print("pipeline.sh created")


# Open template_trio.ped and replace with value from iderare.yml
with open("templates/template_trio.ped") as t:
    template = t.read()

    template = template.replace("{{proband_gender}}", str(proband_gender))
    template = template.replace("{{proband_phen}}", str(proband_phen))
    template = template.replace("{{mother_phen}}", str(mother_phen))
    template = template.replace("{{father_phen}}", str(father_phen))

    # Check if directory exist
    if not os.path.exists(os.path.join(data_dir, "input")):
        os.makedirs(os.path.join(data_dir, "input"))

    # Export the new trio.ped
    with open(os.path.join(data_dir, "input", "trio.ped"), "w+") as p:
        p.write(template)
        print("trio.ped created")

# Create Trio and Solo Exomiser YAML
## Trio
with open("templates/template_exomiser.yml") as i:
    y = yaml.safe_load(i)
    y['analysis']['vcf'] = data_dir + '/annotated/' + proband + '-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepTrio.vcf.gz'
    y['analysis']['ped'] = data_dir + "/input/trio.ped"
    y['analysis']['proband'] = 'child'
    y['analysis']['hpoIds'] = hpo_ids
    y['outputOptions']['outputDirectory'] = data_dir + '/exomiser'
    y['outputOptions']['outputFileName'] = proband + "-exomiser-trio"

    
with open(os.path.join(data_dir, proband + "_exomiser_trio.yml"), "w+") as o:
    yaml.dump(y, o, default_flow_style=False, sort_keys=False)
    print(proband + "_exomiser_trio.yml created")

## Solo
with open("templates/template_exomiser.yml") as i:
    y = yaml.safe_load(i)
    y['analysis']['vcf'] = data_dir + '/annotated/' + proband + '-SnpEff-dbSNP-ClinVar-dbNSFP_annotated-deepVariant.vcf.gz'
    y['analysis']['hpoIds'] = hpo_ids
    y['outputOptions']['outputDirectory'] = data_dir + '/exomiser'
    y['outputOptions']['outputFileName'] = proband + "-exomiser-solo"

with open(os.path.join(data_dir, proband + "_exomiser_solo.yml"), "w+") as o:
    yaml.dump(y, o, default_flow_style=False, sort_keys=False)
    print(proband + "_exomiser_solo.yml created")

## Tiddit
with open("templates/template_exomiser.yml") as i:
    y = yaml.safe_load(i)
    y['analysis']['vcf'] = data_dir + '/sv_tiddit/' + 'output.filtered.dbnsfp.vcf'
    y['analysis']['hpoIds'] = hpo_ids
    y['outputOptions']['outputDirectory'] = data_dir + '/exomiser'
    y['outputOptions']['outputFileName'] = proband + "-tiddit-exomiser-solo"

with open(os.path.join(data_dir, proband + "_tiddit_exomiser_solo.yml"), "w+") as o:
    yaml.dump(y, o, default_flow_style=False, sort_keys=False)
    print(proband + "_tiddit_exomiser_solo.yml created")

# Open template_application.properties and replace {{variable}} with the value from iderare.yml
with open("templates/template_application.properties") as t:
    template = t.read()
    template = template.replace("{{exomiser_dir}}", exomiser_dir)
    template = template.replace("{{exomiser_data_ver}}", str(exomiser_data_ver))

    # Export the new application.properties
    with open(os.path.join(data_dir, "application.properties"), "w+") as p:
        p.write(template)
        print("application.properties created")