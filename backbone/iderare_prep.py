#!/usr/bin/python

import os
import yaml
import argparse

# Define a custom type function to validate the input
def valid_proband_type(proband_type):
    if proband_type not in ['solo', 'trio', 'both']:
        raise argparse.ArgumentTypeError("Invalid proband type. Allowed values are 'solo', 'trio', or 'both'")
    return proband_type

parser = argparse.ArgumentParser("iderare_prep.py", description="Prepare the pipeline.sh, trio.ped, and exomiser YAML files for the iDeRare pipeline")
parser.add_argument("-m", "--mode", help="Mode of analysis ('solo', 'trio', 'both')", type=valid_proband_type, default='both')
parser.add_argument("-t", "--trimming", help="Command to directly add the final HPOset to iderare.yml file", action="store_true", default=False)  # Default value set to False
args = parser.parse_args()

# Check if trimming is provided as an argument
if args.trimming:
    trimming_value = True
else:
    trimming_value = False

# Check if mode of analysis is provided as an argument
if args.mode:
    mode_value = args.mode
else:
    mode_value = 'both'

with open("iderare.yml") as i:
    y = yaml.safe_load(i)

    # Load data for trio analysis
    data_dir = y['analysis']['data_dir'] if y['analysis']['data_dir'] else ''
    proband = y['analysis']['proband'] if y['analysis']['proband'] else ''
    mother = y['analysis']['mother'] if y['analysis']['mother'] else ''
    father = y['analysis']['father'] if y['analysis']['father'] else ''

    proband_gender = y['analysis']['proband_gender'] if y['analysis']['proband_gender'] else ''
    proband_phen = y['analysis']['proband_phen'] if y['analysis']['proband_phen'] else ''
    mother_phen = y['analysis']['mother_phen'] if y['analysis']['mother_phen'] else ''
    father_phen = y['analysis']['father_phen'] if y['analysis']['mother_phen'] else ''
    hpo_ids = y['analysis']['hpo_ids'] if y['analysis']['hpo_ids'] else ''

    library = y['analysis']['library'] if y['analysis']['library'] else ''
    method = y['analysis']['method'] if y['analysis']['method'] else ''

    # Load data for setup database path / source
    dv_version = y['setup']['dv_version'] if y['setup']['dv_version'] else ''
    dv_model = y['setup']['dv_model'] if y['setup']['dv_model'] else ''

    glnexus_version = y['setup']['glnexus_version'] if y['setup']['glnexus_version'] else ''
    tiddit_version = y['setup']['tiddit_version'] if y['setup']['tiddit_version'] else ''
    max_mem = y['setup']['max_mem'] if y['setup']['max_mem'] else ''

    ref_dir = y['setup']['ref_dir'] if y['setup']['ref_dir'] else ''
    ref_fasta = y['setup']['ref_fasta'] if y['setup']['ref_fasta'] else ''

    snpEff_dir = y['setup']['snpEff_dir'] if y['setup']['snpEff_dir'] else ''
    snpEff_ver = y['setup']['snpEff_ver'] if y['setup']['snpEff_ver'] else ''

    exomiser_dir = y['setup']['exomiser_dir'] if y['setup']['exomiser_dir'] else ''
    exomiser_data_ver = y['setup']['exomiser_data_ver'] if y['setup']['exomiser_data_ver'] else ''
    
    dbNSFP_file = y['setup']['dbNSFP_file'] if y['setup']['dbNSFP_file'] else ''
    dbSNP_file = y['setup']['dbSNP_file'] if y['setup']['dbSNP_file'] else ''
    ClinVar_file = y['setup']['ClinVar_file'] if y['setup']['ClinVar_file'] else ''

    chr_rename = y['setup']['chr_rename'] if y['setup']['chr_rename'] else ''

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

    # Check if trimming is provided as an argument
    if trimming_value:
        template = template.replace("{{trimming}}", "true")
    else:
        template = template.replace("{{trimming}}", "false")

    # Check if mode is provided as an argument
    if mode_value == 'solo':
        template = template.replace("{{trio_analysis}}", "false")
        template = template.replace("{{solo_analysis}}", "true")
    elif mode_value == 'trio':
        template = template.replace("{{trio_analysis}}", "true")
        template = template.replace("{{solo_analysis}}", "false")
    else :
        template = template.replace("{{trio_analysis}}", "true")
        template = template.replace("{{solo_analysis}}", "true")

    # Export the new pipeline.yml
    with open("pipeline.sh", "w+") as p:
        p.write(template)
        print("pipeline.sh created")

# If trio analysis is selected, create trio.ped
if args.mode == 'trio' or args.mode == 'both':
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
else :
    print("trio.ped not created")

# Create Trio and Solo Exomiser YAML
## Trio
if args.mode == 'trio' or args.mode == 'both':
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

if args.mode == 'solo' or args.mode == 'both':
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

## Tiddit (Solo only)
if args.mode == 'solo' or args.mode == 'both':
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