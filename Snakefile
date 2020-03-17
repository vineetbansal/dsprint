import pandas as pd

configfile: "dsprint/config.json"
threads: 1

PFAM_VERSIONS = ['32']
CHROMOSOMES = list(range(23)) + ['X', 'Y']

rule all:
    input: []

rule gunzip_exac:
    input: f"{config['input_dir']}/{config['exac_file']}"
    output: f"{config['output_dir']}/{config['exac_file'][:-3]}"
    shell: "gunzip -c --keep {input} > {output}"

# Parse human chromosome data and save useful information in a csv file
rule csq:
    input: f"{config['output_dir']}/{config['exac_file'][:-3]}"
    output: expand(f"{config['output_dir']}/csq/{config['csq']['output_prefix']}{{chromosome}}{config['csq']['output_suffix']}", chromosome=CHROMOSOMES)
    script: "scripts/1.parse_ExAC/ExAC_parser.py"

# Parse pfam data and save useful information in a csv file
rule parse_pfam:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.hmm", pfam_version=PFAM_VERSIONS)
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/pfam.csv", pfam_version=PFAM_VERSIONS)
    script: "scripts/2.parse_Pfam/parse_pfam.py"

rule handle_clans:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/{{filename}}.tsv", pfam_version=PFAM_VERSIONS, filename=['Pfam-A.clans', '9606'])
    output:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/{{pik_name}}.pik",
            pfam_version=PFAM_VERSIONS,
            pik_name=[
                'updated_domain_to_clan_dict',
                'updated_clan_to_domains_dict',
                'updated_domain_to_pfam_acc_dict',
                'domain_to_clan_dict',
                'clan_to_domains_dict',
                'domain_to_pfam_acc_dict',
            ]
        )
    script: "scripts/2.parse_Pfam/map_domain_to_clan.py"

rule emission_prob:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.hmm", pfam_version=PFAM_VERSIONS)
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/{{pik_name}}.pik", pfam_version=PFAM_VERSIONS, pik_name=['domains_hmm_dict', 'domains_hmm_prob_dict'])
    script: "scripts/2.parse_Pfam/domains_emission_prob.py"

rule exon_frameshifts:
    input: f"{config['input_dir']}/exons_seqs"
    output: f"{config['output_dir']}/exons_index_length.pik"
    script: "scripts/3.parse_HMMER/exons_frameshifts.py"

rule process_hmmer_results:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/allhmmresbyprot-v{{pfam_version}}.tsv", pfam_version=PFAM_VERSIONS)
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/allhmm_parsed-v{{pfam_version}}.csv", pfam_version=PFAM_VERSIONS)
    script: "scripts/3.parse_HMMER/process_hmmer_results.py"

rule _get_domain_hmm:
    input:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/allhmm_parsed-v{{pfam_version}}.csv", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/pfam.csv", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_hmm_prob_dict.pik", pfam_version=PFAM_VERSIONS)
    output: dynamic(expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms/{{hmm}}.csv", pfam_version=PFAM_VERSIONS, hmm='{hmm}'))
    script: "scripts/3.parse_HMMER/get_domain_hmm.py"

rule get_domain_hmm:
    input: dynamic(expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms/{{hmm}}.csv", pfam_version=PFAM_VERSIONS, hmm='{hmm}'))

rule canonic_prot_seq:
    input:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot", pfam_version=PFAM_VERSIONS),
        f"{config['input_dir']}/exons_seqs",
        f"{config['input_dir']}/hg19.2bit",
        f"{config['output_dir']}/exons_index_length.pik"
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/all_domains_genes_prot_seq.pik", pfam_version=PFAM_VERSIONS)
    script: "scripts/3.parse_HMMER/get_canonic_prot_seq.py"