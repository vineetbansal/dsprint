import pandas as pd

configfile: "dsprint/config.json"
threads: 1

PFAM_VERSIONS = ['32']
CHROMOSOMES = list(range(1, 23)) + ['X', 'Y']

rule all:
    input: []

rule gunzip_exac:
    input: f"{config['input_dir']}/{config['exac_file']}"
    output: f"{config['output_dir']}/{config['exac_file'][:-3]}"
    shell: "gunzip -c --keep {input} > {output}"

# -----------------------------------------------------------------------------
# Parse human chromosome data and save useful information in a csv file,
# one csv file per chromosome
# -----------------------------------------------------------------------------
rule csq:
    input: f"{config['output_dir']}/{config['exac_file'][:-3]}"
    output: expand(f"{config['output_dir']}/csq/{config['csq']['output_prefix']}{{chromosome}}{config['csq']['output_suffix']}", chromosome=CHROMOSOMES)
    script: "scripts/1.parse_ExAC/ExAC_parser.py"

# -----------------------------------------------------------------------------
# Parse pfam data and save useful information (domain_name, length,
# gathering threshold) in a csv file
# -----------------------------------------------------------------------------
rule parse_pfam:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.hmm", pfam_version=PFAM_VERSIONS)
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/pfam.csv", pfam_version=PFAM_VERSIONS)
    script: "scripts/2.parse_Pfam/parse_pfam.py"

# -----------------------------------------------------------------------------
# Parse a 'clans' tsv and save useful objects
# -----------------------------------------------------------------------------
rule handle_clans:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/{{filename}}.tsv", pfam_version=PFAM_VERSIONS, filename=['Pfam-A.clans', '9606'])
    output:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/{{pik_name}}.pik",
            pfam_version=PFAM_VERSIONS,
            pik_name=[
                'updated_domain_to_clan_dict',       # domain -> clan mapping
                'updated_clan_to_domains_dict',      # clan -> domains mapping
                'updated_domain_to_pfam_acc_dict',   # domain -> pfam accession id mapping
                'domain_to_clan_dict',               # domain -> clan mapping, for human proteome
                'clan_to_domains_dict',              # clan -> domains mapping, for human proteome
                'domain_to_pfam_acc_dict',           # domain -> pfam accession id mapping, for human proteome
            ]
        )
    script: "scripts/2.parse_Pfam/map_domain_to_clan.py"

# -----------------------------------------------------------------------------
# For a Pfam database, save a mapping
#   <domain_name>: [<log_prob1>, <log_prob2>, .. ] for all transition states
# -----------------------------------------------------------------------------
rule emission_prob:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.hmm", pfam_version=PFAM_VERSIONS)
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/{{pik_name}}.pik", pfam_version=PFAM_VERSIONS, pik_name=['domains_hmm_dict', 'domains_hmm_prob_dict'])
    script: "scripts/2.parse_Pfam/domains_emission_prob.py"

# -----------------------------------------------------------------------------
# Save an object mapping
#    <exon_id_file>: [(<position>, <base_pairs_length>, <base_pairs>), (..), ..]
# The values indicate the positions at which 'frame shifts' occur
# -----------------------------------------------------------------------------
rule exon_frameshifts:
    input: f"{config['input_dir']}/exons_seqs"
    output: f"{config['output_dir']}/exons_index_length.pik"
    script: "scripts/3.parse_HMMER/exons_frameshifts.py"

# -----------------------------------------------------------------------------
# Take as input domains (from Hmmer 2.3.2 and 3.1.b2) identified for human protein sequences
# and save in a csv file, with one row per chromosome ('chrom_num')
# The column 'chromosome' has format:
#    GRCh37:4:complement(join(68619532..68620053,68610286..68610505,68606198..68606442))
# -----------------------------------------------------------------------------
rule process_hmmer_results:
    input: expand(f"{config['input_dir']}/pfam/{{pfam_version}}/allhmmresbyprot-v{{pfam_version}}.tsv", pfam_version=PFAM_VERSIONS)
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/allhmm_parsed-v{{pfam_version}}.csv", pfam_version=PFAM_VERSIONS)
    script: "scripts/3.parse_HMMER/process_hmmer_results.py"

# -----------------------------------------------------------------------------
# private, used by get_domain_hmm
# -----------------------------------------------------------------------------
rule _get_domain_hmm:
    input:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/allhmm_parsed-v{{pfam_version}}.csv", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/pfam.csv", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_hmm_prob_dict.pik", pfam_version=PFAM_VERSIONS)
    output: dynamic(expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms/{{hmm}}.csv", pfam_version=PFAM_VERSIONS, hmm='{hmm}'))
    script: "scripts/3.parse_HMMER/get_domain_hmm.py"

# -----------------------------------------------------------------------------
# csv files, one per domain
# after filtering domain data to the domain instances that contain the major allele of 'conserved' states
# with emission probability above 0.99
# -----------------------------------------------------------------------------
rule get_domain_hmm:
    input: dynamic(expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms/{{hmm}}.csv", pfam_version=PFAM_VERSIONS, hmm='{hmm}'))

# -----------------------------------------------------------------------------
# For each domain, for each gene in the domain, find the canonical protein id
# and save as a dictionary <gene_id>: <protein_id> in the file
#   <domain>_canonic_prot.pik
# -----------------------------------------------------------------------------
rule canonical_protein:
    input: dynamic(expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms/{{hmm}}.csv", pfam_version=PFAM_VERSIONS, hmm='{hmm}'))
    output: dynamic(expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot/{{hmm}}_canonic_prot.pik", pfam_version=PFAM_VERSIONS, hmm='{hmm}'))

# -----------------------------------------------------------------------------
# For every gene and its canonical protein id, find out the amino acid sequence
# and save to a dictionary <gene_id>: { <canon_protein_id>: 'MGSRAEL..'}
# -----------------------------------------------------------------------------
rule canonic_prot_seq:
    input:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot", pfam_version=PFAM_VERSIONS),
        f"{config['input_dir']}/exons_seqs",
        f"{config['input_dir']}/hg19.2bit",
        f"{config['output_dir']}/exons_index_length.pik"
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/all_domains_genes_prot_seq.pik", pfam_version=PFAM_VERSIONS)
    script: "scripts/3.parse_HMMER/get_canonic_prot_seq.py"

# -----------------------------------------------------------------------------
# For each domain, find out how many genes, and how many instances of the
# canonical protein id exist, and save in a table - domains_stats_df.csv
# Note: not generating the following files
#   human_domains_list.pik - all values taken by {hmm}
#   domains_stats_dict.pik <domain_name>: (<no_of_genes>, <no_of_instances>) (same info as the df we save here)
#   all_domains_list.pik = index values in our df
#   filtered{INSTANCE_THRESHOLD}_domains_df.csv -> filtered df where 'instances' > INSTANCE_THRESHOLD
#   filtered{INSTANCE_THRESHOLD}_list.pik -> index values of above
# -----------------------------------------------------------------------------
rule domain_statistics:
    input:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot", pfam_version=PFAM_VERSIONS),
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_stats_df.csv", pfam_version=PFAM_VERSIONS)
    script: "scripts/5.domain_stats/domain_statistics.py"

# -----------------------------------------------------------------------------
# <domain_name>: {<gene>: <target_seq_of_canonic_protein_of_gene>, .. }
# Note: Target_Seq is transformed as seq.replace('-', '').replace('X', '').replace('.', ' ').upper()
# -----------------------------------------------------------------------------
rule domain_sequences:
    input:
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_stats_df.csv", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot", pfam_version=PFAM_VERSIONS),
    output: expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_sequences_dict.pik", pfam_version=PFAM_VERSIONS)
    script: "scripts/5.domain_stats/domains_sequences_todict.py"