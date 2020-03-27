configfile: "dsprint/config.json"
threads: 4

PFAM_VERSIONS = ['32']
CHROMOSOMES = list(range(1, 23)) + ['X', 'Y']

wildcard_constraints:
    pfam_version="\d+"

rule all:
    input:
        expand(f"{config['output_dir']}/csq/{config['csq']['output_prefix']}{{chromosome}}{config['csq']['output_suffix']}", chromosome=CHROMOSOMES),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/all_domains_genes_prot_seq.pik", pfam_version=PFAM_VERSIONS),
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_sequences_dict.pik", pfam_version=PFAM_VERSIONS),

        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/updated_domain_to_clan_dict.pik", pfam_version=PFAM_VERSIONS),      # domain -> clan mapping
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/updated_clan_to_domains_dict.pik", pfam_version=PFAM_VERSIONS),     # clan -> domains mapping
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/updated_domain_to_pfam_acc_dict.pik", pfam_version=PFAM_VERSIONS),  # domain -> pfam accession id mapping
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domain_to_clan_dict.pik", pfam_version=PFAM_VERSIONS),              # domain -> clan mapping, for human proteome
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/clan_to_domains_dict.pik", pfam_version=PFAM_VERSIONS),             # clan -> domains mapping, for human proteome
        expand(f"{config['output_dir']}/pfam/{{pfam_version}}/domain_to_pfam_acc_dict.pik", pfam_version=PFAM_VERSIONS)           # domain -> pfam accession id mapping, for human proteome

# -----------------------------------------------------------------------------
# Parse human chromosome data and save useful information in a csv file,
# one csv file per chromosome
# -----------------------------------------------------------------------------
rule csq:
    input: f"{config['input_dir']}/{config['exac_file']}"
    output: f"{config['output_dir']}/csq/{config['csq']['output_prefix']}{{chromosome}}{config['csq']['output_suffix']}"
    script: "scripts/1.parse_ExAC/ExAC_parser.py"

# -----------------------------------------------------------------------------
# Parse pfam data and save useful information (domain_name, length,
# gathering threshold) in a csv file
# -----------------------------------------------------------------------------
rule parse_pfam:
    input: f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.hmm"
    output: f"{config['output_dir']}/pfam/{{pfam_version}}/pfam.csv"
    script: "scripts/2.parse_Pfam/parse_pfam.py"

# -----------------------------------------------------------------------------
# Parse a 'clans' tsv and save useful objects
# -----------------------------------------------------------------------------
rule handle_clans:
    input:
        f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.clans.tsv",
        f"{config['input_dir']}/pfam/{{pfam_version}}/9606.tsv"
    output:
        f"{config['output_dir']}/pfam/{{pfam_version}}/updated_domain_to_clan_dict.pik",      # domain -> clan mapping
        f"{config['output_dir']}/pfam/{{pfam_version}}/updated_clan_to_domains_dict.pik",     # clan -> domains mapping
        f"{config['output_dir']}/pfam/{{pfam_version}}/updated_domain_to_pfam_acc_dict.pik",  # domain -> pfam accession id mapping
        f"{config['output_dir']}/pfam/{{pfam_version}}/domain_to_clan_dict.pik",              # domain -> clan mapping, for human proteome
        f"{config['output_dir']}/pfam/{{pfam_version}}/clan_to_domains_dict.pik",             # clan -> domains mapping, for human proteome
        f"{config['output_dir']}/pfam/{{pfam_version}}/domain_to_pfam_acc_dict.pik"           # domain -> pfam accession id mapping, for human proteome
    script: "scripts/2.parse_Pfam/map_domain_to_clan.py"

# -----------------------------------------------------------------------------
# For a Pfam database, save a mapping
#   <domain_name>: [<log_prob1>, <log_prob2>, .. ] for all transition states
# -----------------------------------------------------------------------------
rule emission_prob:
    input: f"{config['input_dir']}/pfam/{{pfam_version}}/Pfam-A.hmm"
    output:
        f"{config['output_dir']}/pfam/{{pfam_version}}/domains_hmm_dict.pik",
        f"{config['output_dir']}/pfam/{{pfam_version}}/domains_hmm_prob_dict.pik"
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
    input: f"{config['input_dir']}/pfam/{{pfam_version}}/allhmmresbyprot-v{{pfam_version}}.tsv"
    output: f"{config['output_dir']}/pfam/{{pfam_version}}/allhmm_parsed-v{{pfam_version}}.csv"
    script: "scripts/3.parse_HMMER/process_hmmer_results.py"

# -----------------------------------------------------------------------------
# csv files, one per domain
# after filtering domain data to the domain instances that contain the major allele of 'conserved' states
# with emission probability above 0.99
# -----------------------------------------------------------------------------
rule get_domain_hmm:
    input:
        f"{config['output_dir']}/pfam/{{pfam_version}}/allhmm_parsed-v{{pfam_version}}.csv",
        f"{config['output_dir']}/pfam/{{pfam_version}}/pfam.csv",
        f"{config['output_dir']}/pfam/{{pfam_version}}/domains_hmm_prob_dict.pik"
    output: directory(f"{config['output_dir']}/pfam/{{pfam_version}}/hmms")
    script: "scripts/3.parse_HMMER/get_domain_hmm.py"

# -----------------------------------------------------------------------------
# For each domain, for each gene in the domain, find the canonical protein id
# and save as a dictionary <gene_id>: <protein_id> in the file
#   <domain>_canonic_prot.pik
# -----------------------------------------------------------------------------
rule canonical_protein:
    input:
        f"{config['output_dir']}/pfam/{{pfam_version}}/hmms",
        f"{config['input_dir']}/{config['uniprot']['fasta']}",
        f"{config['input_dir']}/{config['uniprot']['idmapping']}"
    output:
        directory(f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot")
    script: "scripts/4.parse_Uniprot/canonical_protein.py"

# -----------------------------------------------------------------------------
# For every gene and its canonical protein id, find out the amino acid sequence
# and save to a dictionary <gene_id>: { <canon_protein_id>: 'MGSRAEL..'}
# -----------------------------------------------------------------------------
rule canonic_prot_seq:
    params:
        hmm_folder=f"{config['output_dir']}/pfam/{{pfam_version}}/hmms"
    input:
        canonic_prot_folder=f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot",
        exons_seqs_folder=f"{config['input_dir']}/exons_seqs",
        hg19_file=f"{config['input_dir']}/hg19.2bit",
        exon_len_file=f"{config['output_dir']}/exons_index_length.pik"
    output: dynamic(f"{config['output_dir']}/pfam/{{pfam_version}}/all_domains_genes_prot_seq.pik")
    script: "scripts/3.parse_HMMER/get_canonic_prot_seq.py"

# -----------------------------------------------------------------------------
# For each domain, find out how many genes, and how many instances of the
# canonical protein id exist, and save in a table - domains_stats_df.csv
# The dataframe is indexed on domain name
# Note: not generating the following files
#   human_domains_list.pik - all values taken by {hmm}
#   domains_stats_dict.pik <domain_name>: (<no_of_genes>, <no_of_instances>) (same info as the df we save here)
#   all_domains_list.pik = index values in our df
#   filtered{INSTANCE_THRESHOLD}_domains_df.csv -> filtered df where 'instances' > INSTANCE_THRESHOLD
#   filtered{INSTANCE_THRESHOLD}_list.pik -> index values of above
#   regular_human_domains_list.pik -> as above, applicable for pfam32, with some threshold TBD
# -----------------------------------------------------------------------------
rule domain_statistics:
    input:
        hmm_folder=f"{config['output_dir']}/pfam/{{pfam_version}}/hmms",
        canonic_prot_folder=f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot"
    output:
        f"{config['output_dir']}/pfam/{{pfam_version}}/domains_stats_df.csv"
    script: "scripts/5.domain_stats/domain_statistics.py"

# -----------------------------------------------------------------------------
# <domain_name>: {<gene>: <target_seq_of_canonic_protein_of_gene>, .. }
# Note: Target_Seq is transformed as seq.replace('-', '').replace('X', '').replace('.', ' ').upper()
# -----------------------------------------------------------------------------
rule domain_sequences:
    params:
        hmm_folder=f"{config['output_dir']}/pfam/{{pfam_version}}/hmms",
        canonic_prot_folder=f"{config['output_dir']}/pfam/{{pfam_version}}/domains_canonic_prot"
    input:
        f"{config['output_dir']}/pfam/{{pfam_version}}/domains_stats_df.csv",
    output:
        f"{config['output_dir']}/pfam/{{pfam_version}}/domains_sequences_dict.pik"
    script: "scripts/5.domain_stats/domains_sequences_todict.py"
