import os.path
import pandas as pd
from collections import defaultdict
import pickle


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 5:
        print('Usage: <script> <domain_stats_csv> <hmms_folder> <canonic_prot_folder> <output_domain_seqs_pik>')
        sys.exit(0)

    DOMAIN_STATS_CSV, HMMS_FOLDER, CANONIC_PROT_FOLDER, OUTPUT_FILE = sys.argv[1:]
else:
    DOMAIN_STATS_CSV = snakemake.input[0]
    HMMS_FOLDER = snakemake.params.hmm_folder
    CANONIC_PROT_FOLDER = snakemake.params.canonic_prot_folder
    OUTPUT_FILE = snakemake.output[0]

INSTANCE_THRESHOLD = 10  # lower bound on instances, exclusive

if __name__ == '__main__':

    domains_stats = pd.read_csv(DOMAIN_STATS_CSV, sep='\t', index_col=0)

    gene_dict = defaultdict(dict)

    for domain_name in domains_stats[domains_stats.instances > INSTANCE_THRESHOLD].index:
        print(domain_name)

        filename = domain_name + ".csv"
        domain_data = pd.read_csv(os.path.join(HMMS_FOLDER, filename), sep='\t', index_col=0, dtype={"chrom_num": str})

        # Sort the domain data
        sorted_domain_data = domain_data.sort_values(by=["chrom_num", "gene", "TargetStart"])
        sorted_domain_data = sorted_domain_data.reset_index(drop=True)

        with open(os.path.join(CANONIC_PROT_FOLDER, f'{domain_name}_canonic_prot.pik'), 'rb') as f:
            canonic_protein = pickle.load(f)

        for gene in sorted_domain_data.loc[:, 'gene']:
            # No need to process a gene twice
            if gene in gene_dict[domain_name]:
                continue
            # Get sequence
            protein = canonic_protein[gene]
            seq = sorted_domain_data.loc[sorted_domain_data.loc[:, 'prot'] == protein, 'Target_Seq'].values[0]
            gene_dict[domain_name][gene] = seq.replace('-', '').replace('X', '').replace('.', ' ').upper()

    with open(os.path.join(OUTPUT_FILE), 'wb') as f:
        pickle.dump(gene_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
