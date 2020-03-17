import os.path
import pandas as pd
from collections import defaultdict
import pickle


domains_path = '/media/vineetb/t5-vineetb/dsprint/out/pfam/32/hmms'
canonic_prot_path = '/media/vineetb/t5-vineetb/dsprint/in/pfam/32/domains_canonic_prot'
OUTPUT_PATH = '/media/vineetb/t5-vineetb/dsprint/out/pfam/32/'
INSTANCE_THRESHOLD = 10


if __name__ == '__main__':

    with open(os.path.join(OUTPUT_PATH, f'filtered{INSTANCE_THRESHOLD}_list.pik'), 'rb') as f:
        filtered_domains_list = pickle.load(f)
    filtered_domains_list.sort()

    gene_dict = defaultdict(dict)
    genes_list = []

    for domain_name in filtered_domains_list:
        print(domain_name)

        filename = domain_name + ".csv"
        domain_data = pd.read_csv(os.path.join(domains_path, filename), sep='\t', index_col=0, dtype={"chrom_num": str})

        # Sort the domain data
        sorted_domain_data = domain_data.sort_values(by=["chrom_num", "gene", "TargetStart"])
        sorted_domain_data = sorted_domain_data.reset_index(drop=True)

        with open(os.path.join(canonic_prot_path, f'{domain_name}_canonic_prot.pik'), 'rb') as f:
            canonic_protein = pickle.load(f)

        for gene in sorted_domain_data.loc[:, 'gene']:
            # No need to process a gene twice
            if gene in gene_dict[domain_name]:
                continue
            # Get sequence
            protein = canonic_protein[gene]
            seq = sorted_domain_data.loc[sorted_domain_data.loc[:, 'prot'] == protein, 'Target_Seq'].values[0]
            gene_dict[domain_name][gene] = seq.replace('-', '').replace('X', '').replace('.', ' ').upper()

            genes_list.append(gene)

    with open(os.path.join(OUTPUT_PATH, 'domains_sequences_dict.pik'), 'wb') as f:
        pickle.dump(gene_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
