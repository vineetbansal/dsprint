import os
import pandas as pd
import pickle


def count_domain_instances(domain_gene_table, count_overlaps=True):
    if count_overlaps:
        return domain_gene_table.shape[0]

    else:
        instance_counter = 0
        last_target_end = 0

        for i, row in domain_gene_table.iterrows():
            curr_target_start = int(row["TargetStart"])
            curr_target_end = int(row["TargetEnd"])
            if curr_target_start > last_target_end:
                instance_counter += 1
                last_target_end = curr_target_end
            # If the instance overlpas the previous one
            else:
                # Updating to the smaller traget end
                if curr_target_end > last_target_end:
                    last_target_end = curr_target_end
                # Continue without incrememnting the counter
                continue

        return instance_counter


domains_path = '/media/vineetb/t5-vineetb/dsprint/out/pfam/32/hmms'
canonic_prot_path = '/media/vineetb/t5-vineetb/dsprint/in/pfam/32/domains_canonic_prot'
OUTPUT_PATH = '/media/vineetb/t5-vineetb/dsprint/out/pfam/32/'
INSTANCE_THRESHOLD = 10

if __name__ == '__main__':

    # domains_stats = {}
    #
    # for dom_filename in os.listdir(domains_path):
    #     print(dom_filename)
    #     curr_domain_stats = []
    #
    #     domain_sym = dom_filename[:dom_filename.find(".")]
    #     domain_data = pd.read_csv(os.path.join(domains_path, dom_filename), sep='\t', index_col=0)
    #     domain_ens_genes = (domain_data["gene"]).unique()
    #
    #     with open(os.path.join(canonic_prot_path, domain_sym + '_canonic_prot.pik'), 'rb') as handle:
    #         canonic_protein = pickle.load(handle)
    #
    #     domain_instance_num = 0
    #     for ens_gene in domain_ens_genes:
    #         # Filtering the domain data for this gene according to the canonical protein id
    #         canonic_prot = canonic_protein[ens_gene]
    #         canonic_prot_t = canonic_prot[:canonic_prot.find(".")]  # Trimming the ".#" at the end
    #         domain_gene_table = domain_data[domain_data["prot"] == canonic_prot]
    #
    #         # Count the number of domain instances in this gene
    #         domain_gene_table = domain_gene_table.sort_values(by=["TargetStart", "BitScore"], ascending=[True, False])
    #         gene_intance_num = count_domain_instances(domain_gene_table, count_overlaps=True)
    #         domain_instance_num += gene_intance_num
    #
    #     # Saving domains stats:
    #     curr_domain_stats.append(len(domain_ens_genes))  # No. of genes
    #     curr_domain_stats.append(domain_instance_num)  # No. of domain instances
    #
    #     # Updating the big dict
    #     domains_stats[domain_sym] = curr_domain_stats
    #
    # with open(os.path.join(OUTPUT_PATH, 'domains_stats_dict.pik'), 'wb') as f:
    #     pickle.dump(domains_stats, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(os.path.join(OUTPUT_PATH, 'domains_stats_dict.pik'), 'rb') as f:
        domains_stats = pickle.load(f)

    domains_stats_df = pd.DataFrame.from_dict(domains_stats, orient='index')
    domains_stats_df.columns = ["genes", "instances"]
    domains_stats_df = domains_stats_df.sort_values(by=["instances", "genes"], ascending=[False, False])
    domains_stats_df.to_csv(os.path.join(OUTPUT_PATH, 'domains_stats_df.csv'), sep='\t')

    domains_stats_df = pd.read_csv(os.path.join(OUTPUT_PATH, 'domains_stats_df.csv'), sep='\t', index_col=0)

    all_domains_list = domains_stats_df.index
    with open(os.path.join(OUTPUT_PATH, 'all_domains_list.pik'), 'wb') as handle:
        pickle.dump(all_domains_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    filtered_domains_df = domains_stats_df[domains_stats_df["instances"] > INSTANCE_THRESHOLD]

    filtered_domains_df.to_csv(os.path.join(OUTPUT_PATH, f'filtered{INSTANCE_THRESHOLD}_domains_df.csv'), sep='\t')

    filtered_domains_list = []
    for domain in domains_stats.keys():
        if domains_stats[domain][1] > INSTANCE_THRESHOLD:
            filtered_domains_list.append(domain)

    with open(os.path.join(OUTPUT_PATH, f'filtered{INSTANCE_THRESHOLD}_list.pik'), 'wb') as f:
        pickle.dump(filtered_domains_list, f, protocol=pickle.HIGHEST_PROTOCOL)
