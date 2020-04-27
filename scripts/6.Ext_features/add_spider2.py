import sys
import os.path
import glob
import pandas as pd
import pickle
from functools import lru_cache

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 7:
        print('Usage: <script> <domain> <hmm_folder> <hmm_states_folder> <canonic_prot_folder> <spider_folder> <all_domains_genes_prot_seq_file>')
        sys.exit(0)

    DOMAIN, HMM_FOLDER, HMM_STATES_FOLDER, CANONIC_PROT_FOLDER, SPIDER_FOLDER, ALL_DOMAINS_GENE_PROT_SEQ_FILE = sys.argv[1:]
else:
    DOMAIN, HMM_FOLDER, HMM_STATES_FOLDER, CANONIC_PROT_FOLDER, SPIDER_FOLDER, ALL_DOMAINS_GENE_PROT_SEQ_FILE = snakemake.input


# file extensions by spider that have csv data
SPIDER_EXTS = 'spd3', 'hsa2', 'hsb2'
# columns in the above found files that give us protein positions, and can thus serve as DataFrame indices
SPIDER_INDEX_COLS = '#', '#index', '#index'


@lru_cache(maxsize=8192)
def spider_dataframes(domain, gene):
    folder = os.path.join(SPIDER_FOLDER, domain)
    files = [os.path.join(folder, f'{gene}.{ext}') for ext in SPIDER_EXTS]
    dfs = [None] * len(files)
    for i, file in enumerate(files):
        if os.path.exists(file):
            dfs[i] = pd.read_csv(file, sep='\t', index_col=SPIDER_INDEX_COLS[i])

    merged = dfs[0].merge(
        dfs[1], left_index=True, right_index=True
    ).merge(
        dfs[2], left_index=True, right_index=True, suffixes=('_hsa2', '_hsb2')
    )

    return merged


if __name__ == '__main__':

    with open(os.path.join(CANONIC_PROT_FOLDER, f'{DOMAIN}_canonic_prot.pik'), 'rb') as f:
        canonic_protein = pickle.load(f)

    with open(ALL_DOMAINS_GENE_PROT_SEQ_FILE, 'rb') as f:
        seqs = pickle.load(f)

    with open(f'{HMM_STATES_FOLDER}/{DOMAIN}.pik', 'rb') as f:
        states_dict = pickle.load(f)

    for state, ds in states_dict.items():

        for d in ds:
            chromosome = d['chrom']
            gene = d['ens_gene']
            protein_position = d['prot_pos']
            states_dict_aa = d['aa_ref_orig']

            protein = canonic_protein[gene]

            seq = seqs[gene][protein]
            seq_ = seq[0:d['prot_pos']]
            skip = seq_.count('X') + seq_.count('*') + seq_.count('.') + seq_.count('-')
            if skip > 0:
                raise RuntimeError('not ready for this shit')
            spider_pos = d['prot_pos'] - skip

            spider_df = spider_dataframes(DOMAIN, gene)

            print('debug')