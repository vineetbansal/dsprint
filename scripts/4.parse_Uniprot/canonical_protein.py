import pandas as pd
import glob
import os.path
from biomart import BiomartServer
import pickle

ENS_GENES = BiomartServer("http://grch37.ensembl.org/biomart").datasets[u'hsapiens_gene_ensembl']


def get_uniprot_id(ensembl_id):
    d = {
        'filters': {'ensembl_gene_id': ensembl_id},
        'attributes': ['uniprotswissprot']
    }

    if ENS_GENES.count(d) == 0:
        return None

    response = ENS_GENES.search(d)
    for line in response.iter_lines():
        line = line.decode('utf-8')
        return line.split("\t")[0]


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <hmmer_results_folder>')
        sys.exit(0)

    HMMS_FOLDER, OUTPUT_FOLDER = sys.argv[1:]
else:
    HMMS_FOLDER = snakemake.input[0]
    OUTPUT_FOLDER = snakemake.output[0]


if __name__ == '__main__':

    domains_files = glob.glob(HMMS_FOLDER + '/*.csv')

    # TODO: Is the 6393 in the notebook an experiment?
    for dom_filename in domains_files:

        domain = os.path.splitext(os.path.basename(dom_filename))[0]
        domain_data = pd.read_csv(dom_filename, sep='\t', index_col=0)
        canonic_protein = {}

        for gene_id, gene_table in domain_data.groupby('gene'):
            protein_ids = gene_table['prot'].unique()

            if len(protein_ids) == 1:
                canonic_protein[gene_id] = protein_ids[0]

            # If there's more then one transcript: finding what's the canonic protein length from Uniprot
            else:
                uniprot_id = get_uniprot_id(gene_id.split('.')[0])
                if uniprot_id is None:  # Uniprot ID unavailable? Just take the longest protein
                    canonic_protein[gene_id] = gene_table.loc[gene_table.length.idxmax()].prot
                else:
                    try:
                        canonic_len = pd.read_json(f'http://togows.dbcls.jp/entry/uniprot/{uniprot_id}.json').aalen[0]
                    except:  # Unable to determine canonical length? Just take the longest protein
                        canonic_protein[gene_id] = gene_table.loc[gene_table.length.idxmax()].prot
                    else:
                        _gene_table = gene_table[gene_table.length == canonic_len]
                        if not _gene_table.empty:
                            canonic_protein[gene_id] = _gene_table.prot.iloc[0]
                        else:  # No protein matching canonical length? Just take the longest protein
                            canonic_protein[gene_id] = gene_table.loc[gene_table.length.idxmax()].prot

        with open(os.path.join(OUTPUT_FOLDER, domain + '_canonic_prot.pik'), 'wb') as f:
            pickle.dump(canonic_protein, f, protocol=pickle.HIGHEST_PROTOCOL)
