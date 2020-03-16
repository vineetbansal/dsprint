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


def find_longest_prot(prot_ids):
    max_len = 0
    max_prot_id = ""
    for prot in prot_ids:
        prot_lens = (gene_table[gene_table["prot"] == prot]["length"]).unique()
        if len(prot_lens) > 1:
            print("Error: more than one length fir protein id: "+prot)  #Sanity check
        if prot_lens[0] > max_len:
            max_len = prot_lens[0]
            max_prot_id = prot
    return max_prot_id


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

        domain_data = pd.read_csv(dom_filename, sep='\t', index_col=0)
        canonic_protein = {}
        no_uniprot = 0
        no_canonic_len = 0

        for gene_id in domain_data['gene'].unique():
            print(gene_id)
            gene_table = domain_data[domain_data["gene"] == gene_id]
            protein_ids = gene_table["prot"].unique()

            if len(protein_ids) == 1:
                # Saving the one protein id available
                canonic_protein[gene_id] = protein_ids[0]

            # If there's more then one transcript: finding what's the canonic protein length from uniprot
            else:
                uniprot_id = get_uniprot_id(gene_id.split('.')[0])
                if uniprot_id is None:
                    canonic_protein[gene_id] = find_longest_prot(protein_ids)
                    no_uniprot += 1
                    continue

                uniprot_url = "http://togows.dbcls.jp/entry/uniprot/" + uniprot_id + ".json"
                try:
                    uniprot_json = pd.read_json(uniprot_url)
                    canonic_len = uniprot_json.aalen[0]
                except:
                    print("Cannot find canonical length for this uniprot id")
                    canonic_protein[gene_id] = find_longest_prot(protein_ids)
                    no_uniprot += 1
                    continue
                else:
                    print("Found canonical length for this uniprot id")

                found = False
                for prot in protein_ids:
                    prot_lens = gene_table[gene_table['prot'] == prot]['length'].unique()
                    if len(prot_lens) > 1:
                        print("Error: more than one length for protein id: " + prot)  # Sanity check

                    # If the length equals the canonical, this is the canonical protein id
                    if prot_lens[0] == canonic_len:
                        found = True
                        canonic_protein[gene_id] = prot
                        break

                if not found:
                    print(gene_id+": Proteins don't match the Uniprot canonic")
                    no_canonic_len += 1
                    canonic_protein[gene_id] = find_longest_prot(protein_ids)

        domain = os.path.splitext(os.path.basename(dom_filename))[0]
        with open(os.path.join(OUTPUT_FOLDER, domain + "_canonic_prot.pik"), 'wb') as handle:
            pickle.dump(canonic_protein, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print("Finished " + dom_filename)

