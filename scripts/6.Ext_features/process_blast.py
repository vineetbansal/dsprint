import os
import pickle
import json
from tempfile import NamedTemporaryFile
from Bio.Blast.Applications import NcbipsiblastCommandline

import dsprint

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <domain_seq_dict> <output_folder>')
        sys.exit(0)

    DOMAIN_SEQUENCES_DICT, OUTPUT_FOLDER = sys.argv[1:]
else:
    DOMAIN_SEQUENCES_DICT = snakemake.input.domain_sequences_dict
    OUTPUT_FOLDER = snakemake.output.output_folder

with open(os.path.join(os.path.dirname(dsprint.__file__), 'config.json'), 'r') as f:
    config = json.load(f)['blast']
    BLAST_PATH = config['path']
    REMOTE = config['remote']

if __name__ == '__main__':

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    with open(DOMAIN_SEQUENCES_DICT, 'rb') as f:
        domain_sequences_dict = pickle.load(f)

    for domain in domain_sequences_dict:

        n_genes = len(domain_sequences_dict[domain])
        for i, (gene, seq) in enumerate(domain_sequences_dict[domain].items(), start=1):
            with NamedTemporaryFile(mode='w') as in_file:
                in_file.write(seq)
                in_file.flush()

                out_file_name = os.path.join(OUTPUT_FOLDER, f'{gene}.pssm')

                if REMOTE:
                    db = config['default_db']
                    kwargs = {}
                else:
                    db = config['dbs'][config['default_db']]
                    kwargs = {'num_iterations': config['num_iterations'], 'num_threads': config['num_threads']}

                cline = NcbipsiblastCommandline(
                    cmd=os.path.join(BLAST_PATH, 'psiblast'),
                    query=in_file.name,
                    db=db,
                    num_alignments=1,
                    out_ascii_pssm=out_file_name,
                    remote=REMOTE,
                    **kwargs
                )
                print(cline)
                stdout, stderr = cline()

            print(f'Finished gene {i}/{n_genes}')
        print(f'Finished domain {domain}')