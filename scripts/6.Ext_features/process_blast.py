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
    if len(sys.argv) != 6:
        print('Usage: <script> <hmm_file> <domain_seq_dict> <blast_path> <blast_db> <output_folder>')
        sys.exit(0)

    HMM, DOMAIN_SEQUENCES_DICT, BLAST_PATH, BLAST_DB, OUTPUT_FOLDER = sys.argv[1:]
else:
    HMM = snakemake.params.hmm
    DOMAIN_SEQUENCES_DICT = snakemake.input.domain_sequences_dict
    OUTPUT_FOLDER = snakemake.output.output_folder

    with open(os.path.join(os.path.dirname(dsprint.__file__), 'config.json'), 'r') as f:
        config = json.load(f)
        BLAST_PATH = config['blast']['path']
        BLAST_DB = config['blast']['db']

if __name__ == '__main__':

    os.mkdir(OUTPUT_FOLDER)
    with open(DOMAIN_SEQUENCES_DICT, 'rb') as f:
        genes = pickle.load(f)[HMM]

    for gene, seq in genes.items():
        with NamedTemporaryFile(mode='w') as in_file:
            in_file.write(seq)
            in_file.flush()

            out_file_name = os.path.join(OUTPUT_FOLDER, f'{gene}.pssm')

            cline = NcbipsiblastCommandline(
                cmd=os.path.join(BLAST_PATH, 'psiblast'),
                query=in_file.name,
                db=BLAST_DB,
                num_iterations=3,
                num_alignments=1,
                num_threads=8,
                out_ascii_pssm=out_file_name
            )
            print(cline)
            stdout, stderr = cline()
            print('done')
