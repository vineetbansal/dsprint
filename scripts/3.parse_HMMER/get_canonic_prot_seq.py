import pandas as pd
import os.path
from collections import defaultdict
import pickle
import glob
import subprocess

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from dsprint.mapping_func import create_exon_pos_table


def retrieve_exon_seq(exon_start, exon_end, chrom, hg19_file, reverse_complement=False):
    # index conversion for twoBitToFa
    # start, 1-indexed inclusive -> 0-indexed inclusive, subtract 1
    # end, 1-indexed inclusive -> 0-indexed exclusive, unchanged
    exon_start = int(exon_start) - 1

    seq = subprocess.check_output(
        f'twoBitToFa {hg19_file} stdout -seq=chr{chrom} -start={exon_start} -end={exon_end}',
        shell=True
    )
    # Remove 1st line, whitespaces and newlines
    seq = ''.join(seq.decode('ascii').split()[1:]).upper()

    if reverse_complement:
        d = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join([d[x] for x in seq])[::-1]
    else:
        return seq


TEST_PROCCESSED_DOMAINS = False

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 7:
        print('Usage: <script> <hmmer_results_folder> <canonic_protein_folder> <exon_seq_folder> <hg19_2bit_file> <frameshift_file> <output_file>')
        sys.exit(0)

    HMMS_FOLDER, CANONIC_PROT_FOLDER, EXON_SEQ_FOLDER, HG19_FILE, FRAMESHIFT_FILE, OUTPUT_FILE = sys.argv[1:]
else:
    HMMS_FOLDER = snakemake.input[0]
    CANONIC_PROT_FOLDER = snakemake.input[1]
    EXON_SEQ_FOLDER = snakemake.input[2]
    HG19_FILE = snakemake.input[3]
    FRAMESHIFT_FILE = snakemake.input[4]
    OUTPUT_FILE = snakemake.output[0]


if __name__ == '__main__':

    if TEST_PROCCESSED_DOMAINS:
        raise NotImplementedError

    chromosome_names = [str(i) for i in range(1, 23)] + ['X', 'Y']
    gene_dict = defaultdict(dict)

    for domain_file in glob.glob(HMMS_FOLDER + '/*.csv'):

        domain = os.path.splitext(os.path.basename(domain_file))[0]
        print(domain)
        domain_data = pd.read_csv(domain_file, sep='\t', index_col=0, dtype={"chrom_num": str})

        with open(os.path.join(CANONIC_PROT_FOLDER, domain + "_canonic_prot.pik"), 'rb') as f:
            canonic_protein = pickle.load(f)

        for gene, _domain_data in domain_data.groupby('gene'):

            prot_id = canonic_protein[gene]
            if gene in gene_dict and prot_id in gene_dict[gene]:
                continue

            chrom = _domain_data['chrom_num'].unique()[0]
            if chrom not in chromosome_names:
                continue

            chromosome = _domain_data[_domain_data["#TargetID"] == prot_id]['chromosome'].unique()
            if len(chromosome) > 1:
                print(" Error: " + gene + ": more than one chromosome raw data")  # sanity check
            chromosome = chromosome[0]

            exon_table = create_exon_pos_table(chromosome, prot_id, FRAMESHIFT_FILE)
            dna_seq = ''.join([
                retrieve_exon_seq(row["start_pos"], row["end_pos"], chrom, HG19_FILE,
                                  reverse_complement=chromosome.find('complement') >= 0)
                for _, row in exon_table.iterrows()
            ])

            gene_dict[gene][prot_id] = str(Seq(dna_seq, generic_dna).translate())

    # Saving one dictionary for all the domains together
    if TEST_PROCCESSED_DOMAINS:
        raise NotImplementedError
    else:
        with open(OUTPUT_FILE, 'wb') as f:
            pickle.dump(gene_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
