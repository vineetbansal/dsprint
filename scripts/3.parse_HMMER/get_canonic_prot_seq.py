import pandas as pd
import os.path
import json
import numpy as np
import math
from collections import defaultdict
import pickle
import sys
import subprocess

import dsprint
from dsprint.calc_exac_freq_func import codon_table, retrieve_codon_seq
from dsprint.mapping_func import create_exon_pos_table

curr_dir = '.'
pfam_version = "32"
domains_th = "10"
TEST_PROCCESSED_DOMAINS = False


def reverse_complement(seq):
    """
    Given a DNA sequence, return the reverse-complement of that sequence.
    """

    # Complement strand - transversing the bp to base-complement
    complement_seq = []
    for c in seq:
        if c.upper() == 'A':
            complement_seq.append('T')
        elif c.upper() == 'T':
            complement_seq.append('A')
        elif c.upper() == 'G':
            complement_seq.append('C')
        else:
            complement_seq.append('G')

    # Reversing the sequence
    comp_seq = ''.join(complement_seq)
    rev_comp_seq = comp_seq[::-1]

    return rev_comp_seq


def retrieve_exon_seq(exon_start, exon_end, chrom):
    """
    Retrieve the exon sequence from the ref sequence, according to exons start and end positions
    """
    chromsome_name = "chr" + chrom
    seq_start = int(exon_start) - 1

    # Calling hg19.2bit to retreive the DNA sequence
    query = subprocess.check_output(
        "../5.HMM_alter_align/twoBitToFa ../5.HMM_alter_align/hg19.2bit stdout -seq=%s -start=%s -end=%s" % (
        chromsome_name, str(seq_start), str(exon_end)), shell=True)
    query = ''.join(query.split()[1:])  # Remove 1st line, whitespaces and newlines

    return query.upper()


def exons_translate_to_prot(exon_table, chrom_raw_data, chrom):
    dna_seq = ""

    # Get all the exons dna sequence
    for index, exon in exon_table.iterrows():
        exon_seq = retrieve_exon_seq(exon["start_pos"], exon["end_pos"], chrom)

        if chrom_raw_data.find("complement") >= 0:
            exon_seq = reverse_complement(exon_seq)

        dna_seq = dna_seq + exon_seq

    # Translate to protein sequence
    prot_seq = []
    next_codon_idx = 0
    while (next_codon_idx + 2 < len(dna_seq)):
        codon = dna_seq[next_codon_idx:next_codon_idx + 3]
        prot_seq.append(codon_table[codon])
        next_codon_idx += 3

    # Convert all codons to one amino acids string
    protein_str = ''.join(prot_seq)
    return protein_str


if __name__ == '__main__':

    with open(os.path.join(os.path.dirname(dsprint.__file__), 'config.json')) as json_file:
        config = json.load(json_file)

    # Read the list of domains
    if TEST_PROCCESSED_DOMAINS:
        with open(curr_dir[0] + "/../13.Process_domains_not_in_training/processed_domains_not_in_pipeline_final_list.pik",
                  'rb') as handle:
            filtered_domains_list = pickle.load(handle)
    else:
        if pfam_version == "32":
            with open(config['pfam'][pfam_version]['human_domains_list'], 'rb') as handle:
                filtered_domains_list = pickle.load(handle)
        else:
            with open(curr_dir[0] + "/../5.domains_stats/pfam-v" + pfam_version + "/filtered" + domains_th + "_list.pik",
                      'rb') as handle:
                filtered_domains_list = pickle.load(handle)
    filtered_domains_list.sort()

    # Read the substitutions table (for the exons translation)
    # with open("/home/anat/Research/ExAC/9.Features_exploration/codon_ns_table.pik", 'rb') as handle:
    #     codon_ns_table = pickle.load(handle)

    chromosome_names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                        "19", "20", "21", "22", "X", "Y"]
    len(filtered_domains_list)

    gene_dict = defaultdict(dict)

    for domain_name in filtered_domains_list:

        in_path = curr_dir[0] + "/hmm_domains/pfam-v" + pfam_version + "/"
        filename = domain_name + ".csv"
        domain_data = pd.read_csv(in_path + filename, sep='\t', index_col=0, dtype={"chrom_num": str})

        # Sort the domain data
        sorted_domain_data = domain_data.sort_values(by=["chrom_num", "gene", "TargetStart"])
        sorted_domain_data = sorted_domain_data.reset_index(drop=True)

        # Get the canonic protein ids file for the domain
        with open(curr_dir[
                      0] + "/../4.parse_Uniprot/domains_canonic_prot/pfam-v" + pfam_version + "/" + domain_name + "_canonic_prot.pik",
                  'rb') as handle:
            canonic_protein = pickle.load(handle)

        for gene in sorted_domain_data["gene"].unique():

            prot_id = canonic_protein[gene]
            if gene in gene_dict and prot_id in gene_dict[gene]:
                continue

            # Get the exons sequence file for the protein
            chrom = sorted_domain_data[sorted_domain_data["gene"] == gene]["chrom_num"].unique()[0]
            if chrom not in chromosome_names:
                continue

            exons_file = pd.read_csv(
                curr_dir[0] + "/from_shilpa/exons_seqs/" + chrom + "/" + gene + "/" + prot_id + ".exons.txt", skiprows=1,
                header=None, names=["pos", "exon_seq"], sep='\t')

            # Get the chrom raw data
            chrom_raw_data = \
            sorted_domain_data[sorted_domain_data["gene"] == gene][sorted_domain_data["#TargetID"] == prot_id][
                "chromosome"].unique()  # there should be only one element here
            if len(chrom_raw_data) > 1:
                print(" Error: " + gene + ": more than one chromosome raw data")  # sanity check
            chrom_raw_data = chrom_raw_data[0]

            # Create exons table
            exon_table = create_exon_pos_table(chrom_raw_data, prot_id)  # exons frameshifts are also fixed here!

            # Translate all the exons dna sequences to one protein sequence
            prot_seq = exons_translate_to_prot(exon_table, chrom_raw_data, chrom)
            gene_dict[gene][prot_id] = prot_seq

        print("Finished domain " + str(domain_name))

    # Saving one dictionary for all the domains together
    if TEST_PROCCESSED_DOMAINS:
        with open(curr_dir[0] + "/canonic_prot_seq/pfam-v" + pfam_version + "/processed_domains_genes_prot_seq.pik", 'wb') as handle:
            pickle.dump(gene_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(curr_dir[0] + "/canonic_prot_seq/pfam-v" + pfam_version + "/all_domains_genes_prot_seq.pik",
                  'wb') as handle:
            pickle.dump(gene_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

