import pandas as pd
import os.path
from collections import defaultdict
import pickle

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from calc_exac_freq_func import create_alt_codon, exac_validation_checks, retrieve_codon_seq, codon_table
from af_format_calc import format_af, calculate_af_adj

from dsprint.core import CHROMOSOMES, get_chromosome_number, retrieve_exon_seq
from dsprint.mapping_func import is_number, find_chrom_bps

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 6:
        print('Usage: <script> <hmm_file> <canonic_prot_folder> <indels_folder> <hg19_file> <output_file>')
        sys.exit(0)

    HMM_FILE, CANONIC_PROT_FOLDER, INDELS_FOLDER, HG19_FILE, OUTPUT_FILE = sys.argv[1:]
else:
    HMM_FILE, CANONIC_PROT_FOLDER, INDELS_FOLDER, HG19_FILE = snakemake.input
    OUTPUT_FILE = str(snakemake.output)


def change_ref_aa(res_dict, alterations_af_dict, alterations_af_adj_dict, aa, aa_sum, aa_adj_sum, bp_af_dict, bp_af_adj_dict):
    """Changing the ref aa and updating all relevant dictionaries"""

    # Adding the refrence allele to the alterations dicts
    old_ref = res_dict["aa_ref"]
    sum_of_all_alt = sum(sum(alterations_af_dict.values(), []))
    sum_of_all_alt_adj = sum(sum(alterations_af_adj_dict.values(), []))
    alterations_af_dict[old_ref] = [1 - sum_of_all_alt]
    alterations_af_adj_dict[old_ref] = [1 - sum_of_all_alt_adj]

    # Updating the aa to be the ref
    res_dict["aa_ref"] = aa

    # Finding the codon that codes aa with highest frequency and update bp_ref
    old_bp_ref = res_dict["bp_ref"]
    max_af = 0
    for codon in (bp_af_adj_dict.keys()):
        if codon_table[codon.upper()] == aa:
            if bp_af_adj_dict[codon] >= max_af:
                max_af = bp_af_adj_dict[codon]
                res_dict["bp_ref"] = codon

    # Calculating where in the codon the change occured
    if old_bp_ref[:1] != res_dict["bp_ref"][:1]:
        pos_of_change = 0
    elif old_bp_ref[1:2] != res_dict["bp_ref"][1:2]:
        pos_of_change = 1
    else:
        pos_of_change = 2

    # Adding the ref bp to the bp dicts: calculating the af to add
    af_sum_snps_at_position_of_change = 0
    af_adj_sum_snps_at_position_of_change = 0
    for bp in bp_af_dict:
        if bp[pos_of_change:pos_of_change + 1] != old_bp_ref[pos_of_change:pos_of_change + 1]:
            af_sum_snps_at_position_of_change += bp_af_dict[bp]
            af_adj_sum_snps_at_position_of_change += bp_af_adj_dict[bp]

    bp_af_dict[old_bp_ref] = format_af(1 - af_sum_snps_at_position_of_change)
    bp_af_adj_dict[old_bp_ref] = format_af(1 - af_adj_sum_snps_at_position_of_change)

    # Updating the Frequencies of ref
    res_dict["af"] = format_af(1 - aa_sum)
    res_dict["af_adj"] = format_af(1 - aa_adj_sum)

    # Deleting from the alterations dict
    del alterations_af_dict[aa]
    del alterations_af_adj_dict[aa]

    # Deleting from the bps dict
    del bp_af_dict[res_dict["bp_ref"]]
    del bp_af_adj_dict[res_dict["bp_ref"]]


def change_syn_ref_bp(new_bp_ref, old_bp_ref, res_dict, bp_af_dict, bp_af_adj_dict):
    # Updating bp_ref with the high freq. codon
    res_dict["bp_ref"] = new_bp_ref

    # Add the old ref to the bp dicts with its freq.
    bp_freq_adj_sum = sum(bp_af_adj_dict.values())
    bp_freq_sum = sum(bp_af_dict.values())
    bp_af_adj_dict[old_bp_ref] = (1 - bp_freq_adj_sum)
    bp_af_dict[old_bp_ref] = (1 - bp_freq_sum)

    # Delete the new ref from the bp dicts
    del bp_af_adj_dict[new_bp_ref]
    del bp_af_dict[new_bp_ref]


# A function that return a dict with the MAF info for the protein position and corresponding chromosomal location
def calc_exac_maf_data(chrom_pos_list, chrom_gene_table, protein_pos, aa, chrom, is_complementary, seq):

    res_dict = {}
    res_dict["chrom"] = chrom
    res_dict["chrom_pos"] = chrom_pos_list
    res_dict["prot_pos"] = protein_pos
    res_dict["bp_ref"] = seq
    # an counters lists
    an_dict = {k: [] for k in ["an", "an_adj", "an_afr", "an_amr", "an_eas", "an_fin", "an_nfe", "an_oth", "an_sas"]}
    res_dict.update(an_dict)
    # ac counters lists
    ac_dict = {k: [] for k in
               ["ac_adj", "ac_afr", "ac_amr", "ac_eas", "ac_fin", "ac_het", "ac_hom", "ac_nfe", "ac_oth", "ac_sas"]}
    res_dict.update(ac_dict)
    res_dict["SIFT"] = []
    res_dict["PolyPhen"] = []
    res_dict["clin_sig"] = []
    res_dict["SwissProt"] = []
    res_dict["ens_prot"] = []
    res_dict["mean_coverage"] = []

    # Initializing more variables
    indels_cnt = 0
    errors_cnt = 0
    alterations_af_dict = defaultdict(list)
    alterations_af_adj_dict = defaultdict(list)
    bp_af_dict = dict()
    bp_af_adj_dict = dict()
    bp_list = []  # Will not update when changing ref aa/bp, has to correspond to the population prob.

    aa_from_codon = str(Seq(res_dict["bp_ref"], generic_dna).translate())
    if aa == 'X':  # If Hmmer couldn't determine the Amino Acid,
        aa = aa_from_codon
    else:
        assert aa == aa_from_codon, "Amino Acid from Hmmer doesn't match that determined from codon sequence"

    res_dict["aa_ref"] = aa

    # Going over all 3 codon positions
    for i in range(len(chrom_pos_list)):
        chrom_pos = chrom_pos_list[i]
        alt_codon_pos = i

        # Retreiving relevant ExAC entry
        chrom_alter_table = chrom_gene_table[chrom_gene_table['POS'] == chrom_pos]

        if chrom_alter_table.empty:
            # No ExAC entry for this chromosome position - not adding alteration data
            continue
        else:
            # In case there are several alterations for that position, iterating
            for index, line in chrom_alter_table.iterrows():
                chrom_alter = line

                # Extracting ref and alt
                exac_ref_bp = chrom_alter['REF']
                exac_alt_bp = chrom_alter['ALT']

                # Extracting coverage
                res_dict["mean_coverage"].append(chrom_alter['COVERAGE'])

                # Check if indel - skip
                if 'ignore' in chrom_alter['COMMENTS']:
                    indels_cnt += 1
                    continue

                # Perform validation checks (comparing ExAC and HMMER data)
                exac_prot_data, exac_alt_aa, exac_alt_codon, errors = exac_validation_checks(chrom_alter, protein_pos,
                                                                                               aa, alt_codon_pos,
                                                                                               chrom_pos,
                                                                                               res_dict["bp_ref"])
                if errors:
                    errors_cnt += 1

                # Extracting ExAC allele frequency data
                af = chrom_alter["AF"]
                an = int(chrom_alter["AN"])
                an_adj = int(chrom_alter["AN_ADJ"])
                ac_adj = chrom_alter["AC_ADJ"]

                # Calculating the alteration relevant data
                alt_codon = create_alt_codon(exac_ref_bp, exac_alt_bp, res_dict["bp_ref"], alt_codon_pos, is_complementary)
                if len(alt_codon) != 3:
                    continue
                else:
                    alt_aa = codon_table[alt_codon.upper()]

                if exac_prot_data and exac_alt_aa != alt_aa:
                    print(f"{chrom_pos} Error: the ExAC alt aa {exac_alt_aa} doesn't match my alt aa calculation {alt_aa}")

                # Calculating the allele frequency adjusted
                af_adj = calculate_af_adj(an_adj, ac_adj)

                # Saving the bp with the frequency (for both syn and nonsyn)
                bp_af_dict[alt_codon] = format_af(af)
                bp_af_adj_dict[alt_codon] = format_af(af_adj)
                bp_list.append(alt_codon)

                res_dict["an"].append(an)
                res_dict["an_adj"].append(chrom_alter["AN_ADJ"])
                # res_dict["an_afr"].append(chrom_alter["AN_AFR"])
                # res_dict["an_amr"].append(chrom_alter["AN_AMR"])
                # res_dict["an_eas"].append(chrom_alter["AN_EAS"])
                # res_dict["an_fin"].append(chrom_alter["AN_FIN"])
                # res_dict["an_nfe"].append(chrom_alter["AN_NFE"])
                # res_dict["an_oth"].append(chrom_alter["AN_OTH"])
                # res_dict["an_sas"].append(chrom_alter["AN_SAS"])

                res_dict["ac_adj"].append(chrom_alter["AC_ADJ"])
                # res_dict["ac_afr"].append(chrom_alter["AC_AFR"])
                # res_dict["ac_amr"].append(chrom_alter["AC_AMR"])
                # res_dict["ac_eas"].append(chrom_alter["AC_EAS"])
                # res_dict["ac_fin"].append(chrom_alter["AC_FIN"])
                # res_dict["ac_het"].append(chrom_alter["AC_Het"])
                # res_dict["ac_hom"].append(chrom_alter["AC_Hom"])
                # res_dict["ac_nfe"].append(chrom_alter["AC_NFE"])
                # res_dict["ac_oth"].append(chrom_alter["AC_OTH"])
                # res_dict["ac_sas"].append(chrom_alter["AC_SAS"])

                # Non-synonymous(!!!) - logging the alteration in the dictionary
                if alt_aa != res_dict["aa_ref"]:
                    alterations_af_dict[alt_aa].append(format_af(float(af)))
                    alterations_af_adj_dict[alt_aa].append(format_af(af_adj))

                # Saving the SIFT and PolyPhen scores (for both syn and nonsyn)
                res_dict["SIFT"].append(chrom_alter["SIFT"])
                # res_dict["PolyPhen"].append(chrom_alter["PolyPhen"])
                # res_dict["clin_sig"].append(chrom_alter["clin_sig"])

                # Saving SwissProt id and Ensembl prot id
                res_dict["SwissProt"].append(chrom_alter["SWISSPROT"])
                res_dict["ens_prot"].append(chrom_alter["ENSP"])

    # Calculating the overall MAF from the alteration dicts
    res_dict["af"] = 0
    res_dict["af_adj"] = 0

    res_dict["aa_ref_orig"] = res_dict["aa_ref"]
    for aa in alterations_af_dict:
        aa_sum = sum(alterations_af_dict[aa])
        aa_adj_sum = sum(alterations_af_adj_dict[aa])

        # Checking if any alteration is above 0.5, and changing the ref accordingly
        if aa_sum > 0.5:
            # Update all relevant information to switching a ref aa
            change_ref_aa(res_dict, alterations_af_dict, alterations_af_adj_dict, aa, aa_sum, aa_adj_sum, bp_af_dict,
                          bp_af_adj_dict)
            break
        else:
            res_dict["af"] += aa_sum
            res_dict["af_adj"] += aa_adj_sum

        # Fix the AF format
        res_dict["af"] = format_af(res_dict["af"])
        res_dict["af_adj"] = format_af(res_dict["af_adj"])

    # Checking if any syn alteration bp is above 0.5 (nonsyn were already checked), and changing ref bp accordingly
    for codon in bp_af_adj_dict:
        if codon_table[codon] == res_dict["aa_ref"] and bp_af_adj_dict[codon] > 0.5:
            new_bp_ref = codon
            old_bp_ref = res_dict["bp_ref"]
            change_syn_ref_bp(new_bp_ref, old_bp_ref, res_dict, bp_af_dict, bp_af_adj_dict)

    res_dict["alterations_af_adj_dict"] = alterations_af_adj_dict
    res_dict["bp_af_dict"] = bp_af_dict
    res_dict["bp_af_adj_dict"] = bp_af_adj_dict
    res_dict["bp_list"] = bp_list

    return res_dict, indels_cnt, errors_cnt


if __name__ == '__main__':

    domain_name = os.path.splitext(os.path.basename(HMM_FILE))[0]

    with open(os.path.join(CANONIC_PROT_FOLDER, f'{domain_name}_canonic_prot.pik'), 'rb') as f:
        canonic_protein = pickle.load(f)

    domain_data = pd.read_csv(HMM_FILE, sep='\t', index_col=0, dtype={"chrom_num": str})

    chrom_gene_path = os.path.join(INDELS_FOLDER, domain_name)
    states_dict = defaultdict(list)

    # For each ensembl gene in the domain data - finding all the ExAC alterations
    ens_genes = domain_data[domain_data['chrom_num'].isin(CHROMOSOMES)].gene.unique().tolist()
    n_ens_genes = len(ens_genes)

    for i, ens_gene in enumerate(ens_genes):

        # Getting gene-chrom (ExAC) data tables
        chrom_gene_table = pd.read_csv(os.path.join(chrom_gene_path, ens_gene, f"chrom_gene_table.csv"), index_col=0,
                                       dtype={"AC": int, "AC_AFR": int, "AC_AMR": int, "AC_ADJ": int, "AC_EAS": int,
                                              "AC_FIN": int, "AC_Het": int, "AC_Hom": int, "AC_NFE": int, "AC_OTH": int,
                                              "AC_SAS": int, "AF": float, "AN": int, "AN_AFR": int,
                                              "AN_AMR": int, "AN_Adj": int, "AN_EAS": int, "AN_FIN": int, "AN_NFE": int,
                                              "AN_OTH": int, "AN_SAS": int, "PROTEIN_POSITION": str})
        chrom_gene_table.fillna('', inplace=True)
        exon_table = pd.read_csv(os.path.join(chrom_gene_path, ens_gene, "exon_table.csv"), index_col=0)

        # Filtering the domain data for this gene according to the canonical protein id
        canonic_prot = canonic_protein[ens_gene]
        canonic_prot_t = canonic_prot[:canonic_prot.find(".")]  # Trimming the ".#" at the end
        domain_gene_table = domain_data[domain_data["prot"] == canonic_prot]
        # Making sure that if two HMM-matches overlaps, the higher bit score will come first in the table
        domain_gene_table = domain_gene_table.sort_values(by="BitScore", ascending=False)
        domain_gene_name = domain_gene_table["hugoSymbol"].unique()[0]
        if len(domain_gene_table["hugoSymbol"].unique()) > 1:
            raise RuntimeError(" Error: " + ens_gene + ": more than one Hugo symbol")

        # Extracting neccessary information from the gene-domain data
        chrom_raw_data = domain_gene_table["chromosome"].unique()[0]  # there should be only one element here
        if len(domain_gene_table["chromosome"].unique()) > 1:
            raise RuntimeError(" Error: " + ens_gene + ": more than one chromosome raw data")
        targetid = domain_gene_table["#TargetID"].unique()[0]

        n_indels = n_errors = 0

        # Query HG19 and fill in a nucleobase sequence column in our exon_table dataframe
        chromosome = get_chromosome_number(chrom_raw_data)
        is_complementary = 'complement' in chrom_raw_data

        exon_table['seq'] = retrieve_exon_seq(exon_table.start_pos, exon_table.end_pos, chromosome, HG19_FILE,
                                              complement=is_complementary)

        # Iterating over the amino-acids of the protein
        prot_len = int(domain_gene_table["length"].unique()[0])

        for index, row in domain_gene_table.sort_values(by='TargetStart').iterrows():

            target_start = row.TargetStart
            target_end = row.TargetEnd
            hmm_pos = (row["HMM_Pos"]).split(",")
            target_seq = list(row["Target_Seq"])

            if '-' in target_seq:
                # Remove deletions from both lists
                indices = [i for i, x in enumerate(target_seq) if x == "-"]
                target_seq = [i for j, i in enumerate(target_seq) if j not in indices]
                hmm_pos = [i for j, i in enumerate(hmm_pos) if j not in indices]

            for protein_pos in range(target_start, target_end + 1):
                index_inside_match = int(protein_pos - target_start)
                hmm_state_text = hmm_pos[index_inside_match]

                # If there's a match to HMM-state: find the corresponding codon bps chromosome positions
                if is_number(hmm_state_text):
                    hmm_state = int(hmm_state_text)
                    aa = (target_seq[index_inside_match]).upper()
                    chrom_pos_list, seq = find_chrom_bps(protein_pos, exon_table, chrom_raw_data)

                    # Analysis of the amino-acid MAF and related data, returned in a dictionary
                    info_dict, indels_cnt, errors_cnt = calc_exac_maf_data(chrom_pos_list, chrom_gene_table,
                                                                           protein_pos,
                                                                           aa, chromosome, is_complementary, seq)
                    info_dict["ens_gene"] = ens_gene
                    info_dict["prot_len"] = prot_len

                    # Adding the dictionary to the HMM-state list
                    states_dict[hmm_state].append(info_dict)
                    n_indels += indels_cnt
                    n_errors += errors_cnt

        print(f"Finished gene {i}/{n_ens_genes}. Indels={n_indels}, Errors={n_errors}")

    with open(OUTPUT_FILE, 'wb') as f:
        pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
