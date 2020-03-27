import pandas as pd
import numpy as np
from indels_func import is_indel, table_editing, indel_type
from mapping_func import create_exon_pos_table
import pickle
from dsprint.core import CHROMOSOMES, INSTANCE_THRESHOLD


DOMAIN_STATES_DF = "/media/vineetb/t5-vineetb/dsprint/out/pfam/32/domains_stats_df.csv"


if __name__ == '__main__':
    in_path = curr_dir[0] + "/../3.parse_HMMER/hmm_domains/pfam-v" + pfam_version + "/"
    chrom_path = curr_dir[0] + "/../1.parse_ExAC/parsed_filtered/"
    chrom_filename = "parsed_filtered_chrom"
    out_path = curr_dir[0] + "/domain_gene_exac/pfam-v" + pfam_version + "/"

    domain_stats_df = pd.read_csv(DOMAIN_STATES_DF, sep='\t', index_col=0)
    all_domains_list = domain_stats_df.index.tolist()
    filtered_domain_list = domain_stats_df[domain_stats_df.instances > INSTANCE_THRESHOLD].index.tolist()

    for chrom in CHROMOSOMES:

        # Loading the ExAC parsed data of this chromosome
        fields = ["chrom", "pos", "ref", "alt", "filter", "AC", "AC_AFR", "AC_AMR", "AC_Adj", "AC_EAS", "AC_FIN", "AC_Het",
                  "AC_Hom", "AC_NFE", "AC_OTH", "AC_SAS",
                  "AF", "AN", "AN_AFR", "AN_AMR", "AN_Adj", "AN_EAS", "AN_FIN", "AN_NFE", "AN_OTH", "AN_SAS", "prot_pos",
                  "amino_acids", "codons", "ENSP",
                  "SWISSPROT", "SIFT", "PolyPhen", "clin_sig", "conseq", "feature", "strand", "coverage"]
        # More fields that might be interesting:  feature_type, domains, , gene

        chrom_csv = pd.read_csv(chrom_path + chrom_filename + chrom + ".csv", sep='\t', usecols=fields,
                                dtype={"AC": int, "AC_AFR": int, "AC_AMR": int, "AC_Adj": int, "AC_EAS": int,
                                       "AC_FIN": int, "AC_Het": int, "AC_Hom": int, "AC_NFE": int, "AC_OTH": int,
                                       "AC_SAS": int, "AF": float, "AN": int, "AN_AFR": int,
                                       "AN_AMR": int, "AN_Adj": int, "AN_EAS": int, "AN_FIN": int, "AN_NFE": int,
                                       "AN_OTH": int, "AN_SAS": int, "prot_pos": str})
        chrom_csv = chrom_csv.sort_values(by=["pos"])
        chrom_csv = chrom_csv.reset_index(drop=True)
        chrom_csv.fillna('', inplace=True)
        chrom_csv["comments"] = ""

        for domain_name in all_domains_list:

            # Get the canonic protein id for the domain
            with open(curr_dir[
                          0] + "/../4.parse_Uniprot/domains_canonic_prot/pfam-v" + pfam_version + "/" + domain_name + "_canonic_prot.pik",
                      'rb') as handle:
                canonic_protein = pickle.load(handle)

            filename = domain_name + ".csv"
            domain_data = pd.read_csv(in_path + filename, sep='\t', index_col=0, dtype={"chrom_num": str})
            # Sort the domain data
            sorted_domain_data = domain_data.sort_values(by=["chrom_num", "gene", "TargetStart"])
            sorted_domain_data = sorted_domain_data.reset_index(drop=True)

            # Filtering the domain data relevant to this chromosome
            domain_chrom_data = sorted_domain_data[sorted_domain_data["chrom_num"] == chrom]

            # Getting a list of all the relevant ensembl gene ids for this chromosome
            domain_ens_genes = (domain_chrom_data["gene"]).unique()
            # domain_ens_genes_all.extend(domain_ens_genes)

            # If this domain have no genes on this chromosome, continue to the next domain
            if (len(domain_ens_genes) == 0):
                print
                "Finished domain " + domain_name
                continue

            for ens_gene in domain_ens_genes:

                # Filtering the domain data for this gene according to the canonical protein id
                canonic_prot = canonic_protein[ens_gene]
                canonic_prot_t = canonic_prot[:canonic_prot.find(".")]  # Trimming the ".#" at the end
                domain_gene_table = domain_chrom_data[domain_chrom_data["prot"] == canonic_prot]
                # Making sure that if two HMM-matches overlaps, the higher bit score will come first in the table
                domain_gene_table = domain_gene_table.sort_values(by="BitScore", ascending=False)
                domain_gene_name = domain_gene_table["hugoSymbol"].unique()[0]
                if (len(domain_gene_table["hugoSymbol"].unique()) > 1):
                    print
                    functionNameAsString + " Error: " + ens_gene + ": more than one Hugo symbol"  # sanity check

                # Creating a table of the exons for this gene, according to the canonical protein
                chrom_raw_data = domain_gene_table["chromosome"].unique()[0]  # there should be only one element here
                if (len(domain_gene_table["chromosome"].unique()) > 1):
                    print
                    functionNameAsString + " Error: " + ens_gene + ": more than one chromosome raw data"  # sanity check
                targetid = domain_gene_table["#TargetID"].unique()[0]
                exon_table = create_exon_pos_table(chrom_raw_data, targetid)

                # Filtering the chromosome data to the gene exons region
                exons_start_pos = min(exon_table["start_pos"][0], exon_table["start_pos"][
                    len(exon_table) - 1])  # in case of complelemt, the minimal position could be at the last row
                exons_end_pos = max(exon_table["end_pos"][0], exon_table["end_pos"][
                    len(exon_table) - 1])  # in case of complelemt, the maximal position could be at the first row
                chrom_gene_table = \
                chrom_csv[chrom_csv["pos"] >= int(exons_start_pos)][chrom_csv["pos"] <= int(exons_end_pos)][
                    chrom_csv["ENSP"] == canonic_prot_t]
                chrom_gene_table = chrom_gene_table.reset_index(drop=True)

                # Adding chrom column to the table
                chrom_gene_table["chrom"] = chrom

                # Handling indels
                indels_table = table_editing(chrom_gene_table)

                !mkdir - p
                domain_gene_exac / pfam - v31 /$domain_name
                !mkdir - p
                domain_gene_exac / pfam - v31 /$domain_name /$ens_gene
                chrom_gene_table.to_csv(out_path + domain_name + "/" + ens_gene + "/chrom_gene_table.csv", sep='\t')
                indels_table.to_csv(out_path + domain_name + "/" + ens_gene + "/indels_table.csv", sep='\t')
                exon_table.to_csv(out_path + domain_name + "/" + ens_gene + "/exon_table.csv", sep='\t')

                print
                "Finished gene " + ens_gene
            print
            "Finished domain " + domain_name
        print
        "Finished chromosome " + chrom