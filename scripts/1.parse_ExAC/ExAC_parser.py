import gzip
import pandas as pd
import os.path
from collections import defaultdict


EXAC_PATH = '/media/vineetb/t5-vineetb/dsprint/ExAC.r1.sites.vep.vcf.gz'
OUTPUT_PATH = '/media/vineetb/t5-vineetb/dsprint/out'

# Positions in the VCF record
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = tuple(range(8))

# fields for extraction (it's important that those names will match *exactly* the dictionary keys in "data_dict")
headlines = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "AC", "AC_Adj", "AF", "AN", "AN_Adj", "DP", "gene",
             "feature", "feature_type", "conseq", "prot_pos",
             "amino_acids", "codons", "strand", "ENSP", "SWISSPROT", "SIFT", "PolyPhen", "exon", "intron", "domains",
             "clin_sig"]

# CSQ positions
# http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html?redirect=no
GENE = 0  # Ensembl stable ID of affected gene
FEATURE = 1  # Ensembl stable ID of feature
FEATURE_TYPE = 2  # type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature
CONSEQ = 3  # consequence type of this variant
PROT_POS = 6  # relative position of amino acid in protein
AMINO_ACIDS = 7  # the change. only given if the variant affects the protein-coding sequence
CODONS = 8  # the alternative codons with the variant base in upper case
ALLELE_NUM = 10  # Allele number from input; 0 is reference, 1 is first alternate etc
STRAND = 12  # the DNA strand (1 or -1) on which the transcript/feature lies
ENSP = 19  # the Ensembl protein identifier of the affected transcript
SWISSPROT = 20  # UniProtKB/Swiss-Prot identifier of protein product
SIFT = 23  # the SIFT prediction and/or score, with both given as prediction(score)
POLYPHEN = 24  # the PolyPhen prediction and/or score
EXON = 25  # the exon number (out of total number)
INTRON = 26  # the intron number (out of total number)
DOMAINS = 27  # the source and identifer of any overlapping protein domains
GMAF = 30  # Non-reference allele and frequency of existing variant in 1000 Genomes
CLIN_SIG = 37  # Clinical significance of variant from dbSNP: http://varianttools.sourceforge.net/Annotation/DbSNP
# Variant Clinical Significance, 0 - unknown, 1 -
# untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic,
# 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other


# Defined here: https://macarthurlab.org/2016/03/17/reproduce-all-the-figures-a-users-guide-to-exac-part-2/#multi-allelic-enriched-regions
multi_allelic_regions = {'14': [[106329000, 106331000], [107178000, 107180000]],
                         '2': [[89160000, 89162000]],
                         '17': [[18967000, 18968000], [19091000, 19092000]],
                         '22': [[23223000, 23224000]],
                         '1': [[152975000, 152976000]]}
# A flag for removing the multi-allelic regions
removeMultiAllelic = True

# Read the file
vcf_file = gzip.open(EXAC_PATH, 'r')

# Process meta-data
metadata_dict = {}
data_flag = False
for line in vcf_file:
    if line[0:2] == "##":
        # assign keys according to the format
        key = line[2:line.index('=')]
        if key == "ALT":
            val = dict.fromkeys(["ID", "Description"])
        elif key == "FILTER":
            val = dict.fromkeys(["ID", "Description"])
        elif key == "FORMAT":
            val = dict.fromkeys(["ID", "Number", "Type", "Description"])
        elif key == "INFO":
            val = dict.fromkeys(["ID", "Number", "Type", "Description"])
        elif key == "contig":
            val = dict.fromkeys(["ID", "length"])
        elif key == "reference":
            val = dict.fromkeys(["file"])
        # Not processing other metadata types
        else:
            continue

        # fill in the data
        for f in val.keys():
            f_key = line.find(f)
            f_beg = line.find("=", f_key)
            if (f_beg < 0):
                f_beg = line.find(":", f_key)  # When parsing reference line
            f_end = line.find(",", f_beg)
            if (f_end < 0):
                f_end = line.find(">")
            if (f_end < 0):  # When parsing reference line
                f_end = line.find("\n")
            val[f] = line[f_beg + 1:f_end]

        # Adding to the metadata dictionary
        if not metadata_dict.has_key(key):
            metadata_dict[key] = [val]
        else:
            metadata_dict[key].append(val)

    # Processing the data starting the next line
    elif line[0:6] == "#CHROM":
        data_flag = True
        break

# Arrange the INFO metadata to a data-frame
info_df = pd.DataFrame(metadata_dict["INFO"])
info_df = info_df.sort_values("ID")


def find_nth(s, x, n=0, overlap=False):
    """A function to find the nth position on a substring x in the string s"""

    l = 1 if overlap else len(x)
    i = -l
    for c in xrange(n + 1):
        i = s.find(x, i + l)
        if i < 0:
            break
    return i


def data_to_df_csv(data_dict, headlines, chrom_num, output_folder):
    """A function that saves the data dictionary to DataFrame and then to .csv"""

    # Creating a data_frame from all the parsed values of the chromosome
    df = pd.DataFrame([data_dict[h] for h in headlines])
    df = df.transpose()
    df.columns = headlines

    # Saving the df to a file
    df.to_csv(os.path.join(output_folder, "parsed_chrom" + chrom_num + ".csv"), sep='\t')


def update_main_fields(line_parts, data_dict):
    """A function that extract the main fields from line_parts and adds them to the data dictionary.
    Returining the extracted info field for further processing"""

    # Extracting Chromosome number
    data_dict["chrom"].append(line_parts[CHROM])

    # Extracting position
    data_dict["pos"].append(int(line_parts[POS]))

    # Extracting id
    data_dict["id"].append(line_parts[ID])

    # Extracting ref
    data_dict["ref"].append(line_parts[REF])

    # Extracting quality
    data_dict["qual"].append(line_parts[QUAL])

    # Extracting filter
    data_dict["filter"].append(line_parts[FILTER])

    # Extracting fields from the info
    info = line_parts[INFO]

    # AC = Allele Count
    AC_beg = info.find("AC=")
    AC_end = info.find(";", AC_beg)
    AC_list = (info[AC_beg + 3:AC_end]).split(",")

    # AC_adjusted = Adjusted Allele Count
    AC_adj_beg = info.find("AC_Adj=")
    AC_adj_end = info.find(";", AC_adj_beg)
    AC_adj_list = (info[AC_adj_beg + 7:AC_adj_end]).split(",")

    # AF = Allele Frequency
    AF_beg = info.find("AF=")
    AF_end = info.find(";", AF_beg)
    AF_list = (info[AF_beg + 3:AF_end]).split(",")

    # AN = Allele Number
    AN_beg = info.find("AN=")
    AN_end = info.find(";", AN_beg)
    data_dict["AN"].append(info[AN_beg + 3:AN_end])

    # AN_adj = Adjusted Allele Number
    AN_adj_beg = info.find("AN_Adj=")
    AN_adj_end = info.find(";", AN_adj_beg)
    data_dict["AN_Adj"].append(info[AN_adj_beg + 7:AN_adj_end])

    # DP = "Approximate read depth
    DP_beg = info.find("DP=")
    DP_end = info.find(";", DP_beg)
    data_dict["DP"].append(info[DP_beg + 3:DP_end])

    return (AC_list, AC_adj_list, AF_list)


def fill_empty_fields(line_parts, alt_list, data_dict):
    """A function that updates empty strings for CSQ fields when there's no CSQ data."""

    # Adding one line per each variation (can be several even when tere's no CSQ)
    for i in range(len(alt_list)):
        # Update the main fields
        (AC_list, AC_adj_list, AF_list) = update_main_fields(line_parts, data_dict)
        data_dict["AC"].append(AC_list[i])
        data_dict["AC_Adj"].append(AC_adj_list[i])
        data_dict["AF"].append(AF_list[i])

        # Update the alt
        data_dict["alt"].append(alt_list[i])

        # Update the rest of fields with empty string
        data_dict["gene"].append("")
        data_dict["feature"].append("")
        data_dict["feature_type"].append("")
        data_dict["conseq"].append("")
        data_dict["prot_pos"].append("")
        data_dict["amino_acids"].append("")
        data_dict["codons"].append("")
        data_dict["strand"].append("")
        data_dict["ENSP"].append("")
        data_dict["SWISSPROT"].append("")
        data_dict["SIFT"].append("")
        data_dict["PolyPhen"].append("")
        data_dict["exon"].append("")
        data_dict["intron"].append("")
        data_dict["domains"].append("")
        data_dict["clin_sig"].append("")


# Process the data records of the vcf and save each chromosome to a seperate file
chromosome_iter = '1'
data_dict = defaultdict(list)
multi_allele_skipped = 0

for line_no, line in enumerate(vcf_file):
    line_parts = line.split("\t")

    # Excluding multi-allelic regions
    if (removeMultiAllelic == True):
        chrom = line_parts[CHROM]
        if (chrom in multi_allelic_regions.keys()):
            pos = int(line_parts[POS])
            regions = multi_allelic_regions[chrom]
            for region in regions:
                if (pos >= region[0] and pos <= region[1]):
                    # Multi-allelic region - excluding from analysis
                    multi_allele_skipped += 1
                    continue

    # If the next line belongs to a different chromosome - saving to file
    if line_parts[CHROM] != chromosome_iter:
        data_to_df_csv(data_dict, headlines, chromosome_iter, OUTPUT_PATH)
        # Initializing the data dictionary
        data_dict = defaultdict(list)
        print "finished chromosome" + chromosome_iter
        chromosome_iter = line_parts[CHROM]

    # Extracting alt
    alt_list = line_parts[ALT].split(",")

    # Extracting fields from the info
    info = line_parts[INFO]

    # CSQ = Consequence type as predicted by VEP
    CSQ_beg = info.find("CSQ=")
    if (CSQ_beg == -1):
        # NO CSQ data: just fill in empty strings instead
        fill_empty_fields(line_parts, alt_list, data_dict)
    else:
        CSQ_data = info[CSQ_beg + 4:]
        CSQ_features = CSQ_data.split(",")

        for CSQ in CSQ_features:
            # Update the main fields for each CSQ feature (so each CSQ will appear in a different line)
            (AC_list, AC_adj_list, AF_list) = update_main_fields(line_parts, data_dict)

            # Allele_num for deciding which alt, AC and AF to add
            allele_num_beg = find_nth(CSQ, "|", ALLELE_NUM)
            allele_num_end = CSQ.find("|", allele_num_beg + 1)
            allele_num = int(CSQ[allele_num_beg + 1:allele_num_end])
            # Adding the corresponding alt
            if (allele_num == 0):
                print "allele num = 0 " + line_parts[POS]  # Making sure all the features correspond to an alt (0 = ref)
            else:
                data_dict["alt"].append(alt_list[allele_num - 1])
                data_dict["AC"].append(AC_list[allele_num - 1])
                data_dict["AC_Adj"].append(AC_adj_list[allele_num - 1])
                data_dict["AF"].append(AF_list[allele_num - 1])

            # Gene
            gene_beg = find_nth(CSQ, "|", GENE)
            gene_end = CSQ.find("|", gene_beg + 1)
            data_dict["gene"].append(CSQ[gene_beg + 1:gene_end])
            # Feature
            feature_beg = find_nth(CSQ, "|", FEATURE)
            feature_end = CSQ.find("|", feature_beg + 1)
            data_dict["feature"].append(CSQ[feature_beg + 1:feature_end])
            # Feature Type
            feature_type_beg = find_nth(CSQ, "|", FEATURE_TYPE)
            feature_type_end = CSQ.find("|", feature_type_beg + 1)
            data_dict["feature_type"].append(CSQ[feature_type_beg + 1:feature_type_end])
            # Consequence
            conseq_beg = find_nth(CSQ, "|", CONSEQ)
            conseq_end = CSQ.find("|", conseq_beg + 1)
            data_dict["conseq"].append(CSQ[conseq_beg + 1:conseq_end])
            # Protein_pos
            prot_pos_beg = find_nth(CSQ, "|", PROT_POS)
            prot_pos_end = CSQ.find("|", prot_pos_beg + 1)
            data_dict["prot_pos"].append(CSQ[prot_pos_beg + 1:prot_pos_end])
            # Amino Acids
            aa_beg = find_nth(CSQ, "|", AMINO_ACIDS)
            aa_end = CSQ.find("|", aa_beg + 1)
            data_dict["amino_acids"].append(CSQ[aa_beg + 1:aa_end])
            # Codons
            codons_beg = find_nth(CSQ, "|", CODONS)
            codons_end = CSQ.find("|", codons_beg + 1)
            data_dict["codons"].append(CSQ[codons_beg + 1:codons_end])
            # Strand
            strand_beg = find_nth(CSQ, "|", STRAND)
            strand_end = CSQ.find("|", strand_beg + 1)
            data_dict["strand"].append(CSQ[strand_beg + 1:strand_end])
            # ENSP
            ENSP_beg = find_nth(CSQ, "|", ENSP)
            ENSP_end = CSQ.find("|", ENSP_beg + 1)
            data_dict["ENSP"].append(CSQ[ENSP_beg + 1:ENSP_end])
            # Swissprot
            swiss_beg = find_nth(CSQ, "|", SWISSPROT)
            swiss_end = CSQ.find("|", swiss_beg + 1)
            data_dict["SWISSPROT"].append(CSQ[swiss_beg + 1:swiss_end])
            # SIFT
            sift_beg = find_nth(CSQ, "|", SIFT)
            sift_end = CSQ.find("|", sift_beg + 1)
            data_dict["SIFT"].append(CSQ[sift_beg + 1:sift_end])
            # PolyPhen
            polyphen_beg = find_nth(CSQ, "|", POLYPHEN)
            polyphen_end = CSQ.find("|", polyphen_beg + 1)
            data_dict["PolyPhen"].append(CSQ[polyphen_beg + 1:polyphen_end])
            # Exon
            exon_beg = find_nth(CSQ, "|", EXON)
            exon_end = CSQ.find("|", exon_beg + 1)
            data_dict["exon"].append(CSQ[exon_beg + 1:exon_end])
            # Intron
            intron_beg = find_nth(CSQ, "|", INTRON)
            intron_end = CSQ.find("|", intron_beg + 1)
            data_dict["intron"].append(CSQ[intron_beg + 1:intron_end])
            # Domains
            domains_beg = find_nth(CSQ, "|", DOMAINS)
            domains_end = CSQ.find("|", domains_beg + 1)
            data_dict["domains"].append(CSQ[domains_beg + 1:domains_end])
            # clin_sig
            clin_sig_beg = find_nth(CSQ, "|", CLIN_SIG)
            clin_sig_end = CSQ.find("|", clin_sig_beg + 1)
            data_dict["clin_sig"].append(CSQ[clin_sig_beg + 1:clin_sig_end])

# Saving the data of the last chromosome
data_to_df_csv(data_dict, headlines, chromosome_iter, OUTPUT_PATH)
print "finished chromosome" + chromosome_iter

vcf_file.close()

print multi_allele_skipped

