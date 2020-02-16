import pandas as pd
from collections import defaultdict


try:
    INPUT_FILE = snakemake.input[0]
    OUTPUT_FILES = snakemake.output
    CHROMOSOMES = [O[O.rindex('_')+1:O.rindex('.')] for O in OUTPUT_FILES]
except NameError:
    import sys
    if len(sys.argv) < 3:
        print('Usage: <script> <input_folder> <chromosome_number>')
        sys.exit(0)

    INPUT_FILE, CHROMOSOMES = sys.argv[1], [sys.argv[2]]
    OUTPUT_FILES = [f"/media/vineetb/t5-vineetb/dsprint/out/parsed_chrom_{C}.csv" for C in CHROMOSOMES]


def update_main_fields(line_parts, d):
    d["chrom"].append(line_parts[CHROM])
    d["pos"].append(int(line_parts[POS]))
    d["id"].append(line_parts[ID])
    d["ref"].append(line_parts[REF])
    d["qual"].append(line_parts[QUAL])
    d["filter"].append(line_parts[FILTER])
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
    d["AN"].append(info[AN_beg + 3:AN_end])

    # AN_adj = Adjusted Allele Number
    AN_adj_beg = info.find("AN_Adj=")
    AN_adj_end = info.find(";", AN_adj_beg)
    d["AN_Adj"].append(info[AN_adj_beg + 7:AN_adj_end])

    # DP = "Approximate read depth
    DP_beg = info.find("DP=")
    DP_end = info.find(";", DP_beg)
    d["DP"].append(info[DP_beg + 3:DP_end])

    return AC_list, AC_adj_list, AF_list


def fill_empty_fields(line_parts, alt_list, d):
    """A function that updates empty strings for CSQ fields when there's no CSQ data."""

    # Adding one line per each variation (can be several even when tere's no CSQ)
    for i in range(len(alt_list)):
        # Update the main fields
        (AC_list, AC_adj_list, AF_list) = update_main_fields(line_parts, d)
        d["AC"].append(AC_list[i])
        d["AC_Adj"].append(AC_adj_list[i])
        d["AF"].append(AF_list[i])

        # Update the alt
        d["alt"].append(alt_list[i])

        # Update the rest of fields with empty string
        d["gene"].append("")
        d["feature"].append("")
        d["feature_type"].append("")
        d["conseq"].append("")
        d["prot_pos"].append("")
        d["amino_acids"].append("")
        d["codons"].append("")
        d["strand"].append("")
        d["ENSP"].append("")
        d["SWISSPROT"].append("")
        d["SIFT"].append("")
        d["PolyPhen"].append("")
        d["exon"].append("")
        d["intron"].append("")
        d["domains"].append("")
        d["clin_sig"].append("")


if __name__ == '__main__':

    # Positions in the VCF record
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = tuple(range(8))

    # CSQ positions
    # http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html?redirect=no
    GENE = 1  # Ensembl stable ID of affected gene
    FEATURE = 2  # Ensembl stable ID of feature
    FEATURE_TYPE = 3  # type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature
    CONSEQ = 4  # consequence type of this variant
    PROT_POS = 7  # relative position of amino acid in protein
    AMINO_ACIDS = 8  # the change. only given if the variant affects the protein-coding sequence
    CODONS = 9  # the alternative codons with the variant base in upper case
    ALLELE_NUM = 11  # Allele number from input; 0 is reference, 1 is first alternate etc
    STRAND = 13  # the DNA strand (1 or -1) on which the transcript/feature lies
    ENSP = 20  # the Ensembl protein identifier of the affected transcript
    SWISSPROT = 21  # UniProtKB/Swiss-Prot identifier of protein product
    SIFT = 24  # the SIFT prediction and/or score, with both given as prediction(score)
    POLYPHEN = 25  # the PolyPhen prediction and/or score
    EXON = 26  # the exon number (out of total number)
    INTRON = 27  # the intron number (out of total number)
    DOMAINS = 28  # the source and identifer of any overlapping protein domains
    GMAF = 31  # Non-reference allele and frequency of existing variant in 1000 Genomes
    CLIN_SIG = 38  # Clinical significance of variant from dbSNP: http://varianttools.sourceforge.net/Annotation/DbSNP
    # Variant Clinical Significance, 0 - unknown, 1 -
    # untested, 2 - non-pathogenic, 3 - probable-non-pathogenic, 4 - probable-pathogenic,
    # 5 - pathogenic, 6 - drug-response, 7 - histocompatibility, 255 - other

    # https://macarthurlab.org/2016/03/17/reproduce-all-the-figures-a-users-guide-to-exac-part-2/#multi-allelic-enriched-regions
    multi_allelic_regions = {'14': [[106329000, 106331000], [107178000, 107180000]],
                             '2': [[89160000, 89162000]],
                             '17': [[18967000, 18968000], [19091000, 19092000]],
                             '22': [[23223000, 23224000]],
                             '1': [[152975000, 152976000]]}
    # A flag for removing the multi-allelic regions
    removeMultiAllelic = True

    for chromosome, output_file in zip(CHROMOSOMES, OUTPUT_FILES):

        # Read the file
        vcf_file = open(INPUT_FILE, 'r')

        # Process meta-data
        metadata_dict = {}
        data_flag = False
        for line in vcf_file:
            if line.startswith("##"):
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
                    if f_beg < 0:
                        f_beg = line.find(":", f_key)  # When parsing reference line
                    f_end = line.find(",", f_beg)
                    if f_end < 0:
                        f_end = line.find(">")
                    if f_end < 0:  # When parsing reference line
                        f_end = line.find("\n")
                    val[f] = line[f_beg + 1:f_end]

                # Adding to the metadata dictionary
                if key not in metadata_dict:
                    metadata_dict[key] = [val]
                else:
                    metadata_dict[key].append(val)

            # Processing the data starting the next line
            elif line.startswith("#CHROM"):
                data_flag = True
                break

        # Arrange the INFO metadata to a data-frame
        info_df = pd.DataFrame(metadata_dict["INFO"])
        info_df = info_df.sort_values("ID")

        d = defaultdict(list)
        for line_no, line in enumerate(vcf_file):

            line_parts = line.split("\t")

            if line_parts[CHROM] != chromosome:
                continue

            # Excluding multi-allelic regions
            if removeMultiAllelic:
                chrom = line_parts[CHROM]
                if chrom in multi_allelic_regions.keys():
                    pos = int(line_parts[POS])
                    regions = multi_allelic_regions[chrom]
                    for region in regions:
                        if region[0] <= pos <= region[1]:
                            # Multi-allelic region - excluding from analysis
                            continue

            # Extracting alt
            alt_list = line_parts[ALT].split(",")

            # Extracting fields from the info
            info = line_parts[INFO]

            # CSQ = Consequence type as predicted by VEP
            CSQ_beg = info.find("CSQ=")
            if CSQ_beg == -1:
                # NO CSQ data: just fill in empty strings instead
                fill_empty_fields(line_parts, alt_list, d)
            else:
                CSQ_data = info[CSQ_beg + 4:]
                CSQ_features = CSQ_data.split(",")

                for CSQ in CSQ_features:
                    CSQ_PARTS = CSQ.split('|')
                    # Update the main fields for each CSQ feature (so each CSQ will appear in a different line)
                    AC_list, AC_adj_list, AF_list = update_main_fields(line_parts, d)

                    # Allele_num for deciding which alt, AC and AF to add
                    allele_num = int(CSQ_PARTS[ALLELE_NUM])
                    # Adding the corresponding alt
                    if allele_num == 0:
                        print("allele num = 0 " + line_parts[POS])  # Making sure all the features correspond to an alt (0 = ref)
                    else:
                        d["alt"].append(alt_list[allele_num - 1])
                        d["AC"].append(AC_list[allele_num - 1])
                        d["AC_Adj"].append(AC_adj_list[allele_num - 1])
                        d["AF"].append(AF_list[allele_num - 1])

                    d["gene"].append(CSQ_PARTS[GENE])
                    d["feature"].append(CSQ_PARTS[FEATURE])
                    d["feature_type"].append(CSQ_PARTS[FEATURE_TYPE])
                    d["conseq"].append(CSQ_PARTS[CONSEQ])
                    d["prot_pos"].append(CSQ_PARTS[PROT_POS])
                    d["amino_acids"].append(CSQ_PARTS[AMINO_ACIDS])
                    d["codons"].append(CSQ_PARTS[CODONS])
                    d["strand"].append(CSQ_PARTS[STRAND])
                    d["ENSP"].append(CSQ_PARTS[ENSP])
                    d["SWISSPROT"].append(CSQ_PARTS[SWISSPROT])
                    d["SIFT"].append(CSQ_PARTS[SIFT])
                    d["PolyPhen"].append(CSQ_PARTS[POLYPHEN])
                    d["exon"].append(CSQ_PARTS[EXON])
                    d["intron"].append(CSQ_PARTS[INTRON])
                    d["domains"].append(CSQ_PARTS[DOMAINS])
                    d["clin_sig"].append(CSQ_PARTS[CLIN_SIG])

        # Only needed for checksum
        columns = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "AC", "AC_Adj", "AF", "AN", "AN_Adj", "DP",
                     "gene",
                     "feature", "feature_type", "conseq", "prot_pos",
                     "amino_acids", "codons", "strand", "ENSP", "SWISSPROT", "SIFT", "PolyPhen", "exon", "intron",
                     "domains",
                     "clin_sig"]

        pd.DataFrame(d).to_csv(output_file, sep='\t', columns=columns)

        vcf_file.close()

