import os.path
import re
import gzip
import pandas as pd
from collections import defaultdict

import dsprint


# Positions in the VCF record - these are fixed as per .vcf format and can/should be hardcoded
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = tuple(range(8))

# The 17 CSQ Fields we are interested in
CSQs = ('GENE', 'FEATURE', 'FEATURE_TYPE', 'CONSEQUENCE', 'PROTEIN_POSITION', 'AMINO_ACIDS', 'CODONS', 'ALLELE_NUM',
        'STRAND', 'ENSP', 'SWISSPROT', 'SIFT', 'POLYPHEN', 'EXON', 'INTRON', 'DOMAINS', 'CLIN_SIG')

# https://macarthurlab.org/2016/03/17/reproduce-all-the-figures-a-users-guide-to-exac-part-2/#multi-allelic-enriched-regions
MULTI_ALLELIC_REGIONS = {
    # Keep the keys strings to account for x/y chromosomes
    '1':  [(152975000, 152976000)],
    '2':  [(89160000, 89162000)],
    '14': [(106329000, 106331000), (107178000, 107180000)],
    '17': [(18967000, 18968000), (19091000, 19092000)],
    '22': [(23223000, 23224000)],
}

REMOVE_MULTI_ALLELIC = True

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 4:
        print('Usage: <script> <input_file> <chromosome_number> <output_file>')
        sys.exit(0)

    INPUT_FILE, CHROMOSOMES, OUTPUT_FILES = sys.argv[1:]
    CHROMOSOMES, OUTPUT_FILES = [CHROMOSOMES], [OUTPUT_FILES]
else:
    INPUT_FILE = snakemake.input[0]
    OUTPUT_FILES = snakemake.output
    CHROMOSOMES = [os.path.splitext(os.path.basename(o)[len('parsed_chrom')])[0] for o in OUTPUT_FILES]


def update_main_fields(line_parts, d):

    for k in ('CHROM', 'POS', 'ID', 'REF', 'QUAL', 'FILTER'):
        d[k].append(line_parts[globals()[k]])

    info = line_parts[INFO]

    # AC = allele count in genotypes, for each ALT allele, in the same order as listed
    AC_beg = info.find('AC=')
    AC_end = info.find(';', AC_beg)
    AC_list = (info[AC_beg + 3:AC_end]).split(',')

    # AC_adjusted = Adjusted Allele Count
    AC_adj_beg = info.find('AC_Adj=')
    AC_adj_end = info.find(';', AC_adj_beg)
    AC_adj_list = (info[AC_adj_beg + 7:AC_adj_end]).split(',')

    # AF = allele frequency for each ALT allele in the same order as listed (use this when estimated from primary data,
    # not called genotypes)
    AF_beg = info.find('AF=')
    AF_end = info.find(';', AF_beg)
    AF_list = (info[AF_beg + 3:AF_end]).split(',')

    # AN = total number of alleles in called genotypes
    AN_beg = info.find('AN=')
    AN_end = info.find(';', AN_beg)
    d['AN'].append(info[AN_beg + 3:AN_end])

    # AN_adj = Adjusted Allele Number
    AN_adj_beg = info.find('AN_Adj=')
    AN_adj_end = info.find(';', AN_adj_beg)
    d['AN_ADJ'].append(info[AN_adj_beg + 7:AN_adj_end])

    # DP = combined depth across samples, e.g. DP=154
    DP_beg = info.find('DP=')
    DP_end = info.find(';', DP_beg)
    d['DP'].append(info[DP_beg + 3:DP_end])

    return AC_list, AC_adj_list, AF_list


def fill_empty_fields(line_parts, alt_list, d):
    for i in range(len(alt_list)):
        AC_list, AC_adj_list, AF_list = update_main_fields(line_parts, d)
        d['AC'].append(AC_list[i])
        d['AC_ADJ'].append(AC_adj_list[i])
        d['AF'].append(AF_list[i])
        d['ALT'].append(alt_list[i])

        for k in CSQs:
            d[k].append('')


if __name__ == '__main__':

    # The 0-indexed positions of desired CSQs in the data - will be determined once we parse the metadata
    CSQ_I = {}

    _open = gzip.open if INPUT_FILE.endswith('.gz') else open

    for chromosome, output_file in zip(CHROMOSOMES, OUTPUT_FILES):
        with _open(INPUT_FILE, 'rt') as vcf_file:

            info_list = []
            d = defaultdict(list)

            for line_no, line in enumerate(vcf_file):

                if line.startswith('##INFO='):
                    line = line[7:].strip().lstrip('<').rstrip('>')
                    _id, _desc = re.match(r'ID=(\w+),.*Description=\"(.*)\"', line).groups()
                    info_list.append({'ID': _id, 'Description': _desc})

                    continue

                elif line.startswith('#'):
                    continue

                if info_list and not CSQ_I:
                    # Converting INFO entries to a DataFrame for ease of use, but not saving yet
                    info_df = pd.DataFrame(info_list).set_index('ID')
                    assert 'CSQ' in info_df.index, \
                        'CSQ key not found in INFO. Regenerate .vcf using --vcf_info_field CSQ'
                    try:
                        _desc = info_df.loc['CSQ'].Description
                        csq_format = _desc[_desc.index('Format: ') + 8:].strip('\"\'').upper()
                        csq_fields = csq_format.split('|')
                        CSQ_I = {CSQ: csq_fields.index(CSQ) for CSQ in CSQs}
                    except (IndexError, ValueError):
                        raise RuntimeError('Unable to determine CSQ positions for desired CSQ fields')

                line_parts = line.split('\t')

                if line_parts[CHROM] != chromosome:
                    if d:
                        break  # Chromosome entry changed, and we have data for this chromosome - no need to continue
                    else:
                        continue

                if REMOVE_MULTI_ALLELIC:
                    pos = int(line_parts[POS])
                    for _start, _end in MULTI_ALLELIC_REGIONS.get(chromosome, []):
                        if _start <= pos <= _end:
                            continue

                alt_list = line_parts[ALT].split(',')
                info = line_parts[INFO]

                CSQ_start = info.find('CSQ=')
                if CSQ_start == -1:
                    fill_empty_fields(line_parts, alt_list, d)
                else:
                    CSQ_end = info.find(';', CSQ_start)
                    if CSQ_end == -1:
                        CSQ_end = None
                    CSQ_features = info[CSQ_start + 4: CSQ_end].split(',')

                    for CSQ in CSQ_features:
                        CSQ_PARTS = CSQ.split('|')
                        # Update the main fields for each CSQ feature (so each CSQ will appear in a different line)
                        AC_list, AC_adj_list, AF_list = update_main_fields(line_parts, d)

                        # Allele_num for deciding which alt, AC and AF to add - 1-indexed in .vcf file
                        allele_num = int(CSQ_PARTS[CSQ_I['ALLELE_NUM']]) - 1
                        assert allele_num >= 0, 'Unexpected Condition'

                        d['ALT'].append(alt_list[allele_num])
                        d['AC'].append(AC_list[allele_num])
                        d['AC_ADJ'].append(AC_adj_list[allele_num])
                        d['AF'].append(AF_list[allele_num])

                        for k in CSQs:
                            d[k].append(CSQ_PARTS[CSQ_I[k]])

            pd.DataFrame(d).to_csv(output_file, sep='\t')
