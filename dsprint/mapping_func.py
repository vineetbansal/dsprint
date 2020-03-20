import pandas as pd
import numpy as np
import pickle
import math


def correct_exons_frameshift(exon_df, targetid, frameshift_file):
    """
    Correct the exons end/start positions if there's a frameshift.
    Frameshift are extracted from Shilpa's exons sequences.
    This function corrects the indices of the exons table accordingly.
    """

    # Get the frameshifts index and length of the exons
    with open(frameshift_file, 'rb') as f:
        exons_frameshifts = pickle.load(f)

    for frameshift in exons_frameshifts[targetid + ".exons.txt"]:
        idx, length, _ = frameshift

        # Find the exon we need to add bps to
        first_bp_count = 1
        for index, exon in exon_df.iterrows():
            ex_start = int(exon[0])
            ex_end = int(exon[1])
            exon_len = (ex_end - ex_start + 1)

            # Fixing start pos of the exon
            if idx <= first_bp_count:
                exon_df.loc[index].start_pos = ex_start - length
                break
            # Fixing end pos of the exon
            elif idx <= (first_bp_count + exon_len):
                exon_df.loc[index].end_pos = ex_end + length
                break
            first_bp_count += exon_len


# -------------------------------------------------------------------------------------------#

def create_exon_pos_table(chrom_raw, targetid, frameshift_file):
    """
    A function that get chromosome raw data from the hmmer results and return a data-frame of the exons.
    """
    exons_raw = chrom_raw

    # Removing the complement bracates if exist
    if exons_raw.find("complement(") >= 0:
        exons_raw = exons_raw[exons_raw.find("complement(") + 11:-1]

    # Removing the join bracates if exist
    if exons_raw.find("join(") >= 0:
        exons_raw = exons_raw[exons_raw.find("join(") + 5:-1]

    # In case there's only one exon, take everything after the second ":"
    else:
        exons_raw = exons_raw[exons_raw.find(":", chrom_raw.find(":") + 1) + 1:]

    exons_list = exons_raw.split(",")
    exon_pos = []
    frameshift_flag = False
    for ex in exons_list:
        # flag cases where Shilpa added "-" to a position number to signify frameshift in the sequences
        if ex[0] == "-":
            frameshift_flag = True
            continue

        # Adding the real exons to exons_pos list
        exon_pos.append(ex.split(".."))

    # Creating a table for the start and end of exons
    exon_df = pd.DataFrame(exon_pos)
    exon_df.columns = ["start_pos", "end_pos"]

    # Correct frameshift if frameshift exist
    if frameshift_flag:
        correct_exons_frameshift(exon_df, targetid, frameshift_file)

    exon_len = []
    for index, exon in exon_df.iterrows():
        exon_len.append(int(exon[1]) - int(exon[0]) + 1)
    exon_df["length"] = exon_len
    first_bp_count = 1
    first_bp_list = []
    for index, exon in exon_df.iterrows():
        first_bp_list.append(first_bp_count)
        first_bp_count += int(exon[2])
    exon_df["first_bp_count"] = first_bp_list
    return exon_df


# -------------------------------------------------------------------------------------------#

def find_protein_pos(chrom_pos, exon_df, chrom_raw):
    """
    A function that get chromosome position and data-frame of exons,
    and return the protein position or -1 if it's not within any exon.
    Not currently beeing used in the pipeline.
    """
    # Go over all the exons to search for the chrom position there
    for index, exon in exon_df.iterrows():
        start_pos = int(exon[0])
        end_pos = int(exon[1])
        first_bp_count = int(exon[3])

        # If the chrom position is inside this exon
        if chrom_pos >= start_pos and chrom_pos <= end_pos:

            # Calculate position for reverse complement strand:
            # the protein is translated from the end position towards the start position of the exon
            if chrom_raw.find("complement") >= 0:
                len_from_exon_start = end_pos - chrom_pos
            # Calculate position for forward strand
            else:
                len_from_exon_start = chrom_pos - start_pos

            # Calculate the position on the mRNA transcript
            transcript_pos = len_from_exon_start + first_bp_count

            # Calculate the position on the protein sequence
            protein_pos = int(math.ceil(float(transcript_pos) / 3))

            return protein_pos

    # If the position wasn't in the regions of any exon
    return -1


# -------------------------------------------------------------------------------------------#

def find_chrom_bps(protein_pos, exon_table, chrom_raw_data):
    """
    A function that get protein position and data-frame of exons,
    and return the chromosome positions of the corresponding codon.
    """

    # calculate the mRNA transcript index of this protein position (the 1st bp in the triplet)
    transcript_pos = (protein_pos * 3) - 2

    # Iterating over all the gene exons
    for index, exon in exon_table.iterrows():
        first_bp_count = int(exon["first_bp_count"])
        exon_length = int(exon["length"])
        last_bp_count = first_bp_count + exon_length - 1

        # Checking if the transcript position is within this exon
        if first_bp_count <= transcript_pos and transcript_pos <= last_bp_count:

            start_pos = int(exon["start_pos"])
            end_pos = int(exon["end_pos"])

            len_from_exon_start = transcript_pos - first_bp_count

            # Calculate bps position for reverse complement strand:
            # the protein is translated from the end position towards the start position of the exon
            if chrom_raw_data.find("complement") >= 0:
                chrom_pos_1st = end_pos - len_from_exon_start

                chrom_pos_2nd = chrom_pos_1st - 1
                # If the exons end here: move to the next exon
                if chrom_pos_2nd < start_pos:
                    index += 1
                    chrom_pos_2nd = int(exon_table["end_pos"][index])
                    start_pos = int(exon_table["start_pos"][index])
                    end_pos = int(exon_table["end_pos"][index])

                # If the exons ends here: move to the next exon
                chrom_pos_3rd = chrom_pos_2nd - 1
                if chrom_pos_3rd < start_pos:
                    index += 1
                    chrom_pos_3rd = int(exon_table["end_pos"][index])

            # Calculate position for forward strand
            else:
                chrom_pos_1st = start_pos + len_from_exon_start

                chrom_pos_2nd = chrom_pos_1st + 1
                # If the exons end here: move to the next exon
                if chrom_pos_2nd > end_pos:
                    index += 1
                    chrom_pos_2nd = int(exon_table["start_pos"][index])
                    start_pos = int(exon_table["start_pos"][index])
                    end_pos = int(exon_table["end_pos"][index])

                # If the exons end here: move to the next exon
                chrom_pos_3rd = chrom_pos_2nd + 1
                if chrom_pos_3rd > end_pos:
                    index += 1
                    chrom_pos_3rd = int(exon_table["start_pos"][index])

            return chrom_pos_1st, chrom_pos_2nd, chrom_pos_3rd


# -------------------------------------------------------------------------------------------#

def is_number(s):
    """
    Boolean function - determine if a given text can be converted to a number
    """
    try:
        float(s)
        return True
    except ValueError:
        return False


def protein_pos_to_hmm_state_and_aa(protein_pos, domain_gene_table):
    """
    A function that return the hmm state of that protein position, and the amino acid.
    return -1 for positions outside of domains regions, -2 for matching insertion
    #TODO: do we need also transcript id? do we want to consider more than 1 transcript per gene?
    """
    for index, row in domain_gene_table.iterrows():
        target_start = row["TargetStart"]
        target_end = row["TargetEnd"]
        aa = "-"

        # Check if the position is inside this domain instance of the gene
        if protein_pos >= target_start and protein_pos <= target_end:

            hmm_pos = (row["HMM_Pos"]).split(",")
            target_seq = list(row["Target_Seq"])
            index_inside_match = int(protein_pos - target_start)

            # Get deletions indices
            indices = [i for i, x in enumerate(target_seq) if x == "-"]

            # Remove deletions from both lists
            target_seq_no_del = [i for j, i in enumerate(target_seq) if j not in indices]
            hmm_pos_no_del = [i for j, i in enumerate(hmm_pos) if j not in indices]

            # Get the aa
            aa = (target_seq_no_del[index_inside_match]).upper()

            # Find the HMM match state
            hmm_state_text = hmm_pos_no_del[index_inside_match]
            if (is_number(hmm_state_text) == True):
                hmm_state = int(hmm_state_text)
            else:
                # the position match insertion
                hmm_state = -2

            # Returning hmm_state and aa for match inside a domain's regions
            return (hmm_state, aa)

    # The protein position isn't in any domain region
    return -1, '-'
