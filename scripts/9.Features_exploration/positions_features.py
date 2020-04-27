import pandas as pd
import numpy as np
import pickle
from collections import defaultdict
import math

from dnds_func import seq_ns
from aa_chemical_properties import aa_charge, aa_charge_dict, aa_functional_group, aa_functional_group_dict, hindex_Kyte_Doolitle, aa_propensity,\
                                    propensity_chou_fasman, aa_volume_group, aa_volume, aa_volume_group_dict, aa_h_bond_donor, aa_h_bond_acceptor
from ext_predictors_codes import sift_codes, polyphen_codes, clinvar_codes
from calc_exac_freq_func import codon_table
from entropy_func import SE_hist, JSD_background, JSD_hist
from go_groups import go_term_group


def ExAC_MAF_features(features_dict, state_id, table_columns, sites_aa_num, sites_aa_alter_num, maf_list):
    # Feature: avg MAF
    if (sites_aa_num == 0):
        avg_maf_overall = 0
    else:
        avg_maf_overall = np.sum(maf_list) / float(sites_aa_num)
    features_dict[state_id].append(avg_maf_overall)
    table_columns.append("avg_maf_all")

    # Feature: avg MAF of all the altered sites
    if sites_aa_alter_num == 0:
        avg_maf_only_altered = 0
    else:
        avg_maf_only_altered = np.sum(maf_list) / float(sites_aa_alter_num)
    features_dict[state_id].append(avg_maf_only_altered)
    table_columns.append("avg_maf_altered")

    bins = [0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.5]
    non_zero_maf_lst = np.array(maf_list)[np.nonzero(maf_list)[0].tolist()]

    maf_hist = np.histogram(non_zero_maf_lst, bins)[0]

    features_dict[state_id].extend(maf_hist)
    for i in range(len(bins) - 1):
        hist_col_title = "maf_hist_" + str(bins[i]) + "-" + str(bins[i + 1])
        table_columns.append(hist_col_title)


def ExAC_population_features(features_dict, state_id, table_columns, ac_sum, ac_sum_syn, ac_sum_nonsyn,
                             an_list, pop_maf_list, pop_maf_syn_list, pop_maf_nonsyn_list):

    # Feature: populations total maf avg
    for i in range(len(an_str)):
        if (len(pop_maf_list[i]) == 0):
            avg_pop_maf = 0
        else:
            avg_pop_maf = np.average(pop_maf_list[i])
        features_dict[state_id].append(avg_pop_maf)
        table_columns.append("maf_" + an_str[i][3:])

    # Feature: populations syn maf avg
    for i in range(len(an_str)):
        if (len(pop_maf_syn_list[i]) == 0):
            avg_pop_maf_syn = 0
        else:
            avg_pop_maf_syn = np.average(pop_maf_syn_list[i])
        features_dict[state_id].append(avg_pop_maf_syn)
        table_columns.append("maf_syn_" + an_str[i][3:])

    # Feature: populations non-syn maf avg
    for i in range(len(an_str)):
        if (len(pop_maf_nonsyn_list[i]) == 0):
            avg_pop_maf_nonsyn = 0
        else:
            avg_pop_maf_nonsyn = np.average(pop_maf_nonsyn_list[i])
        features_dict[state_id].append(avg_pop_maf_nonsyn)
        table_columns.append("maf_nonsyn_" + an_str[i][3:])


def ExAC_count_features(features_dict, state_id, table_columns, sites_aa_num, sites_aa_alter_num,
                        sites_snp_num, sites_snp_alter_num):
    # Feature: number of alterations - aa level (raw and normalized by total number of matched positions)
    if (sites_aa_num == 0):
        norm_aa_alter_num = 0
    else:
        norm_aa_alter_num = sites_aa_alter_num / float(sites_aa_num)
    features_dict[state_id].append(sites_aa_alter_num)
    table_columns.append("alter_num_aa")
    features_dict[state_id].append(norm_aa_alter_num)
    table_columns.append("alter_num_aa_norm")

    # Feature: number of alterations - DNA level (raw and normalized by total number of matched positions)
    if (sites_snp_num == 0):
        norm_snp_alter_num = 0
    else:
        norm_snp_alter_num = sites_snp_alter_num / float(sites_snp_num)
    features_dict[state_id].append(sites_snp_alter_num)
    table_columns.append("alter_num_snp")
    features_dict[state_id].append(norm_snp_alter_num)
    table_columns.append("alter_num_snp_norm")

    # Feature: average number of poymorphisms at one site
    if sites_aa_alter_num == 0:
        avg_poly_aa = 0
    else:
        avg_poly_aa = sites_poly_aa_num / float(sites_aa_alter_num)
    features_dict[state_id].append(avg_poly_aa)
    table_columns.append("avg_aa_polymorphisms")

    # Feature: fraction of altered sites with more than 1 polymorphism
    if sites_aa_alter_num == 0:
        frac_poly_several = 1
    else:
        frac_poly_several = sites_poly_aa_several / float(sites_aa_alter_num)
    features_dict[state_id].append(frac_poly_several)
    table_columns.append("frac_poly_aa")


def ExAC_rareSNP_features(features_dict, state_id, table_columns, sites_snp_alter_num, rare_5_num,
                          rare_05_num, rare_005_num):
    # Feature: fraction of rare SNPs (0.5%, 0.05%, 0.005%)
    if (sites_snp_alter_num == 0):
        frac_rare_5 = 0
        frac_rare_05 = 0
        frac_rare_005 = 0
    else:
        frac_rare_5 = rare_5_num / float(sites_snp_alter_num)
        frac_rare_05 = rare_05_num / float(sites_snp_alter_num)
        frac_rare_005 = rare_005_num / float(sites_snp_alter_num)

    features_dict[state_id].append(frac_rare_5)
    table_columns.append("rare_poly_0.5")
    features_dict[state_id].append(frac_rare_05)
    table_columns.append("rare_poly_0.05")
    features_dict[state_id].append(frac_rare_005)
    table_columns.append("rare_poly_0.005")


def conservation_features(features_dict, state_id, table_columns, phastCons_dict, phyloP_dict):
    # Features: conservation scores avg for each codon position - phastCons
    features_dict[state_id].append(np.nanmean(phastCons_dict[1]))
    table_columns.append("phastCons1_avg")
    features_dict[state_id].append(np.nanmean(phastCons_dict[2]))
    table_columns.append("phastCons2_avg")
    features_dict[state_id].append(np.nanmean(phastCons_dict[3]))
    table_columns.append("phastCons3_avg")

    # Features: conservation scores avg for each codon position - phyloP
    features_dict[state_id].append(np.nanmean(phyloP_dict[1]))
    table_columns.append("phyloP1_avg")
    features_dict[state_id].append(np.nanmean(phyloP_dict[2]))
    table_columns.append("phyloP2_avg")
    features_dict[state_id].append(np.nanmean(phyloP_dict[3]))
    table_columns.append("phyloP3_avg")

    # Features: conservation scores histograms for each codon position - phastCons
    phastCons_bins = np.concatenate((np.linspace(0, 0.75, 4), np.linspace(0.8, 1.0, 5)), axis=0)
    phastCons1_hist = np.histogram(np.array(phastCons_dict[1])[~np.isnan(phastCons_dict[1])], phastCons_bins)[0]
    phastCons2_hist = np.histogram(np.array(phastCons_dict[2])[~np.isnan(phastCons_dict[2])], phastCons_bins)[0]
    phastCons3_hist = np.histogram(np.array(phastCons_dict[3])[~np.isnan(phastCons_dict[3])], phastCons_bins)[0]

    features_dict[state_id].extend(phastCons1_hist)
    features_dict[state_id].extend(phastCons2_hist)
    features_dict[state_id].extend(phastCons3_hist)
    for i in range(len(phastCons_bins) - 1):
        hist_col_title = "phastCons1_hist_" + str(phastCons_bins[i]) + "-" + str(phastCons_bins[i + 1])
        table_columns.append(hist_col_title)
    for i in range(len(phastCons_bins) - 1):
        hist_col_title = "phastCons2_hist_" + str(phastCons_bins[i]) + "-" + str(phastCons_bins[i + 1])
        table_columns.append(hist_col_title)
    for i in range(len(phastCons_bins) - 1):
        hist_col_title = "phastCons3_hist_" + str(phastCons_bins[i]) + "-" + str(phastCons_bins[i + 1])
        table_columns.append(hist_col_title)

    # Features: conservation scores histograms for each codon position - phyloP
    phyloP_bins = np.concatenate((np.array([-14, -1]), np.linspace(0, 3, 4), np.linspace(3.5, 6, 6)), axis=0)
    phyloP_hist1 = np.histogram(np.array(phyloP_dict[1])[~np.isnan(phyloP_dict[1])], phyloP_bins)[0]
    phyloP_hist2 = np.histogram(np.array(phyloP_dict[2])[~np.isnan(phyloP_dict[2])], phyloP_bins)[0]
    phyloP_hist3 = np.histogram(np.array(phyloP_dict[3])[~np.isnan(phyloP_dict[3])], phyloP_bins)[0]

    features_dict[state_id].extend(phyloP_hist1)
    features_dict[state_id].extend(phyloP_hist2)
    features_dict[state_id].extend(phyloP_hist3)
    for i in range(len(phyloP_bins) - 1):
        hist_col_title = "phyloP1_hist_" + str(phyloP_bins[i]) + "-" + str(phyloP_bins[i + 1])
        table_columns.append(hist_col_title)
    for i in range(len(phyloP_bins) - 1):
        hist_col_title = "phyloP2_hist_" + str(phyloP_bins[i]) + "-" + str(phyloP_bins[i + 1])
        table_columns.append(hist_col_title)
    for i in range(len(phyloP_bins) - 1):
        hist_col_title = "phyloP3_hist_" + str(phyloP_bins[i]) + "-" + str(phyloP_bins[i + 1])
        table_columns.append(hist_col_title)

    # Features: histogram of avg in each codon
    phastCons_codons_avg = []
    phyloP_codons_avg = []
    for i in range(len(phastCons_dict[1])):
        phastCons_score_avg = np.nanmean([phastCons_dict[1][i], phastCons_dict[2][i], phastCons_dict[3][i]])
        phastCons_codons_avg.append(phastCons_score_avg)
        phyloP_score_avg = np.nanmean([phyloP_dict[1][i], phyloP_dict[2][i], phyloP_dict[3][i]])
        phyloP_codons_avg.append(phyloP_score_avg)

    phastCons_codons_hist = np.histogram(phastCons_codons_avg, phastCons_bins)[0]
    phyloP_codons_hist = np.histogram(phyloP_codons_avg, phyloP_bins)[0]

    features_dict[state_id].extend(phastCons_codons_hist)
    features_dict[state_id].extend(phyloP_codons_hist)
    for i in range(len(phastCons_bins) - 1):
        hist_col_title = "phastCons_codons_hist_" + str(phastCons_bins[i]) + "-" + str(phastCons_bins[i + 1])
        table_columns.append(hist_col_title)
    for i in range(len(phyloP_bins) - 1):
        hist_col_title = "phyloP_codons_hist_" + str(phyloP_bins[i]) + "-" + str(phyloP_bins[i + 1])
        table_columns.append(hist_col_title)


def sub_matrix_features(features_dict, state_id, table_columns, sub_list, weigted_sub_list, sub_name):
    if len(sub_list) == 0:
        sub_avg = 0
        weigted_sub_avg = 0
        sub_postivies = 0
        sub_negatives = 0
        sub_ratio = 1
    else:
        # Feature: BLOSUM62 average and frequency weighted-average
        sub_avg = sum(sub_list) / float(len(sub_list))
        weigted_sub_avg = sum(weigted_sub_list) / float(len(weigted_sub_list))

        # Feature: BLOSUM62 count of positives and negatives
        sub_postivies = sum(1 for x in sub_list if x > 0)
        sub_negatives = sum(1 for x in sub_list if x < 0)

        # Feature: BLOSUM62 positives/negatives ratio
        if (sub_postivies == 0 or sub_negatives == 0):
            sub_ratio = 0
        else:
            sub_ratio = sub_postivies / float(sub_negatives)

    features_dict[state_id].append(sub_avg)
    table_columns.append(sub_name + "_avg")
    features_dict[state_id].append(weigted_sub_avg)
    table_columns.append(sub_name + "_avg_weighted")
    features_dict[state_id].append(sub_postivies)
    table_columns.append(sub_name + "_positive_num")
    features_dict[state_id].append(sub_negatives)
    table_columns.append(sub_name + "_negative_num")
    features_dict[state_id].append(sub_ratio)
    table_columns.append(sub_name + "_ratio")


def SIFT_features(features_dict, state_id, table_columns, sift_scores_list, weighted_sift_scores_list):
    if len(sift_scores_list) > 0:
        # Feature: SIFT average
        sift_avg = np.mean(sift_scores_list)

        # Feature: weighted (by frequency) SIFT average
        sift_w_avg = np.mean(weighted_sift_scores_list)

        # Feature: SIFT number of deleterious (score <=0.05)
        sift_deleterious_num = sum(1 for x in sift_scores_list if x <= SIFT_THRESHOLD)

        # Feature: SIFT number of tolerated (score > 0.05)
        sift_tolerated_num = sum(1 for x in sift_scores_list if x > SIFT_THRESHOLD)

        # Feature: deleterious/tolerated ratio
        if sift_tolerated_num == 0 or sift_deleterious_num == 0:
            sift_ratio = 0
        else:
            sift_ratio = sift_deleterious_num / float(sift_tolerated_num)

        # Feature: SIFT "majority-decision" (deleterious/tolerated)
        if sift_deleterious_num > sift_tolerated_num:
            sift_majority = sift_codes.SIFT_DELETERIOUS.value
        elif sift_tolerated_num > sift_deleterious_num:
            sift_majority = sift_codes.SIFT_TOLERATED.value
        else:
            sift_majority = sift_codes.SIFT_TIE.value

    else:
        sift_avg = sift_w_avg = -1
        sift_deleterious_num = 0
        sift_tolerated_num = 0
        sift_ratio = 1
        sift_majority = sift_codes.SIFT_TIE.value

    features_dict[state_id].append(sift_avg)
    table_columns.append("sift_avg")
    features_dict[state_id].append(sift_w_avg)
    table_columns.append("sift_avg_weighted")
    features_dict[state_id].append(sift_deleterious_num)
    table_columns.append("sift_deleterious_num")
    features_dict[state_id].append(sift_tolerated_num)
    table_columns.append("sift_tolerated_num")
    features_dict[state_id].append(sift_ratio)
    table_columns.append("sift_ratio")
    features_dict[state_id].append(sift_majority)
    table_columns.append("sift_majority")


def PolyPhen_features(features_dict, state_id, table_columns, polyphen_scores_list, polyphen_pred_list,
                      weighted_polyphen_scores_list):
    if (len(polyphen_scores_list) > 0):
        # Feature: PolyPhen average
        polyphen_avg = np.mean(polyphen_scores_list)

        # Feature: weighted (by frequency) PolyPhen average
        polyphen_w_avg = np.mean(weighted_polyphen_scores_list)

        # Feature: polyPhen number of benign
        polyphen_benign_num = polyphen_pred_list.count("benign")

        # Feature: polyPhen number of possibly_damaging
        polyphen_possibly_num = polyphen_pred_list.count("possibly_damaging")

        # Feature: polyPhen number of probably_damaging
        polyphen_probably_num = polyphen_pred_list.count("probably_damaging")

        # Feature: polyPhen "majority-decision" (benign/possibly_damaging/probably_damaging/unknown)
        if ((polyphen_benign_num > polyphen_probably_num and polyphen_benign_num > polyphen_possibly_num) or
                (polyphen_benign_num > polyphen_probably_num and polyphen_benign_num == polyphen_possibly_num)):
            polyphen_majority = polyphen_codes.POLYPHEN_BENIGN.value

        elif ((polyphen_probably_num > polyphen_benign_num and polyphen_probably_num > polyphen_possibly_num) or
              (polyphen_probably_num > polyphen_benign_num and polyphen_probably_num == polyphen_possibly_num)):
            polyphen_majority = polyphen_codes.POLYPHEN_PROBABLY.value

        elif (polyphen_possibly_num > polyphen_benign_num and polyphen_possibly_num > polyphen_probably_num):
            polyphen_majority = polyphen_codes.POLYPHEN_POSSIBLY.value

        elif (polyphen_benign_num == polyphen_probably_num == polyphen_possibly_num):
            polyphen_majority = polyphen_codes.PLOYPHEN_EQUAL.value

        else:
            polyphen_majority = polyphen_codes.POLYPHEN_UNKNOWN.value

    else:
        polyphen_avg = polyphen_w_avg = -1
        polyphen_benign_num = 0
        polyphen_possibly_num = 0
        polyphen_probably_num = 0
        polyphen_majority = polyphen_codes.POLYPHEN_UNKNOWN.value

    features_dict[state_id].append(polyphen_avg)
    table_columns.append("polyphen_avg")
    features_dict[state_id].append(polyphen_w_avg)
    table_columns.append("polyphen_avg_weighted")
    features_dict[state_id].append(polyphen_benign_num)
    table_columns.append("polyphen_benign_num")
    features_dict[state_id].append(polyphen_possibly_num)
    table_columns.append("polyphen_possibly_num")
    features_dict[state_id].append(polyphen_probably_num)
    table_columns.append("polyphen_probably_num")
    features_dict[state_id].append(polyphen_majority)
    table_columns.append("polyphen_majority")


def ClinVar_scores(features_dict, state_id, table_columns, clinsig_list, clinsig_af):
    valid_scores = []
    valid_scores_weighted = []

    for i in range(len(clinsig_list)):
        sig = clinsig_list[i]
        sig_list = pd.Series(sig.split("&")).unique().tolist()
        # Skipping
        if ("not" in sig_list or "" in sig_list):
            continue

        # Determine the alteration clinvar score
        if len(sig_list) == 1 and sig_list[0] == "pathogenic":
            score = clinvar_codes.CLINVAR_PATHOGENIC.value
        elif len(sig_list) == 1 and sig_list[0] == "benign":
            score = clinvar_codes.CLINVAR_BENIGN.value
        elif len(sig_list) == 2 and "benign" in sig_list and "likely" in sig_list:
            score = clinvar_codes.CLINVAR_LIKELY_BENIGN.value
        elif len(sig_list) == 2 and "pathogenic" in sig_list and "uncertain" in sig_list:
            score = clinvar_codes.CLINVAR_LIKELY_PATHOGENIC.value
        elif len(sig_list) == 2 and "pathogenic" in sig_list and "other" in sig_list:
            score = clinvar_codes.CLINVAR_PATHOGENIC_OTHER.value
        else:
            score = clinvar_codes.CLINVAR_UNCERTAIN.value  # value of 0

        valid_scores.append(score)
        score_af = clinsig_af[i]
        valid_scores_weighted.append(score * score_af)

    # ===Feature: Avg. and weighted avg. ClinVar score===#
    if len(valid_scores) == 0:
        avg_clinvar_score = 0
        avg_w_clinvar_score = 0
    else:
        avg_clinvar_score = np.mean(valid_scores)
        avg_w_clinvar_score = np.mean(valid_scores_weighted)

    features_dict[state_id].append(avg_clinvar_score)
    table_columns.append("avg_clinvar_score")
    features_dict[state_id].append(avg_w_clinvar_score)
    table_columns.append("avg_clinvar_weighted")


# Calculates a normalized Shannon entropy (from Miller et al, 2015)
def entropy(a):
    if len(a) == 1:
        return 0  # Min entropy - all the change is in one value

    a = np.asarray(a) / float(sum(a))
    entropy = 0

    for val in a:
        if val == 0 or np.isnan(val):
            continue
        if val < 0:
            print(a)
        entropy += val * math.log(val)

    entropy_adj = -entropy / math.log(len(a))  # To account for different size input

    return entropy_adj


def entropy_features(features_dict, state_id, table_columns, maf_list):
    # Feature: entropy of nonsyn SNPs distributed across instances
    if np.sum(maf_list) == 0:
        instances_entropy = math.log(len(maf_list))  # if no SNPs- each instance has the prob. = max. entropy ln(n)
    else:
        instances_entropy = entropy(maf_list)
        if (len(maf_list) == 0):
            print
            "maf_list empty"
        if (np.isnan(instances_entropy)):
            print
            "entropy nan"
            print
            maf_list

    features_dict[state_id].append(instances_entropy)
    table_columns.append("snp_nonsyn_entropy")


def pseudo_dNdS_features(features_dict, state_id, table_columns, ref_seq, Nd, Sd):
    N, S = seq_ns(ref_seq)  # Reference expected syn/nonsyn per site
    PN = 0 if N == 0 else Nd / float(N)  # Proportion of nonsyn
    PS = 0 if S == 0 else Sd / float(S)  # Proportion of syn

    # num of nonsyn substitutions per nonsyn site
    dN = -0.75 * (np.log(1 - 4 * PN / float(3)))
    features_dict[state_id].append(dN)
    table_columns.append("pseudo_nonsyn")

    # num of syn substitutions per syn site
    if 4 * PS / float(3) >= 1:
        dS = 1
    else:
        dS = -0.75 * (np.log(1 - 4 * PS / float(3)))
    features_dict[state_id].append(dS)
    table_columns.append("pseudo_syn")

    if dN == 0 or dS == 0:
        dN_dS = 1  # There isn't enough information to calculate dN/dS (1 is a neutral value)
    else:
        dN_dS = dN / dS
        if dN_dS == np.nan:
            print("dN = " + str(dN))
            print("dS = " + str(dS))
            dN_dS = 1  # There isn't enough information to calculate dN/dS (1 is a neutral value)

    features_dict[state_id].append(dN_dS)
    table_columns.append("pseudo_dNdS")


def pfam_emission_prob_features(features_dict, state_id, table_columns, domain_name, state):
    # Feature: Max. emission probability
    state_max_emiss_prob = max(hmm_prob_dict[domain_name][state])
    features_dict[state_id].append(state_max_emiss_prob)
    table_columns.append("pfam_prob_max")

    # Features: emission prob. for each amino acid
    for i in range(len(hmm_prob_dict[domain_name][state])):
        features_dict[state_id].append(hmm_prob_dict[domain_name][state][i])
        prob_aa_title = "pfam_prob_" + str(pfam_aa_order[i])
        table_columns.append(prob_aa_title)


def pfam_conserved_state_feature(features_dict, state_id, table_columns, state, con_states_dict):
    # Feature: is state is conserved according to Pfam?
    con_state = False
    if (state in con_states_dict.keys()):
        con_state = True

    features_dict[state_id].append(con_state)
    table_columns.append("is_pfam_conserved")


def instance_individuals_100way_change_features(features_dict, state_id, table_columns, maf_list,
                                                aa_ref_hist, jsd100way_list):
    # Computing Orthologus conservation in different ways (from 100way-ucsc alignment)
    # Computing Paralogus conservartion in different ways (from different instances)
    # Combining both to measurments that maximize ortho. con. and minimize para. con.

    ##Paralogus##
    # ===Feature: fraction of change across instances===#

    # determine majority aa (index of one of the majority)
    minor_counts = 0
    max_pos = aa_ref_hist.index(max(aa_ref_hist))
    for i in range(len(aa_ref_hist)):
        if i == max_pos:
            continue
        minor_counts += aa_ref_hist[i]

    instances_change_frac = minor_counts / float(np.sum(aa_ref_hist))
    features_dict[state_id].append(instances_change_frac)
    table_columns.append("instances_change_frac")

    # ===Feature: entropy of ref AA===#
    aa_ref_entropy = SE_hist(aa_ref_hist)
    features_dict[state_id].append(aa_ref_entropy)
    table_columns.append("aa_ref_SE")

    # ===Feature: JSD of ref AA===#
    aa_ref_jsd = JSD_hist(aa_ref_hist, background=JSD_background.BLOSUM62)
    features_dict[state_id].append(aa_ref_jsd)
    table_columns.append("aa_ref_jsd")

    ##Orthologus##

    # first remove -1 illegal scores of JSD mismatch (positions where JSD alignment didn't match, I added -1):
    jsd100way_list_no_mismatch = [i for i in jsd100way_list if i != -1]

    # ===Feature: median JSD score across 100way vertbrates===#
    if (len(jsd100way_list_no_mismatch) == 0):
        med_jsd = 0
    else:
        med_jsd = np.median(jsd100way_list_no_mismatch)
    features_dict[state_id].append(med_jsd)
    table_columns.append("med_jsd_100way_blosum")

    # ===Feature: Histogram of JSD score across 100way vertebrates===#
    jsd_median_bins = [0, 0.5, 0.6, 0.7, 0.8, 1]
    jsd_median_hist = np.histogram(jsd100way_list_no_mismatch, bins=jsd_median_bins)[0]

    features_dict[state_id].extend(jsd_median_hist)
    for i in range(len(jsd_median_bins) - 1):
        hist_col_title = "jsd_median_hist_" + str(jsd_median_bins[i]) + "-" + str(jsd_median_bins[i + 1])
        table_columns.append(hist_col_title)

    ##Functional measurments of both##
    # ===Feature: ratio: change across instances / change across individuals(MAF)===#
    # idea: low MAF (orthologues), high instances change (paralogous) = SDPs
    if np.sum(maf_list) == 0:
        avg_maf_overall = 0.0000001  # set the minimal non-zero in our data
    else:
        avg_maf_overall = np.sum(maf_list) / float(len(maf_list))

    instances_individuals_ratio = instances_change_frac / float(avg_maf_overall)
    features_dict[state_id].append(instances_individuals_ratio)

    table_columns.append("instances_individuals_change_ratio")

    # ===Feature: ratio: med JSD across 100way vertebrates / instances major allele freq.===#
    # idea: high JSD (orthologues), high instances change (paralogous) = SDPs
    instances_major_frac = 1 - instances_change_frac  # We want high MAF -> small 1-MAF
    jsd_instances_major_ratio = med_jsd / float(instances_major_frac)  # We want high JSD
    features_dict[state_id].append(jsd_instances_major_ratio)
    table_columns.append("jsd_100way_instances_major_ratio")


    # ===Feature: multiplication: med JSD * shannon entropy of ref aa===#
    # idea: high JSD (orthologues), high shannon entropy (paralogous) = SDPs
    jsd_mul_se = (med_jsd * aa_ref_entropy)
    features_dict[state_id].append(jsd_mul_se)
    table_columns.append("jsd_mul_aa_ref_SE")

    # ===Feature: ratio: med JSD / (max. entropy - shannon entropy of ref aa)===#
    # idea: high JSD (orthologues), low diff. of max SE to shannon entropy (paralogous) = SDPs
    max_entropy = SE_hist([0] * len(amino_acids_sym))
    entropy_diff = max_entropy - aa_ref_entropy
    jsd_SE_diff_ratio = med_jsd / float(entropy_diff)
    features_dict[state_id].append(jsd_SE_diff_ratio)
    table_columns.append("jsd_SE_diff_ratio")

    # ===Feature: sum: med JSD + normalized shannon entropy of ref aa ===#
    # idea: high JSD (orthologues), high shannon entropy (paralogous) = SDPs
    norm_SE = aa_ref_entropy / float(max_entropy)
    jsd_SE_sum = med_jsd + norm_SE
    features_dict[state_id].append(jsd_SE_sum)
    table_columns.append("jsd_SE_sum")

    # ===Feature: ratio:shannon entropy of ref aa / (max. JSD - med JSD)===#
    # idea: high shannon entropy (paralogous), low diff. of max JSD to avg JSD (orthologues) = SDPs
    max_jsd = 1
    jsd_diff = (max_jsd - med_jsd)
    SE_jsd_diff_ratio = aa_ref_entropy / float(jsd_diff)
    features_dict[state_id].append(SE_jsd_diff_ratio)
    table_columns.append("SE_jsd_diff_ratio")

    # ===Feature: ratio: med 100way-JSD / (aa ref JSD)===#
    # idea: high JSD (orthologues), low JSD (paralogoues) = SDPs
    jsds_ratio = med_jsd / float(aa_ref_jsd)
    features_dict[state_id].append(jsds_ratio)
    table_columns.append("jsds_ratio")

    # ===Feature: subtraction: (avg 100way-JSD) - (aa ref JSD)===#
    # idea: high difference between orthoulogus (more conserved) and paralogous (less conserved)
    jsds_subtraction = med_jsd - aa_ref_jsd
    features_dict[state_id].append(jsds_subtraction)
    table_columns.append("jsds_subtraction")


def aa_identity_features(features_dict, state_id, table_columns, aa_ref_hist, type_str):
    # ===Features: aa identity histogram===#
    for i in range(len(amino_acids_sym)):
        features_dict[state_id].append(aa_ref_hist[i])
        table_columns.append(type_str + "_hist_" + str(amino_acids_sym[i]))

    # ===Features: aa identity prob. vector===#
    if np.sum(aa_ref_hist) == 0:
        aa_ref_prob = aa_ref_hist
    else:
        aa_ref_prob = np.asarray(aa_ref_hist) / float(np.sum(aa_ref_hist))
    for i in range(len(amino_acids_sym)):
        features_dict[state_id].append(aa_ref_prob[i])
        table_columns.append(type_str + "_prob_" + str(amino_acids_sym[i]))


def major_allele_charge(features_dict, state_id, table_columns, aa_ref_hist):
    # ===Feature: major allele aa charge counts===#
    charge_positive_count = 0
    charge_negative_count = 0
    charge_neutral_count = 0
    for i in range(len(amino_acids_sym)):
        aa_count = aa_ref_hist[i]
        if (aa_count > 0):
            charge = aa_charge_dict[amino_acids_sym[i]]
            if (charge.value == 0):
                charge_neutral_count += aa_count
            elif (charge.value == 1):
                charge_positive_count += aa_count
            else:
                charge_negative_count += aa_count

    features_dict[state_id].append(charge_positive_count)
    table_columns.append("aa_ref_charge_positive_count")
    features_dict[state_id].append(charge_negative_count)
    table_columns.append("aa_ref_charge_negative_count")
    features_dict[state_id].append(charge_neutral_count)
    table_columns.append("aa_ref_charge_neutral_count")

    # ===Feature: major allele majority charge===#
    charge_majority = aa_charge.NEUTRAL.value
    if (charge_positive_count > charge_neutral_count and charge_positive_count > charge_negative_count):
        charge_majority = aa_charge.POSITIVE.value
    elif (charge_negative_count > charge_neutral_count and charge_negative_count > charge_positive_count):
        charge_majority = aa_charge.NEGATIVE.value

    features_dict[state_id].append(charge_majority)
    table_columns.append("aa_ref_charge_majority")


def major_allele_functional_group(features_dict, state_id, table_columns, aa_ref_hist):
    # ===Feature: major allele aa functional group counts===#
    func_counters = [0] * (len(aa_functional_group) - 1)  # Major allele is never a stop codon
    for i in range(len(amino_acids_sym)):
        aa_count = aa_ref_hist[i]
        if (aa_count > 0):
            func_group_num = aa_functional_group_dict[
                amino_acids_sym[i]].value  # getting numeric functional group value
            if (func_group_num == aa_functional_group.STOP.value):  # Major allele is never a stop codon
                continue
            func_counters[func_group_num] += aa_count

    features_dict[state_id].extend(func_counters)
    for group in aa_functional_group:
        if group == aa_functional_group.STOP:  # Major allele is never a stop codon
            continue
        func_str = "aa_ref_" + str(group) + "_count"
        table_columns.append(func_str)


def sub_diff_functional_group(features_dict, state_id, table_columns, ref_alt_pairs):
    # ===Features: count and frequency staying in functional group Vs. moving to other group===#
    stay_cnt = 0
    stay_cnt_freq = 0
    move_cnt = 0
    move_cnt_freq = 0

    for (ref, alt, af) in ref_alt_pairs:
        ref_func_group = aa_functional_group_dict[ref].value
        alt_func_group = aa_functional_group_dict[alt].value
        if (ref_func_group == alt_func_group):
            stay_cnt += 1
            stay_cnt_freq += af
        else:
            move_cnt += 1
            move_cnt_freq += af

    features_dict[state_id].append(stay_cnt)
    table_columns.append("sub_func_group_stay_cnt")
    features_dict[state_id].append(stay_cnt_freq)
    table_columns.append("sub_func_group_stay_freq")
    features_dict[state_id].append(move_cnt)
    table_columns.append("sub_func_group_move_cnt")
    features_dict[state_id].append(move_cnt_freq)
    table_columns.append("sub_func_group_move_freq")

    # ===Features: functional groups transitions counts===#
    transitions_vec_size = (len(aa_functional_group) - 1) * len(
        aa_functional_group)  # excluding transitions from STOP codons
    transitions_vec = [0] * transitions_vec_size

    for (ref, alt, af) in ref_alt_pairs:
        ref_func_group = aa_functional_group_dict[ref].value
        alt_func_group = aa_functional_group_dict[alt].value
        # Calculate counter position on the vector (ref_func_group is never STOP = 5)
        trans_vec_i = ref_func_group * (len(aa_functional_group) - 1)
        trans_vec_i += alt_func_group
        transitions_vec[trans_vec_i] += 1

    features_dict[state_id].extend(transitions_vec)
    for i in range(len(aa_functional_group) - 1):  # -1 for excluding transitions from STOP
        for j in range(len(aa_functional_group)):
            trans_col_title = "sub_func_group_trans_" + str(i) + "-" + str(j)
            table_columns.append(trans_col_title)


def major_allele_hydrophobicity(features_dict, state_id, table_columns, aa_ref_hist):
    # ===Feature: major allele hydrophicity average, hydrophobic and polar counts===#
    h_sum = 0
    h_cnt = 0
    hydrophobic_cnt = 0
    polar_charge_cnt = 0
    for i in range(len(amino_acids_sym)):
        aa_count = aa_ref_hist[i]
        if (aa_count > 0):
            hindex = hindex_Kyte_Doolitle[amino_acids_sym[i]]
            h_sum += hindex * aa_count
            h_cnt += aa_count

            if (hindex > 0):
                hydrophobic_cnt += aa_count
            else:
                polar_charge_cnt += aa_count

    h_avg = 0 if h_cnt == 0 else h_sum / float(h_cnt)

    features_dict[state_id].append(h_avg)
    table_columns.append("hindex_avg")
    features_dict[state_id].append(hydrophobic_cnt)
    table_columns.append("hindex_pos_cnt")
    features_dict[state_id].append(polar_charge_cnt)
    table_columns.append("hindex_neg_cnt")


def sub_diff_hydrophobicity(features_dict, state_id, table_columns, ref_alt_pairs):
    # ===Feature: hydrophicity difference average and weighted average===#
    hindex_diff_sum = 0
    hindex_diff_sum_weighted = 0
    hindex_diff_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        ref_hindex = hindex_Kyte_Doolitle[ref]
        alt_hindex = hindex_Kyte_Doolitle[alt]
        hindex_diff = (alt_hindex - ref_hindex)
        hindex_diff_sum += hindex_diff
        hindex_diff_sum_weighted += hindex_diff * af
        hindex_diff_cnt += 1

    if hindex_diff_cnt == 0:
        hindex_diff_avg = hindex_diff_avg_weighted = 0
    else:
        hindex_diff_avg = hindex_diff_sum / float(hindex_diff_cnt)
        hindex_diff_avg_weighted = hindex_diff_sum_weighted / float(hindex_diff_cnt)

    features_dict[state_id].append(hindex_diff_avg)
    table_columns.append("sub_diff_hindex_avg")
    features_dict[state_id].append(hindex_diff_avg_weighted)
    table_columns.append("sub_diff_hindex_avg_weighted")


def major_allele_volume(features_dict, state_id, table_columns, aa_ref_hist):
    # ===Feature: major allele volume average, tiny, small and big counts===#
    vol_sum = 0
    vol_cnt = 0
    tiny_cnt = 0
    small_cnt = 0
    big_cnt = 0
    for i in range(len(amino_acids_sym)):
        aa_count = aa_ref_hist[i]
        if (aa_count > 0):
            volume = aa_volume[amino_acids_sym[i]]
            vol_sum += volume * aa_count
            vol_cnt += aa_count

            vol_group = aa_volume_group_dict[amino_acids_sym[i]]
            if vol_group == aa_volume_group.TINY:
                tiny_cnt += aa_count
            elif vol_group == aa_volume_group.SMALL:
                small_cnt += aa_count
            elif vol_group == aa_volume_group.BIG:
                big_cnt += aa_count

    if vol_cnt == 0:
        vol_avg = 0
    else:
        vol_avg = vol_sum / float(vol_cnt)

    features_dict[state_id].append(vol_avg)
    table_columns.append("vol_avg")
    features_dict[state_id].append(tiny_cnt)
    table_columns.append("vol_tiny_cnt")
    features_dict[state_id].append(small_cnt)
    table_columns.append("vol_small_cnt")
    features_dict[state_id].append(big_cnt)
    table_columns.append("vol_big_cnt")


def sub_diff_volume(features_dict, state_id, table_columns, ref_alt_pairs):
    # ===Feature: volume difference average and weighted average===#
    volume_diff_sum = 0
    volume_diff_sum_weighted = 0
    volume_diff_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        ref_vol = aa_volume[ref]
        alt_vol = aa_volume[alt]
        vol_diff = (ref_vol - alt_vol)
        volume_diff_sum += vol_diff
        volume_diff_sum_weighted += vol_diff * af
        volume_diff_cnt += 1

    if volume_diff_cnt == 0:
        volume_diff_avg = volume_diff_avg_weighted = 0
    else:
        volume_diff_avg = volume_diff_sum / float(volume_diff_cnt)
        volume_diff_avg_weighted = volume_diff_sum_weighted / float(volume_diff_cnt)

    features_dict[state_id].append(volume_diff_avg)
    table_columns.append("sub_diff_vol_avg")
    features_dict[state_id].append(volume_diff_avg_weighted)
    table_columns.append("sub_diff_vol_avg_weighted")


def major_allele_propensity(features_dict, state_id, table_columns, aa_ref_hist):
    prop_sum = [0, 0, 0]
    prop_cnt = 0
    prop_majority_counts = [0, 0, 0]
    for i in range(len(amino_acids_sym)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            curr_prop = propensity_chou_fasman[amino_acids_sym[i]]
            mul_curr_prop = [x * aa_count for x in curr_prop]
            prop_sum = [sum(x) for x in zip(prop_sum, mul_curr_prop)]
            prop_cnt += aa_count

            if (curr_prop[aa_propensity.ALPHA_HELIX.value] == max(curr_prop)):
                prop_majority_counts[aa_propensity.ALPHA_HELIX.value] += 1
            if (curr_prop[aa_propensity.BETA_SHEET.value] == max(curr_prop)):
                prop_majority_counts[aa_propensity.BETA_SHEET.value] += 1
            if (curr_prop[aa_propensity.TURN.value] == max(curr_prop)):
                prop_majority_counts[aa_propensity.TURN.value] += 1

    # ===Feature: major allele propensity avgs===#
    if prop_cnt == 0:
        prop_avg = [0, 0, 0]
    else:
        prop_avg = [x / float(prop_cnt) for x in prop_sum]

    features_dict[state_id].extend(prop_avg)
    table_columns.extend(["aa_ref_alpha_prop_avg", "aa_ref_beta_prop_avg", "aa_ref_turn_prop_avg"])

    # ===Feature: major allele majority propensity===#
    max_idx = np.where(np.array(prop_majority_counts) == max(prop_majority_counts))[0]
    majority_vec = [0, 0, 0]
    for i in max_idx:
        majority_vec[i] = 1  # put 1 in the propensities that has max. count

    features_dict[state_id].extend(majority_vec)
    table_columns.extend(["aa_ref_alpha_is_majority", "aa_ref_beta_is_majority", "aa_ref_turn_is_majority"])


def sub_diff_propensity(features_dict, state_id, table_columns, ref_alt_pairs):
    # ===Feature: propensity difference average===#
    prop_vec_sum = [0, 0, 0]
    prop_vec_sum_weighted = [0, 0, 0]
    prop_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        ref_struct = propensity_chou_fasman[ref]
        alt_struct = propensity_chou_fasman[alt]
        prop_diff = [(x - y) for (x, y) in zip(ref_struct, alt_struct)]
        prop_diff_weighted = [(x - y) * af for (x, y) in zip(ref_struct, alt_struct)]
        prop_vec_sum = [(x + y) for (x, y) in zip(prop_vec_sum, prop_diff)]
        prop_vec_sum_weighted = [(x + y) for (x, y) in zip(prop_vec_sum_weighted, prop_diff_weighted)]

        prop_cnt += 1

    if (prop_cnt == 0):
        prop_vec_avg = prop_vec_avg_weighted = [0, 0, 0]
    else:
        prop_vec_avg = [(x / float(prop_cnt)) for x in prop_vec_sum]
        prop_vec_avg_weighted = [(x / float(prop_cnt)) for x in prop_vec_sum_weighted]

    features_dict[state_id].extend(prop_vec_avg)
    table_columns.append("sub_diff_prop_avg_alpha")
    table_columns.append("sub_diff_prop_avg_beta")
    table_columns.append("sub_diff_prop_avg_turn")
    features_dict[state_id].extend(prop_vec_avg_weighted)
    table_columns.append("sub_diff_prop_avg_alpha_weighed")
    table_columns.append("sub_diff_prop_avg_beta_weighed")
    table_columns.append("sub_diff_prop_avg_turn_weighed")


def major_allele_h_bonds(features_dict, state_id, table_columns, aa_ref_hist):
    # ===Feature: avg donor and acceptor H-bond potential===#
    donor_sum = 0
    acceptor_sum = 0
    bonds_cnt = 0
    for i in range(len(amino_acids_sym)):
        aa_count = aa_ref_hist[i]
        if (aa_count > 0):
            donor_sum += (aa_h_bond_donor[amino_acids_sym[i]] * aa_count)
            acceptor_sum += (aa_h_bond_acceptor[amino_acids_sym[i]] * aa_count)
            bonds_cnt += aa_count

    if bonds_cnt == 0:
        donor_avg = 0
        acceptor_avg = 0
    else:
        donor_avg = donor_sum / float(bonds_cnt)
        acceptor_avg = acceptor_sum / float(bonds_cnt)

    features_dict[state_id].append(donor_avg)
    table_columns.append("H_bond_donor_avg")
    features_dict[state_id].append(acceptor_avg)
    table_columns.append("H_bond_acceptor_avg")


def sub_diff_h_bonds(features_dict, state_id, table_columns, ref_alt_pairs):
    # ===Feature: acceptor and donor diff average and weighted average===#
    donor_diff_sum = 0
    donor_diff_sum_weighted = 0
    acceptor_diff_sum = 0
    acceptor_diff_sum_weighted = 0
    diff_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        ref_donor = aa_h_bond_donor[ref]
        alt_donor = aa_h_bond_donor[alt]
        donor_diff = (ref_donor - alt_donor)
        donor_diff_sum += donor_diff
        donor_diff_sum_weighted += donor_diff * af

        ref_acceptor = aa_h_bond_acceptor[ref]
        alt_acceptor = aa_h_bond_acceptor[alt]
        acceptor_diff = (ref_acceptor - alt_acceptor)
        acceptor_diff_sum += acceptor_diff
        acceptor_diff_sum += acceptor_diff * af

        diff_cnt += 1

    if diff_cnt == 0:
        donor_diff_avg = donor_diff_avg_weighted = 0
        acceptor_diff_avg = acceptor_diff_avg_weighted = 0
    else:
        donor_diff_avg = donor_diff_sum / float(diff_cnt)
        donor_diff_avg_weighted = donor_diff_sum_weighted / float(diff_cnt)
        acceptor_diff_avg = acceptor_diff_sum / float(diff_cnt)
        acceptor_diff_avg_weighted = acceptor_diff_sum_weighted / float(diff_cnt)

    features_dict[state_id].append(donor_diff_avg)
    table_columns.append("donor_diff_avg")
    features_dict[state_id].append(donor_diff_avg_weighted)
    table_columns.append("donor_diff_avg_weighted")
    features_dict[state_id].append(acceptor_diff_avg)
    table_columns.append("acceptor_diff_avg")
    features_dict[state_id].append(acceptor_diff_avg_weighted)
    table_columns.append("acceptor_diff_avg_weighted")


def spider_solvent_acc_pred(features_dict, state_id, table_columns, spider_dict):
    # ===Feature: Accessible Surface Area (solvent accessibility) average===#
    asa_avg = np.nanmean(spider_dict["spider2-ASA"])
    features_dict[state_id].append(asa_avg)
    table_columns.append("solvent_acc_avg")

    # ===Feature: Accessible Surface Area (solvent accessibility) std===#
    asa_std = np.nanstd(spider_dict["spider2-ASA"])
    features_dict[state_id].append(asa_std)
    table_columns.append("solvent_acc_std")


def spider_contact_number_pred(features_dict, state_id, table_columns, spider_dict):
    # ===Feature: contanct number for Cα-Cα average===#
    hsa2_cn_avg = np.nanmean(spider_dict["spider2-hsa2_CN"])
    features_dict[state_id].append(hsa2_cn_avg)
    table_columns.append("hsa2_cn_avg")

    # ===Feature: contanct number for Cα-Cα std===#
    hsa2_cn_std = np.nanstd(spider_dict["spider2-hsa2_CN"])
    features_dict[state_id].append(hsa2_cn_std)
    table_columns.append("hsa2_cn_std")

    # ===Feature: contanct number for Cα-Cβ average===#
    hsb2_cn_avg = np.nanmean(spider_dict["spider2-hsb2_CN"])
    features_dict[state_id].append(hsb2_cn_avg)
    table_columns.append("hsb2_cn_avg")

    # ===Feature: contanct number for Cα-Cβ std===#
    hsb2_cn_std = np.nanstd(spider_dict["spider2-hsb2_CN"])
    features_dict[state_id].append(hsb2_cn_std)
    table_columns.append("hsb2_cn_std")


def spider_angles_pred(features_dict, state_id, table_columns, spider_dict):
    # ===Feature: backbone Phi angle average===#
    Phi_angle_avg = np.nanmean(spider_dict["spider2-angle_Phi"])
    features_dict[state_id].append(Phi_angle_avg)
    table_columns.append("backbone_Phi_angle_avg")

    # ===Feature: backbone Phi angle std===#
    Phi_angle_std = np.nanstd(spider_dict["spider2-angle_Phi"])
    features_dict[state_id].append(Phi_angle_std)
    table_columns.append("backbone_Phi_angle_std")

    # ===Feature: backbone Psi angle average===#
    Psi_angle_avg = np.nanmean(spider_dict["spider2-angle_Psi"])
    features_dict[state_id].append(Psi_angle_avg)
    table_columns.append("backbone_Psi_angle_avg")

    # ===Feature: backbone Psi angle std===#
    Psi_angle_std = np.nanstd(spider_dict["spider2-angle_Psi"])
    features_dict[state_id].append(Psi_angle_std)
    table_columns.append("backbone_Psi_angle_std")

    # ===Feature: c-alpha angle (i-2=>i+1) average===#
    tau_angle_avg = np.nanmean(spider_dict["spider2-angle_tau"])
    features_dict[state_id].append(tau_angle_avg)
    table_columns.append("c-alpha_tau_angle_avg")

    # ===Feature: c-alpha angle (i-2=>i+1) std===#
    tau_angle_std = np.nanstd(spider_dict["spider2-angle_tau"])
    features_dict[state_id].append(tau_angle_std)
    table_columns.append("c-alph_tau_angle_std")

    # ===Feature: c-alpha angle (i-1=>i+1) average===#
    theta_angle_avg = np.nanmean(spider_dict["spider2-angle_theta"])
    features_dict[state_id].append(theta_angle_avg)
    table_columns.append("c-alpha_theta_angle_avg")

    # ===Feature: c-alpha angle (i-1=>i+1) std===#
    theta_angle_std = np.nanstd(spider_dict["spider2-angle_theta"])
    features_dict[state_id].append(theta_angle_std)
    table_columns.append("c-alph_theta_angle_std")


def spider_struct_pred(features_dict, state_id, table_columns, spider_dict):
    # ===Feature: helix prob. avg===#
    helix_prob_avg = np.nanmean(spider_dict["spider2-helix_prob"])
    features_dict[state_id].append(helix_prob_avg)
    table_columns.append("helix_prob_avg")

    # ===Feature: helix prob. std===#
    helix_prob_std = np.nanstd(spider_dict["spider2-helix_prob"])
    features_dict[state_id].append(helix_prob_std)
    table_columns.append("helix_prob_std")

    # ===Feature: sheet prob. avg===#
    sheet_prob_avg = np.nanmean(spider_dict["spider2-sheet_prob"])
    features_dict[state_id].append(sheet_prob_avg)
    table_columns.append("sheet_prob_avg")

    # ===Feature: sheet prob. std===#
    sheet_prob_std = np.nanstd(spider_dict["spider2-sheet_prob"])
    features_dict[state_id].append(sheet_prob_std)
    table_columns.append("sheet_prob_std")

    # ===Feature: turn prob. avg===#
    turn_prob_avg = np.nanmean(spider_dict["spider2-turn_prob"])
    features_dict[state_id].append(turn_prob_avg)
    table_columns.append("turn_prob_avg")

    # ===Feature: turn prob. std===#
    turn_prob_std = np.nanstd(spider_dict["spider2-turn_prob"])
    features_dict[state_id].append(turn_prob_std)
    table_columns.append("turn_prob_std")

    # ===Feature: major allele majority propensity===#
    struct_majority_counts = []
    struct_majority_counts.append(spider_dict["spider2-2nd_struct"].count('H'))
    struct_majority_counts.append(spider_dict["spider2-2nd_struct"].count('E'))
    struct_majority_counts.append(spider_dict["spider2-2nd_struct"].count('C'))

    max_idx = np.where(np.array(struct_majority_counts) == max(struct_majority_counts))[0]
    majority_vec = [0, 0, 0]
    for i in max_idx:
        majority_vec[i] = 1  # put 1 in the struct that has max. count

    features_dict[state_id].extend(majority_vec)
    table_columns.extend(["spd_helix_is_majority", "spd_sheet_is_majority", "spd_turn_is_majority"])


def spider_half_sphere_exposure_pred(features_dict, state_id, table_columns, spider_dict):
    # ===Feature: half-sphere exposure Cα-Cα vectors (HSEα-up) average===#
    hsa2_HSEu_avg = np.mean(spider_dict["spider2-hsa2_HSEu"])
    features_dict[state_id].append(hsa2_HSEu_avg)
    table_columns.append("hsa2_HSE-up_avg")

    # ===Feature:half-sphere exposure Cα-Cα vectors (HSEα-up) std===#
    hsa2_HSEu_std = np.std(spider_dict["spider2-hsa2_HSEu"])
    features_dict[state_id].append(hsa2_HSEu_std)
    table_columns.append("hsa2_HSE-up_std")

    # ===Feature: half-sphere exposure Cα-Cα vectors (HSEα-down) average===#
    hsa2_HSEd_avg = np.mean(spider_dict["spider2-hsa2_HSEu"])
    features_dict[state_id].append(hsa2_HSEd_avg)
    table_columns.append("hsa2_HSE-down_avg")

    # ===Feature:half-sphere exposure Cα-Cα vectors (HSEα-down) std===#
    hsa2_HSEd_std = np.std(spider_dict["spider2-hsa2_HSEd"])
    features_dict[state_id].append(hsa2_HSEd_std)
    table_columns.append("hsa2_HSE-down_std")

    # ===Feature: half-sphere exposure Cα-Cβ vectors (HSEβ-up) average===#
    hsb2_HSEu_avg = np.mean(spider_dict["spider2-hsb2_HSEu"])
    features_dict[state_id].append(hsb2_HSEu_avg)
    table_columns.append("hsb2_HSE-up_avg")

    # ===Feature:half-sphere exposure Cα-Cα vectors (HSEβ-up) std===#
    hsb2_HSEu_std = np.std(spider_dict["spider2-hsb2_HSEu"])
    features_dict[state_id].append(hsb2_HSEu_std)
    table_columns.append("hsb2_HSE-up_std")

    # ===Feature: half-sphere exposure Cα-Cβ vectors (HSEβ-down) average===#
    hsb2_HSEd_avg = np.mean(spider_dict["spider2-hsb2_HSEd"])
    features_dict[state_id].append(hsb2_HSEd_avg)
    table_columns.append("hsb2_HSE-down_avg")

    # ===Feature:half-sphere exposure Cα-Cα vectors (HSEβ-down) std===#
    hsb2_HSEd_std = np.std(spider_dict["spider2-hsb2_HSEd"])
    features_dict[state_id].append(hsb2_HSEd_std)
    table_columns.append("hsb2_HSE-down_std")


def whole_domain_conservation(features_dict, state_id, table_columns, domain_name):
    # ===Features: phastCons and PhyloP whole-domain average and std===#
    features_dict[state_id].append(whole_domain_con_dict[domain_name]["phastCons_mean"])
    table_columns.append("whole_domain_phastCons_avg")
    features_dict[state_id].append(whole_domain_con_dict[domain_name]["phastCons_std"])
    table_columns.append("whole_domain_phastCons_std")
    features_dict[state_id].append(whole_domain_con_dict[domain_name]["phyloP_mean"])
    table_columns.append("whole_domain_phyloP_avg")
    features_dict[state_id].append(whole_domain_con_dict[domain_name]["phyloP_std"])
    table_columns.append("whole_domain_phyloP_std")


def whole_domain_GO_group(features_dict, state_id, table_columns, domain_name):
    # ===Feature: domain groups (based on GO analysis)===#
    GO_groups_vec = [0] * len(go_term_group)

    if whole_domain_GO_dict.has_key(domain_name):
        domain_GO_groups = whole_domain_GO_dict[domain_name]
        for group in domain_GO_groups:
            group_val = group.value
            GO_groups_vec[group_val] = 1
    # If the domain doesn't have any of the GO terms we defined, adding +1 to NO_TERM
    else:
        GO_groups_vec[go_term_group.NO_TERM.value] = 1

    features_dict[state_id].extend(GO_groups_vec)
    for term in go_term_group:
        table_columns.append("GO:" + term.name)


def domain_location_features(features_dict, state_id, table_columns, state, max_state):
    # ===Feature: the position in the domain===#
    features_dict[state_id].append(state)
    table_columns.append("domain_pos")

    # ===Feature: the domain total length===#
    features_dict[state_id].append(max_state)
    table_columns.append("domain_length")

    # ===Feature: the location in the domain: beginning/middle/end===#
    location_list = [0, 0, 0]
    BEGIN_POS = 0
    MIDDLE_POS = 1
    END_POS = 2
    domain_location_bins = np.histogram(np.arange(1, max_state), bins=3)[1]
    if (state < domain_location_bins[1]):
        location_list[BEGIN_POS] = 1
    elif (state > domain_location_bins[2]):
        location_list[END_POS] = 1
    else:
        location_list[MIDDLE_POS] = 1

    features_dict[state_id].extend(location_list)
    table_columns.extend(["domain_pos_location_begin", "domain_pos_location_middle", "domain_pos_location_end"])


def protein_location_features(features_dict, state_id, table_columns, protein_pos_list, protein_len_list):
    # ===Feature: the protein total length (average)===#
    protein_len_avg = np.mean(protein_len_list)
    features_dict[state_id].append(protein_len_avg)
    table_columns.append("prot_avg_length")

    # ===Feature: counts of he location in the protein: beginning/middle/end===#
    location_list = [0, 0, 0]
    BEGIN_POS = 0
    MIDDLE_POS = 1
    END_POS = 2
    for i in range(len(protein_pos_list)):
        prot_location_bins = np.histogram(np.arange(1, protein_len_list[i]), bins=3)[1]
        if (protein_pos_list[i] < prot_location_bins[1]):
            location_list[BEGIN_POS] += 1
        elif (protein_pos_list[i] > prot_location_bins[2]):
            location_list[END_POS] += 1
        else:
            location_list[MIDDLE_POS] += 1

    # Normalize to ratios
    location_list_norm = np.array(location_list) / sum(location_list)

    features_dict[state_id].extend(location_list_norm)
    table_columns.extend(["prot_pos_location_begin", "prot_pos_location_middle", "prot_pos_location_end"])


if __name__ == '__main__':

    table_columns = []
    HMM_STATES_FOLDER = '/media/vineetb/t5-vineetb/dsprint/out/pfam/32/hmm_states'
    domain = 'ig'

    SIFT_THRESHOLD = 0.05

    # Rare SNP thresholds
    MAFT_5 = 0.005
    MAFT_05 = 0.0005
    MAFT_005 = 0.00005

    hmm_filename = "domains_hmm_prob_dict.pik"
    pfam_aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    amino_acids_sym = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
                       "*"]
    ligands = ["dna", "dnabase", "dnabackbone", "rna", "rnabase", "rnabackbone", "peptide", "ion", "metabolite", "sm",
               "druglike", "all"]

    with open("../substitution_matrices/BLOSUM62_dict.pik", 'rb') as f:
        blosum62_dict = pickle.load(f)

    with open("../substitution_matrices/PAM40_dict.pik", 'rb') as f:
        pam40_dict = pickle.load(f)

    # Read the whole-domain conservation dict
    with open("8.Whole_domain_analysis/pfam/32/domain_conservation_dict.pik", 'rb') as f:
        whole_domain_con_dict = pickle.load(f)

    # Read the whole_domain GO classification dict
    with open("8.Whole_domain_analysis/pfam/32/domain_go_dict.pik", 'rb') as f:
        whole_domain_GO_dict = pickle.load(f)

    # Open the HMM dict - takes some time
    with open(hmm_filename, 'rb') as f:
        hmm_prob_dict = pickle.load(f)

    features_dict = defaultdict(list)
    an_str = ["an_afr", "an_amr", "an_eas", "an_fin", "an_nfe", "an_oth", "an_sas"]
    ac_str = ["ac_afr", "ac_amr", "ac_eas", "ac_fin", "ac_nfe", "ac_oth", "ac_sas"]
    domains = [domain]

    # Randomize numbers to create id feature
    np.random.seed(0)
    random_ids = np.random.permutation(len(domains))

    # for domain_name in domains:
    for i in range(len(domains)):
        domain_name = domains[i]

        dirfiles = []
        filename = dirfiles[0]
        with open(HMM_STATES_FOLDER, f'{domain_name}.pik', 'rb') as handle:
            states_dict = pickle.load(handle)

        # Create af_adj flat dict
        states_af_adj_dict = defaultdict(list)
        for state in states_dict.keys():
            for d in states_dict[state]:
                states_af_adj_dict[state].append(d["af_adj"])

        # scale the af_dict
        states_MAF_adj_dict_scaled = defaultdict(list)
        for state in states_dict.keys():
            state_len = len(states_dict[state])
            for d in states_dict[state]:
                states_MAF_adj_dict_scaled[state].append(float(d["af_adj"] / state_len))

        # Create a dict of conserved states
        con_states_dict = {}
        con_threshold = 0.5
        for state in hmm_prob_dict[domain_name].keys():
            prob_list = hmm_prob_dict[domain_name][state]
            for i in range(len(prob_list)):
                p = prob_list[i]
                if (p > con_threshold):
                    major_allele = pfam_aa_order[i]
                    con_states_dict[state] = major_allele

        # Adding states features
        for state in states_dict.keys():

            state_id = domain_name + "_" + str(state)

            # Init counters & paramters
            maf_list = []
            sites_aa_alter_num = 0
            sites_snp_alter_num = 0
            sites_aa_num = len(states_dict[state])
            sites_snp_num = 3 * sites_aa_num
            sites_poly_aa_num = 0  # The number of different aa in all the altered sites (most are 1)
            sites_poly_aa_several = 0

            # Rare-poly-counters
            rare_5_num = 0
            rare_05_num = 0
            rare_005_num = 0

            # Conservation params
            phastCons_dict = defaultdict(list)
            phyloP_dict = defaultdict(list)
            jsd100way_list = []

            # SPIDER params
            spider_dict = defaultdict(list)

            # BLOSUM62_params
            blosum62_list = []
            weigted_blosum62_list = []

            # PAM40_params
            pam40_list = []
            weigted_pam40_list = []

            # dn/ds counters and variables
            ref_seq = ""
            Nd = 0
            Sd = 0

            # SIFT params
            sift_scores_list = []
            weighted_sift_scores_list = []

            # PolyPhen params
            polyphen_scores_list = []
            weighted_polyphen_scores_list = []
            polyphen_pred_list = []

            # clinVar params
            clinsig_list = []
            clinsig_af = []

            # Major allele params
            aa_ref_hist = [0] * len(amino_acids_sym)

            # Substitution params
            aa_alt_hist = [0] * len(amino_acids_sym)
            aa_alt_prob = [0] * len(amino_acids_sym)
            aa_alt_prob_avg = [0] * len(amino_acids_sym)
            ref_alt_pairs = []

            # protein position params
            protein_pos_list = []
            protein_len_list = []

            # Populations variables
            ac_sum = [0] * len(ac_str)
            ac_sum_syn = [0] * len(ac_str)
            ac_sum_nonsyn = [0] * len(ac_str)
            an_list = [[] for i in range(len(an_str))]
            pop_maf_list = [[] for i in range(len(an_str))]
            pop_maf_syn_list = [[] for i in range(len(an_str))]
            pop_maf_nonsyn_list = [[] for i in range(len(an_str))]

            # Iterating the state dict to get properties
            for d in states_dict[state]:

                # a list of all maf per instance
                maf_list.append(d["af_adj"])

                # Creating a position pseudo-ref sequence
                ref_codon = d["bp_ref"]
                ref_seq = ref_seq + ref_codon

                # Calculating frequency-based N/S
                bp_af_adj_dict = d["bp_af_adj_dict"]
                for alt_codon in bp_af_adj_dict.keys():
                    alt_aa = codon_table[alt_codon]
                    # syn
                    if alt_aa == d["aa_ref"]:
                        Sd += bp_af_adj_dict[alt_codon]
                    # Non-syn
                    else:
                        Nd += bp_af_adj_dict[alt_codon]

                # Major allele parameters
                aa_ref = d["aa_ref"]
                aa_ref_pos = amino_acids_sym.index(aa_ref)
                aa_ref_hist[aa_ref_pos] += 1

                # Conservation scores
                phastCons_curr_list = d["phastCons"]
                if len(phastCons_curr_list) > 0:
                    phastCons_dict[1].append(phastCons_curr_list[0])
                if len(phastCons_curr_list) > 1:
                    phastCons_dict[2].append(phastCons_curr_list[1])
                else:
                    phastCons_dict[2].append(np.nan)
                if len(phastCons_curr_list) > 2:
                    phastCons_dict[3].append(phastCons_curr_list[2])
                else:
                    phastCons_dict[3].append(np.nan)

                phyloP_curr_list = d["phyloP"]
                if len(phyloP_curr_list) > 0:
                    phyloP_dict[1].append(phyloP_curr_list[0])
                if len(phyloP_curr_list) > 1:
                    phyloP_dict[2].append(phyloP_curr_list[1])
                else:
                    phyloP_dict[2].append(np.nan)
                if len(phyloP_curr_list) > 2:
                    phyloP_dict[3].append(phyloP_curr_list[2])
                else:
                    phyloP_dict[3].append(np.nan)

                jsd100way_list.append(d["100-way-BLOSUM_JSD"])

                # SPIDER parameters (add only if exist)
                if "spider2-2nd_struct" in d:
                    spider_dict["spider2-2nd_struct"].append(d["spider2-2nd_struct"])
                    spider_dict["spider2-helix_prob"].append(float(d["spider2-helix_prob"]))
                    spider_dict["spider2-sheet_prob"].append(float(d["spider2-sheet_prob"]))
                    spider_dict["spider2-turn_prob"].append(float(d["spider2-turn_prob"]))
                    spider_dict["spider2-angle_Phi"].append(float(d["spider2-angle_Phi"]))
                    spider_dict["spider2-angle_Psi"].append(float(d["spider2-angle_Psi"]))
                    spider_dict["spider2-angle_tau"].append(float(d["spider2-angle_tau"]))
                    spider_dict["spider2-angle_theta"].append(float(d["spider2-angle_theta"]))
                    spider_dict["spider2-ASA"].append(float(d["spider2-ASA"]))
                    spider_dict["spider2-hsa2_HSEu"].append(float(d["spider2-hsa2_HSEu"]))
                    spider_dict["spider2-hsa2_HSEd"].append(float(d["spider2-hsa2_HSEd"]))
                    spider_dict["spider2-hsb2_HSEu"].append(float(d["spider2-hsb2_HSEu"]))
                    spider_dict["spider2-hsb2_HSEd"].append(float(d["spider2-hsb2_HSEd"]))
                    spider_dict["spider2-hsa2_CN"].append(float(d["spider2-hsa2_CN"]))
                    spider_dict["spider2-hsb2_CN"].append(float(d["spider2-hsb2_CN"]))

                protein_pos_list.append(d["prot_pos"])
                protein_len_list.append(d["prot_len"])

                if d["af_adj"] > 0:
                    sites_aa_alter_num += 1
                    sites_snp_alter_num += len(d["an_adj"])

                    # Number of different polymorphisms at this site
                    site_poly_num = len(d["alterations_af_adj_dict"].keys())
                    sites_poly_aa_num += site_poly_num
                    if site_poly_num > 1:
                        sites_poly_aa_several += 1

                    # Rare poly features

                    for alt_codon in bp_af_adj_dict.keys():
                        # Add to counters only nonsyn SNPs
                        if codon_table[alt_codon] != codon_table[ref_codon]:
                            if bp_af_adj_dict[alt_codon] < MAFT_005:
                                rare_005_num += 1
                                rare_05_num += 1
                                rare_5_num += 1
                            elif (bp_af_adj_dict[alt_codon] < MAFT_05):
                                rare_05_num += 1
                                rare_5_num += 1
                            elif (bp_af_adj_dict[alt_codon] < MAFT_5):
                                rare_5_num += 1

                    # Alt, BLOSUM62 and PAM40 features
                    ref = d["aa_ref"]
                    for alt in d["alterations_af_adj_dict"].keys():
                        af_adj = np.mean(d["alterations_af_adj_dict"][alt])
                        # BLOSUM
                        blosum_val = blosum62_dict[ref][alt]
                        blosum62_list.append(blosum_val)
                        weigted_blosum62_list.append(blosum_val * af_adj)
                        # PAM
                        pam_val = pam40_dict[ref][alt]
                        pam40_list.append(pam_val)
                        weigted_pam40_list.append(pam_val * af_adj)
                        # Alt aa counts
                        aa_alt_pos = amino_acids_sym.index(alt)
                        aa_alt_hist[aa_alt_pos] += 1
                        # Alt aa prob.
                        aa_alt_prob[aa_alt_pos] += af_adj
                        # ref-alt pairs
                        ref_alt_pairs.append((ref, alt, af_adj))

                    # SIFT
                    sift_list = d["SIFT"]
                    for i in range(len(sift_list)):
                        s = sift_list[i]
                        if s != "":
                            try:
                                s_af = bp_af_adj_dict[d["bp_list"][i]]
                            except:
                                # The major allele was replaced, no score available for the correct substitution
                                continue
                            sift_score = float(s[s.find("(") + 1:s.find(")")])
                            sift_scores_list.append(sift_score)
                            weighted_sift_scores_list.append(sift_score * s_af)

                    # PolyPhen
                    polyphen_list = d["PolyPhen"]
                    for i in range(len(polyphen_list)):
                        s = polyphen_list[i]
                        if s != "":
                            try:
                                s_af = bp_af_adj_dict[d["bp_list"][i]]
                            except:
                                # The major allele was replaced, no score available for the correct substitution
                                continue
                            polyphen_score = float(s[s.find("(") + 1:s.find(")")])
                            polyphen_scores_list.append(polyphen_score)
                            weighted_polyphen_scores_list.append(polyphen_score * s_af)
                            polyphen_pred_list.append(s[:s.find("(")])

                    # clinVar
                    curr_clinsig_list = d["clin_sig"]
                    for i in range(len(curr_clinsig_list)):
                        s = curr_clinsig_list[i]
                        if (s != ""):
                            try:
                                s_af = bp_af_adj_dict[d["bp_list"][i]]
                            except:
                                # The major allele was replaced, no score available for the correct substitution
                                continue
                            clinsig_list.append(s)
                            clinsig_af.append(s_af)

                    # Saving indices of syn and non-syn bps
                    syn_idx = []
                    nonsyn_idx = []
                    for i in range(len(d["bp_list"])):
                        ref_aa = d["aa_ref"]
                        alt_bp = d["bp_list"][i]
                        alt_aa = codon_table[alt_bp.upper()]
                        if (alt_aa == ref_aa):
                            syn_idx.append(i)
                        else:
                            nonsyn_idx.append(i)

                    # Summing the AC per population
                    for i in range(len(ac_str)):
                        ac = ac_str[i]
                        ac_sum[i] += sum(d[ac])
                        # Summing syn and non-syn separately
                        ac_sum_syn[i] += sum(np.array(d[ac])[syn_idx])
                        ac_sum_nonsyn[i] += sum(np.array(d[ac])[nonsyn_idx])

                    # Averaging the AN per population, to do that, gathering all an to a list
                    for i in range(len(an_str)):
                        an = an_str[i]
                        (an_list[i]).extend(d[an])

                    # Averaging the MAF per population, to do that: gathering all maf!=0 to a list
                    for i in range(len(an_str)):
                        ac = ac_str[i]
                        an = an_str[i]
                        for j in range(len(d[ac])):
                            if d[an][j] != 0:
                                pop_maf = d[ac][j] / float(d[an][j])
                                if pop_maf != 0:
                                    if j in syn_idx:
                                        pop_maf_syn_list[i].append(pop_maf)
                                    else:
                                        pop_maf_nonsyn_list[i].append(pop_maf)
                                    pop_maf_list[i].append(pop_maf)

            # ===domain_regular_features===#
            features_dict[state_id].append(domain_name)
            table_columns.append("domain_name")

            domain_idx = (domains.index(domain_name))
            features_dict[state_id].append(random_ids[domain_idx])
            table_columns.append("domain_id")

            # ===ExAC MAF Features===#
            ExAC_MAF_features(features_dict, state_id, table_columns, sites_aa_num, sites_aa_alter_num, maf_list)
            ExAC_population_features(features_dict, state_id, table_columns, ac_sum, ac_sum_syn, ac_sum_nonsyn,
                                     an_list, pop_maf_list, pop_maf_syn_list, pop_maf_nonsyn_list)
            ExAC_count_features(features_dict, state_id, table_columns, sites_aa_num, sites_aa_alter_num,
                                sites_snp_num, sites_snp_alter_num)
            ExAC_rareSNP_features(features_dict, state_id, table_columns, sites_snp_alter_num, rare_5_num,
                                  rare_05_num, rare_005_num)

            # ===Conservation scores features===#
            conservation_features(features_dict, state_id, table_columns, phastCons_dict, phyloP_dict)

            # ===Substitution matrix Features===#
            sub_matrix_features(features_dict, state_id, table_columns, blosum62_list, weigted_blosum62_list, "blosum")
            sub_matrix_features(features_dict, state_id, table_columns, pam40_list, weigted_pam40_list, "pam")

            # ===pseudo-sequence dN/dS feature===#
            pseudo_dNdS_features(features_dict, state_id, table_columns, ref_seq, Nd, Sd)

            # ===Pfam HMM-emission probabilities features===#
            pfam_emission_prob_features(features_dict, state_id, table_columns, domain_name, state)

            pfam_conserved_state_feature(features_dict, state_id, table_columns, state, con_states_dict)

            # ===SIFT score features===#
            SIFT_features(features_dict, state_id, table_columns, sift_scores_list, weighted_sift_scores_list)

            # ===Polyphen score features===#
            PolyPhen_features(features_dict, state_id, table_columns, polyphen_scores_list, polyphen_pred_list,
                              weighted_polyphen_scores_list)

            # ===ClinVar score features===#
            ClinVar_scores(features_dict, state_id, table_columns, clinsig_list, clinsig_af)

            # ===Entropy features===#
            entropy_features(features_dict, state_id, table_columns, maf_list)

            # ===instances-change to individuals-change ratios & instance-change to 100way vertbrate ratio===#
            instance_individuals_100way_change_features(features_dict, state_id, table_columns, maf_list,
                                                        aa_ref_hist, jsd100way_list)

            # ===Major allele aa chemical features===#
            aa_identity_features(features_dict, state_id, table_columns, aa_ref_hist, "aa_ref")

            major_allele_charge(features_dict, state_id, table_columns, aa_ref_hist)
            major_allele_hydrophobicity(features_dict, state_id, table_columns, aa_ref_hist)
            major_allele_volume(features_dict, state_id, table_columns, aa_ref_hist)
            major_allele_functional_group(features_dict, state_id, table_columns, aa_ref_hist)
            major_allele_propensity(features_dict, state_id, table_columns, aa_ref_hist)
            major_allele_h_bonds(features_dict, state_id, table_columns, aa_ref_hist)

            # ===Substitution features===#
            aa_identity_features(features_dict, state_id, table_columns, aa_alt_hist, "aa_alt_cnt")

            for i in range(len(amino_acids_sym)):
                if aa_alt_prob[i] > 0:
                    aa_alt_prob_avg[i] = aa_alt_prob[i] / float(aa_alt_hist[i])
            aa_identity_features(features_dict, state_id, table_columns, aa_alt_prob_avg, "aa_alt_avg_freq")

            # ===Substitution chemical features===#
            sub_diff_hydrophobicity(features_dict, state_id, table_columns, ref_alt_pairs)
            sub_diff_volume(features_dict, state_id, table_columns, ref_alt_pairs)
            sub_diff_functional_group(features_dict, state_id, table_columns, ref_alt_pairs)
            sub_diff_propensity(features_dict, state_id, table_columns, ref_alt_pairs)
            sub_diff_h_bonds(features_dict, state_id, table_columns, ref_alt_pairs)

            # ===SPIDER - secondary structure and solvent accessibility predictions baxsed on the protein sequence===#
            spider_solvent_acc_pred(features_dict, state_id, table_columns, spider_dict)
            spider_contact_number_pred(features_dict, state_id, table_columns, spider_dict)
            spider_angles_pred(features_dict, state_id, table_columns, spider_dict)
            spider_struct_pred(features_dict, state_id, table_columns, spider_dict)
            spider_half_sphere_exposure_pred(features_dict, state_id, table_columns, spider_dict)

            # ===Whole-domain aggregated features===#
            whole_domain_conservation(features_dict, state_id, table_columns, domain_name)
            whole_domain_GO_group(features_dict, state_id, table_columns, domain_name)

            # ===Position and location features===#

            domain_location_features(features_dict, state_id, table_columns, state, max(states_dict.keys()))
            protein_location_features(features_dict, state_id, table_columns, protein_pos_list, protein_len_list)

    # Exporting to data-frames table
    domains_features_df = pd.DataFrame.from_dict(features_dict, orient='index')
    domains_features_df.columns = table_columns
    domains_features_df = domains_features_df.sort_index()

    # Once for ALL domains, save df
    domains_features_df.to_csv("positions_features_fixed.csv", sep='\t')
