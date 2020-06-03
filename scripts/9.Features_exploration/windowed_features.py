import pandas as pd
import numpy as np
import pickle
from collections import defaultdict
import sys
import math
import datetime
import random
from enum import Enum


def calc_relevant_idx(pos_idx, window_size, max_idx):
    "Calculate the relevant domain positions for the window size"

    idx_list = [pos_idx]

    for i in range(1, window_size + 1):

        if ((pos_idx - i) > 0):
            idx_list.append(pos_idx - i)
        if ((pos_idx + i) <= max_idx):
            idx_list.append(pos_idx + i)

    idx_list.sort()
    return idx_list


def calc_windowed_feature(domain_name, feature_name, window_size, domains_features_df):
    """Calculate a windowed (taking the mean across the window) feature for the input domain.
    Returning a list of the windowed feature in the order of the domain positions in the input table."""

    curr_domain_table = domains_features_df[domains_features_df["domain_name"] == domain_name]
    max_pos = max([int(index[index.rfind("_") + 1:]) for index in curr_domain_table.index.tolist()])

    # init features_lists
    feature_domain_mean_values = []
    feature_domain_std_values = []

    for index, row in curr_domain_table.iterrows():
        curr_pos = int(index[index.rfind("_") + 1:])
        window_idx = calc_relevant_idx(curr_pos, window_size, max_pos)

        # Add relevant feature values from the positions in the window
        curr_pos_feature_list = []
        for pos in window_idx:
            idx = domain_name + "_" + str(pos)
            try:
                feature_val = curr_domain_table.loc[idx, :][feature_name]
                curr_pos_feature_list.append(feature_val)
            except:
                # The relevant idx isn't present, don't add it's value
                continue

        feature_domain_mean_values.append(np.mean(curr_pos_feature_list))
        feature_domain_std_values.append(np.std(curr_pos_feature_list))

    return feature_domain_mean_values, feature_domain_std_values, max_pos


if __name__ == '__main__':

    domains_features_df = pd.read_csv('positions_features_fixed.csv', sep='\t', index_col=0)
    domains_list = domains_features_df["domain_name"].unique().tolist()

    features_windows_dict = [["avg_maf_altered", [1,3,5,10]],
                             ["phyloP1_avg", [1,3,5,10]],
                             ["phyloP2_avg", [1,3,5,10]],
                             ["phyloP3_avg", [1,3,5,10]],
                             ["blosum_avg", [1,3,5,10]],
                             ["pam_avg", [1,3,5,10]],
                             ["pseudo_dNdS", [1,3,5,10]],
                             ["pfam_prob_max",[1,3,5,10]],
                             ["sift_avg", [1,3,5,10]],
                             ["polyphen_avg", [1,3,5,10]],
                             ["avg_clinvar_score", [1,3,5,10]],
                             ["med_jsd_100way_blosum", [1,3,5,10]],
                             ["jsds_ratio", [1,3,5,10]],
                             ["hindex_avg", [1,3,5,10]],
                             ["vol_avg", [1,3,5,10]],
                             ["aa_ref_charge_majority", [1,3,5,10]],
                             ["aa_ref_alpha_prop_avg", [1,3,5,10]],
                             ["aa_ref_beta_prop_avg", [1,3,5,10]],
                             ["aa_ref_turn_prop_avg" , [1,3,5,10]],
                             ["H_bond_donor_avg", [1,3,5,10]],
                             ["H_bond_acceptor_avg", [1,3,5,10]],
                             ["sub_diff_hindex_avg_weighted", [1,3,5,10]],
                             ["sub_diff_vol_avg_weighted", [1,3,5,10]],
                             ["sub_func_group_stay_freq", [1,3,5,10]],
                             ["sub_func_group_move_freq", [1,3,5,10]],
                             ["solvent_acc_avg", [1,3,5,10]],
                             ["solvent_acc_std", [1,3,5,10]],
                             ["hsa2_cn_avg", [1,3,5,10]],
                             ["hsb2_cn_avg", [1,3,5,10]],
                             ["backbone_Phi_angle_avg", [1,3,5,10]],
                             ["backbone_Psi_angle_avg", [1,3,5,10]],
                             ["c-alpha_tau_angle_avg", [1,3,5,10]],
                             ["c-alpha_theta_angle_avg", [1,3,5,10]],
                             ["helix_prob_avg", [1,3,5,10]],
                             ["sheet_prob_avg", [1,3,5,10]],
                             ["turn_prob_avg", [1,3,5,10]],
                             ["hsa2_HSE-up_avg", [1,3,5,10]],
                             ["hsa2_HSE-down_avg", [1,3,5,10]],
                             ["hsb2_HSE-up_avg", [1,3,5,10]],
                             ["hsb2_HSE-down_avg", [1,3,5,10]]]

    for i in range(len(features_windows_dict)):

        feature_name = features_windows_dict[i][0]
        for window_size in features_windows_dict[i][1]:
            feature_vals_mean_list = []
            feature_vals_std_list = []
            len_vals_list = []

            for domain_name in domains_list:
                (domain_feature_mean_vals, domain_feature_std_vals, max_pos) = calc_windowed_feature(domain_name,
                                                                                                     feature_name,
                                                                                                     window_size,
                                                                                                     domains_features_df)
                feature_vals_mean_list.extend(domain_feature_mean_vals)
                feature_vals_std_list.extend(domain_feature_std_vals)

                # Adding the windowed feature to the features table
            new_wmean_name = "wm-" + str(window_size) + "-" + feature_name
            domains_features_df[new_wmean_name] = feature_vals_mean_list
            new_wstd_name = "ws-" + str(window_size) + "-" + feature_name
            domains_features_df[new_wstd_name] = feature_vals_std_list
            print("Finished " + feature_name + " window size: " + str(window_size))

    domains_features_df.to_csv("windowed_positions_features.csv", sep='\t')
