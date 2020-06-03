#General functions and constants to use in 10.Prediction code files
import pandas as pd
import numpy as np
from os import getcwd
import pickle
import sys

#ML framework imports
from sklearn.metrics import auc, roc_auc_score, precision_recall_curve, precision_score, fbeta_score
from sklearn.preprocessing import StandardScaler

#Import from other utils files
from CV_funcs import add_domain_name_from_table_idx, calc_CV_idx_iterative
from NN_classes import curr_device

ligands = ["dna", "dnabase", "dnabackbone", "rna", "rnabase", "rnabackbone", "peptide", "ion", "metabolite", "druglike", "sm", "all"]
score_cols_suffix = ["_propensity", "_prop_th_0.1", "_prop_th_0.25", "_prop_th_0.5", "_prop_th_0.75"]

#CV splits dictionary
curr_dir = getcwd()
exac_dir = curr_dir[:curr_dir.find("ExAC")]
pfam_version = "31"
folds_num = 5

#====================================================================================================================#

def get_features_cols(features_all):
    """
    Returning a list of features column names
    """
    
    features_cols = features_all.columns.tolist()
    #removing binding scores and domain name
    for ligand in ligands:
        for suffix in score_cols_suffix:
            if ligand+suffix in features_cols:
                features_cols.remove(ligand+suffix)
    features_cols.remove("domain_name")
    
    return features_cols
#====================================================================================================================#

def remove_unimportant_features(features_table, features_cols, additional_removal_features = [], update_features_cols=False):
    """
    Removing features that aren't useful for the prediction.
    """
    
    GO_BEG1 = 425
    GO_END1 = 432
    
    #Remove domain id feature that was added just as a sanity check
    features_for_removal = ["domain_id"]
    
    #Remove GO terms as these features are missing got most domains and incomplete
    GO_features = features_cols[GO_BEG1:GO_END1]
    features_for_removal.extend(GO_features)
    
    #Removing also features from the input
    features_for_removal.extend(additional_removal_features)
    
    for feature in features_for_removal:
        del features_table[feature]   
    
    #Remove the features from the featues_cols list
    if (update_features_cols):
        for feature in features_for_removal:
                features_cols.remove(feature)
#====================================================================================================================#

def compute_per_domain_auc(y_test, pred_probs, domain_pred_dict, pred_idx, classifier):
    """
    Compute the average per_domain auc and auprc for the test set
    """
    
    y_test_copy = y_test.copy(deep=True)
    y_test_copy["pred_probs"] = pred_probs
    
    domain_auc_list = []
    domain_auprc_list = []
    domain_auprc_ratio_list = []
    domain_name_list = []
    domain_pos_num_list = []
    domain_neg_num_list = []
    
    idx = y_test.index
    y_test_copy["domain_name"] = [x[:x.rfind("_")] for x in idx]
    domains_list = y_test_copy["domain_name"].unique().tolist()
        
    for domain_name in domains_list:
        
        #Get only the domain positions
        domain_df = y_test_copy[y_test_copy["domain_name"] == domain_name]

        #Find the binding and non-binding positions of this domain 
        bind_list = domain_df[domain_df["label"] == 1].index
        bind_idx = [int(x[len(domain_name)+1:]) for x in bind_list]
        bind_num = len(bind_idx)
        non_bind_list = domain_df[domain_df["label"] == 0].index
        non_bind_idx = [int(x[len(domain_name)+1:]) for x in non_bind_list]
        non_bind_num = len(non_bind_idx)
        if (bind_num == 0 or non_bind_num == 0):
            #No positions of one of the classes "binding/non-binding" - skipping"
            continue
      
        domain_pred_dict["obs"].extend(domain_df["label"])
        domain_pred_dict["prob"].extend(domain_df["pred_probs"])
        fold_list = [pred_idx] * len(domain_df["pred_probs"])
        domain_pred_dict["fold"].extend(fold_list)
        model_list = [classifier] * len(domain_df["pred_probs"])
        domain_pred_dict["model"].extend(model_list)
        domain_str_list = [domain_name] * len(domain_df["pred_probs"])
        domain_pred_dict["domain"].extend(domain_str_list)
        
        #Add number of positives and number of negatives
        domain_pos_num_list.append(bind_num)
        domain_neg_num_list.append(non_bind_num)
        #Compute domain AUC
        domain_auc = roc_auc_score(domain_df["label"], domain_df["pred_probs"])
        domain_auc_list.append(domain_auc)
        #Compute domain AUPRC
        precision, recall, thresholds = precision_recall_curve(domain_df["label"], domain_df["pred_probs"])
        domain_auprc = auc(recall, precision)
        domain_auprc_list.append(domain_auprc)
        #Add positives fraction to list
        pos_frac_ratio = bind_num/float(domain_df.shape[0])
        #Add ratio of AUPRC and positives fraction to list
        domain_auprc_ratio_list.append(domain_auprc/float(pos_frac_ratio))
        #Add domain name for AUC/AUPRC/Ratio tables
        domain_name_list.append(domain_name)
        
    #Compute the means for the lists 
    domain_auc_mean = np.mean(domain_auc_list)
    domain_auprc_mean = np.mean(domain_auprc_list)
    domain_auprc_ratio_mean = np.mean(domain_auprc_ratio_list)
    
    return (domain_auc_mean, domain_auprc_mean, domain_auprc_ratio_mean, domain_auc_list, domain_auprc_list, domain_auprc_ratio_list, domain_name_list, domain_pos_num_list, domain_neg_num_list)
#====================================================================================================================#

def area_under_precision_prob_curve(y_true, y_probs):
    
    #probs_list = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0]
    probs_list = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0]
    probs_vals = []
    precision_vals = []
    
    for prob in probs_list:
        binary_decision = [1 if x >= prob else 0 for x in y_probs]
        if (np.count_nonzero(binary_decision) == 0):
            continue
        precision_vals.append(fbeta_score(y_true, binary_decision, 0.001))
        probs_vals.append(prob)
    
    return auc(probs_vals, precision_vals)
