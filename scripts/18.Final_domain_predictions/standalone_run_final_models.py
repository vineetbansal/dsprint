import pandas as pd
import numpy as np
import pickle
from collections import defaultdict
import sys
from os import getcwd
import subprocess

from prediction_general_funcs import get_features_cols, remove_unimportant_features
from tuning_helper_functions import models_req_scaling

curr_dir = '.'


def get_pickled_model(ens_model, ens_ligand, stacked=False, ens_dir=""):
    # Get model filename
    if stacked:
        models_path = curr_dir + "/stacked_pik_models/" + ens_dir
    else:
        models_path = curr_dir + "/pik_models"

    with open(models_path + "/" + ens_ligand + "_" + ens_model + ".pik", 'rb') as handle:
        pik_model = pickle.load(handle)

    return pik_model


def predcit_using_pickeled_model(pred_dict, pik_model, classifier_method, data, stacked=False):
    # Scale the data if the classifier is one of: SVM, LG, NN
    data_index = data.index
    if classifier_method in models_req_scaling:
        cols = data.columns
        # Read the saved Scaler
        if (stacked):
            with open(curr_dir + "/stacked_pik_models/" + ens + "/scaler.pik", 'rb') as handle:
                scaler = pickle.load(handle)
        else:
            with open(curr_dir + "/pik_models/scaler.pik", 'rb') as handle:
                scaler = pickle.load(handle)
        # apply same transformation to data
        data = pd.DataFrame(scaler.transform(data))
        # Restoring indices after scaling
        data.index = data_index
        # Restoring features names
        data.columns = cols

    # Predict using the pickeled model

    probs = pik_model.predict_proba(data)
    if classifier_method == "NN":
        probs_list = probs
    else:
        probs_list = []
        for l in probs:
            probs_list.append(l[1])

    # Arrange predictions in the output dictionary
    pred_dict["idx"].extend(data_index)
    pred_dict["prob"].extend(probs_list)


def create_stacked_dataset(stacking_path, stacking_ligands, stacking_models, features_data, stacking_filename,
                           keep_original_features=True):
    df_stacking_combined = pd.DataFrame()

    for stack_ligand in stacking_ligands:
        for stack_model in stacking_models:

            # Read the stacking-1st level probs of all the ligands
            staking1_filename = stack_ligand + "_" + stack_model + "_" + stacking_filename + ".csv"
            stacking1_df = pd.read_csv(stacking_path + staking1_filename, sep='\t', index_col=0)
            stacking1_df.index = stacking1_df["idx"]
            stacking1_df.columns = ["idx", stack_model + "_" + stack_ligand + "_prob"]

            # Add to the combined df tables
            if df_stacking_combined.shape[0] == 0:
                df_stacking_combined = stacking1_df
            else:
                df_stacking_combined = pd.merge(df_stacking_combined, stacking1_df, on="idx")

    # Remving the idx column after all the merging
    df_stacking_combined.index = stacking1_df["idx"]
    del df_stacking_combined["idx"]

    # Adding the original features
    if (keep_original_features):
        df_stacking_combined = pd.concat([df_stacking_combined, features_data], axis=1)

    print("#(features) = " + str(df_stacking_combined.shape[1]))

    return df_stacking_combined


if __name__ == '__main__':

    all_models_list = ["XGB", "RF", "SVM", "Logistic", "NN"]
    ligands = ["dna", "dnabase", "dnabackbone", "rna", "rnabase", "rnabackbone", "peptide", "ion", "metabolite",
               "druglike", "sm", "all"]

    ligand = "sm"
    ens = "ALL"
    layer = "2"

    # Determine models and ligands based on ensemble type
    hyperparameters = dict()
    # all_models_list_change_order = ["SVM", "XGB", "NN", "RF", "Logistic"]
    if ens == "LIGAND":
        hyperparameters["models"] = all_models_list
        hyperparameters["ligands"] = [ligand]
        out_dir = "ligand_features_probs"
    elif ens == "MODEL":
        hyperparameters["models"] = ["XGB"]
        hyperparameters["ligands"] = ligands
        out_dir = "model_features"
    elif ens == "ALL":
        hyperparameters["models"] = all_models_list
        hyperparameters["ligands"] = ligands
        out_dir = "all_features_probs"
    else:
        hyperparameters["models"] = all_models_list
        hyperparameters["ligands"] = ligands
        out_dir = "just_probs"

    print(hyperparameters)

    features_data = pd.read_csv("../9.Features_exploration/windowed_positions_features.csv", sep='\t', index_col=0)
    print("test samples positions #: "+str(features_data.shape[0]))

    # Get list of features that we use
    features_cols = get_features_cols(features_data)
    print("# of features before removal: "+str(len(features_cols)))
    remove_unimportant_features(features_data, features_cols, update_features_cols=True)
    print("# of features after removal: "+str(len(features_cols)))

    # Filter data to just these features
    features_data = features_data.loc[:,features_cols]

    for col in features_data.columns:
        if col == "domain_name":
            continue
        nan_idx = np.where(np.isnan(features_data[col].tolist()) == True)[0]
        if len(nan_idx) > 0:
            print(col+" has NaNs")

    # ---------------------------------
    # 1st layer predictions
    # ---------------------------------
    # Get Models
    if layer == "1":
        first_layer_models = dict()
        for ens_model in hyperparameters["models"]:
            if ens_model == "NN":
                continue  # This should be run on the gpu
            for ens_ligand in hyperparameters["ligands"]:
                key_str = ens_ligand + "_" + ens_model
                first_layer_models[key_str] = get_pickled_model(ens_model, ens_ligand)

    # Predict using all the models for this ligand
    if layer == "1":
        for ligand_model_key in first_layer_models.keys():

            classifier_method = ligand_model_key.split("_")[1]
            pred_dict = defaultdict(list)
            predcit_using_pickeled_model(pred_dict, first_layer_models[ligand_model_key], classifier_method, features_data)
            pred_df = pd.DataFrame.from_dict(pred_dict)

            pred_df.to_csv(curr_dir+"/1st_level_pred/"+ligand_model_key+"_"+filename_number+".csv", sep='\t')
            print("Finished predicting using "+ligand_model_key)

    # Get the stacked model
    if layer == "2":
        second_layer_model = get_pickled_model("XGB", ligand, stacked=True, ens_dir=ens)

    # Create the 2nd layer features table
    if layer == "2":
        stacking_path = curr_dir+"/1st_level_pred/"
        stacked_features_df = create_stacked_dataset(stacking_path, hyperparameters["ligands"], hyperparameters["models"], features_data, filename_number)

    # Predict using the stacked model
    if layer == "2":
        pred_dict = defaultdict(list)
        predcit_using_pickeled_model(pred_dict, second_layer_model, "XGB", stacked_features_df, stacked=True)
        pred_df = pd.DataFrame.from_dict(pred_dict)

        pred_df.to_csv(curr_dir + "/2nd_level_pred/" + ligand + "_" + ens + "_" + filename_number + ".csv", sep='\t')
        print("Finished predicting using " + ligand + "_" + ens)
