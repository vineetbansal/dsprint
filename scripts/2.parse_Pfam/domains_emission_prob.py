import gzip
import json
import os.path
import pickle
from collections import defaultdict
import numpy as np
import dsprint

EOF = '//'

if __name__ == '__main__':

    with open(os.path.join(os.path.dirname(dsprint.__file__), 'config.json')) as json_file:
        config = json.load(json_file)

    output_folder = config['parse_pfam']['output_folder']
    for pfam_version in config['pfam']:

        pfam_hmm_file = config['pfam'][pfam_version]
        _open = gzip.open if pfam_hmm_file.endswith('.gz') else open

        pfam_file = _open(pfam_hmm_file, 'rt')

        # ----------TO FIX BELOW---------------

        domains_hmm_dict = {}
        print_flag = True

        while print_flag:
            print_flag = False
            hmm_log_prob_dict = defaultdict(list)
            for line in pfam_file:
                print_flag = True
                if line.startswith('NAME'):
                    domain_name = line[6:-1]
                if line.startswith('HMM '):
                    aa = line.split()
                    aa.remove('HMM')
                    break

            states_cnt = 1
            for line in pfam_file:
                print_flag = True
                line_list = line.split()
                if line_list[0] == EOF:
                    break
                if line_list[0] == str(states_cnt):
                    # Saving the probalities as a list for the corresponding HMM state
                    prob_strs = line.split()[1:1 + len(aa)]
                    for i in prob_strs:
                        if i == "*":
                            hmm_log_prob_dict[states_cnt].append(float('inf'))
                        else:
                            hmm_log_prob_dict[states_cnt].append(float(i))
                    states_cnt += 1

            # Checking if its the end of the entire file
            if len(hmm_log_prob_dict.keys()) == 0:
                break
            domains_hmm_dict[domain_name] = hmm_log_prob_dict

        # Saving to file
        with open("v_" + pfam_version + "_domains_hmm_dict.pik", 'wb') as handle:
            pickle.dump(domains_hmm_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Converting to probabilities, and saving a vector for each HMM match state
        domains_hmm_prob_dict = {}
        for domain_name in domains_hmm_dict.keys():
            hmm_prob_dict = {}
            for state in domains_hmm_dict[domain_name].keys():
                hmm_prob_dict[state] = 1 / np.exp(domains_hmm_dict[domain_name][state])
            domains_hmm_prob_dict[domain_name] = hmm_prob_dict

        # Saving to file
        with open("v_" + pfam_version + "_domains_hmm_prob_dict.pik", 'wb') as handle:
            pickle.dump(domains_hmm_prob_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # Reading the dictionary of HMM probabilities
        with open("v_" + pfam_version + "_domains_hmm_prob_dict.pik", 'rb') as handle:
            domains_hmm_prob_dict = pickle.load(handle)

        emission_prob = []
        for domain in domains_hmm_prob_dict.keys():
            for state in domains_hmm_prob_dict[domain]:
                emission_prob.extend(domains_hmm_prob_dict[domain][state])
