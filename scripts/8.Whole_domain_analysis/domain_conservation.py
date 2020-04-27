from collections import defaultdict
import os.path
import glob
import numpy as np
import pickle

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <hmm_states_folder> <output_pik_file>')
        sys.exit(0)

    HMM_STATES_FOLDER, OUTPUT_FILE = sys.argv[1:]
else:
    HMM_STATES_FOLDER = snakemake.input
    OUTPUT_FILE = snakemake.output


if __name__ == '__main__':

    conservation_dict = defaultdict(dict)

    for domain_file in glob.glob(HMM_STATES_FOLDER + '/*.pik'):

        domain = os.path.splitext(os.path.basename(domain_file))[0]

        with open(domain_file, 'rb') as handle:
            states_dict = pickle.load(handle)

        phyloP_list = []
        phastCons_list = []

        for state in states_dict:
            for d in states_dict[state]:
                for score in d["phyloP"]:
                    phyloP_list.append(score)
                for score in d["phastCons"]:
                    phastCons_list.append(score)

        conservation_dict[domain]["phyloP_mean"] = np.mean(phyloP_list)
        conservation_dict[domain]["phyloP_std"] = np.std(phyloP_list)
        conservation_dict[domain]["phastCons_mean"] = np.mean(phastCons_list)
        conservation_dict[domain]["phastCons_std"] = np.std(phastCons_list)

    with open(OUTPUT_FILE, 'wb') as f:
        pickle.dump(conservation_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
