import os
import os.path
import glob
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

from dsprint.conservation_score import ConservationScore

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 5:
        print('Usage: <script> <input_pik_folder> <chrom_score_gz_pattern> <cons_name> <output_pik_folder>')
        sys.exit(0)

    INPUT_PIK_FOLDER, CHROM_SCORE_GZ_PATTERN, CONSERVATION_NAME, OUTPUT_PIK_FOLDER = sys.argv[1:]
else:
    CONSERVATION_NAME = snakemake.params.conservation_name
    CHROM_SCORE_GZ_PATTERN = snakemake.params.chrom_score_gz_pattern
    INPUT_PIK_FOLDER, = snakemake.input
    OUTPUT_PIK_FOLDER, = snakemake.output


def add_conservation_score(input_pik, conservation, output_folder):

    domain = os.path.basename(input_pik).split('.')[0]
    with open(input_pik, 'rb') as f:
        states_dict = pickle.load(f)

    all_scores = []
    for state, ds in states_dict.items():
        # underscored keys indicate an attribute applicable to the whole domain, not a numerical 'state'
        if str(state).startswith('_'):
            continue
        for d in ds:
            chromosome = d['chrom']
            scores = []
            for position in d['chrom_pos']:
                try:
                    score = conservation[chromosome][position]
                except (KeyError, IndexError):
                    pass
                else:
                    scores.append(score)
                    all_scores.append(score)
            d[CONSERVATION_NAME] = scores

    states_dict[f'_{CONSERVATION_NAME}_mean'] = np.mean(all_scores)
    states_dict[f'_{CONSERVATION_NAME}_std'] = np.std(all_scores)

    with open(os.path.join(output_folder, f'{domain}.pik'.format(domain=domain)), 'wb') as f:
        pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':

    if not os.path.exists(OUTPUT_PIK_FOLDER):
        os.mkdir(OUTPUT_PIK_FOLDER)

    conservation_score = ConservationScore(CHROM_SCORE_GZ_PATTERN)
    conservation_score.stage1()

    with ThreadPoolExecutor(max_workers=1) as executor:
        futures = {
            executor.submit(add_conservation_score, input_pik, conservation_score, OUTPUT_PIK_FOLDER):
                input_pik
            for input_pik in glob.glob(f'{INPUT_PIK_FOLDER}/*.pik')
        }
        for future in as_completed(futures):
            print(f'Finished processing {futures[future]}')