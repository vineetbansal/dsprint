import gzip
import json
import os.path
import pandas as pd
import dsprint

# Fields that we are interested in - for multiple valued fields (e.g. GA), the first component is extracted
FIELDS = ('NAME', 'LENG', 'GA')


if __name__ == '__main__':

    with open(os.path.join(os.path.dirname(dsprint.__file__), 'config.json')) as json_file:
        config = json.load(json_file)

    output_folder = config['parse_pfam']['output_folder']
    for pfam_version in config['pfam']:

        pfam_hmm_file = config['pfam'][pfam_version]
        _open = gzip.open if pfam_hmm_file.endswith('.gz') else open

        pfam_file = _open(pfam_hmm_file, 'rt')

        d = {k: [] for k in FIELDS}
        for line in pfam_file:
            parts = line.strip().split()
            first = parts[0]
            if first in FIELDS:
                d[first].append(parts[1])

        pd.DataFrame(d).to_csv(os.path.join(output_folder, f'pfam_{pfam_version}_data.csv'), sep='\t')


