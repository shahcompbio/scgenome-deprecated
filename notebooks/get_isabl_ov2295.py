""" Scripts for generating a list of paths for OV2295 from isabl

Requirements:
    - pandas
    - pyyaml
    - shahlabdata
"""

import shahlabdata.isabl
import pandas as pd
import yaml

applications = [
    'SCDNA-ALIGNMENT',
    'SCDNA-HMMCOPY',
    'SCDNA-ANNOTATION',
    'SCDNA-COUNTHAPS',
    'SCDNA-BREAKPOINTCALLING',
    'SCDNA-VARIANTCALLING',
    'SCDNA-SNVGENOTYPING',
]

info = pd.concat([shahlabdata.isabl.get_results(a, patient_id='OV2295') for a in applications])

info = info.query('result_type == "metadata"')

info['analysis_type'] = info['result_filepath'].apply(lambda f: yaml.safe_load(open(f))['meta']['type'])

info.set_index(['isabl_patient_id', 'isabl_sample_id', 'analysis_type'])[['result_filepath']].to_dict()

info = (info
    .groupby('isabl_sample_id')[['analysis_type', 'result_filepath', 'version']]
    .apply(lambda x: x.set_index('analysis_type').to_dict(orient='index'))
    .to_dict())

with open('OV2295.yaml', 'w') as f:
    yaml.dump(info, f)
