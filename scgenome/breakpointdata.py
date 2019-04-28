import pandas as pd


def load_breakpoint_data(dataset_filepaths, sample_libraries):
    breakpoint_data = []

    for sample_id, library_id in sample_libraries:
        for filepath in dataset_filepaths:
            if filepath.endswith('{}_{}_destruct.tsv'.format(sample_id, library_id)):
                data = pd.read_csv(filepath, sep='\t')
                data['library_id'] = library_id
                data['sample_id'] = sample_id
                breakpoint_data.append(data)

    if len(breakpoint_data) == 0:
        return pd.DataFrame(), pd.DataFrame()

    breakpoint_data = pd.concat(breakpoint_data, ignore_index=True)

    breakpoint_count_data = []

    for sample_id, library_id in sample_libraries:
        for filepath in dataset_filepaths:
            if filepath.endswith('{}_{}_cell_counts_destruct.csv'.format(sample_id, library_id)):
                data = pd.read_csv(filepath)
                data['library_id'] = library_id
                data['sample_id'] = sample_id
                breakpoint_count_data.append(data)

    breakpoint_count_data = pd.concat(breakpoint_count_data, ignore_index=True)
    breakpoint_count_data = breakpoint_count_data.rename(columns={'cluster_id': 'prediction_id'})

    # KLUDGE: normal reads are not filtered properly, filter by their prefix
    breakpoint_count_data = breakpoint_count_data[~breakpoint_count_data['cell_id'].str.startswith('HS')]

    return breakpoint_data, breakpoint_count_data


