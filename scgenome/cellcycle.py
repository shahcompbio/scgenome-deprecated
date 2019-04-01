import pickle
import pkg_resources
import scipy.stats
import numpy as np
import pandas as pd


classifier_filename = pkg_resources.resource_filename('scgenome', 'data/cell_state_classifier')

    
def predict(cn_data):
    corr_data = []
    
    for library_id, library_cn_data in cn_data.groupby('library_id'):
        library_cn_data = library_cn_data.copy()
        library_cn_data = library_cn_data[library_cn_data['gc'] < 1.]
        library_cn_data = library_cn_data[library_cn_data['gc'] > 0.]
        library_cn_data = library_cn_data[library_cn_data['state'] < 9]

        if 'total_reads' in library_cn_data:
            library_cn_data = library_cn_data.drop('total_reads', axis=1)
        library_cn_data = library_cn_data.merge(
            library_cn_data.groupby('cell_id')['reads'].sum().rename('total_reads').reset_index())
        library_cn_data['norm_reads'] = 1e6 * library_cn_data['reads'] / library_cn_data['total_reads']
        library_cn_data = library_cn_data.query('state > 0')
        library_cn_data['norm_reads'] = library_cn_data['norm_reads'] / library_cn_data['state']
    
        #
        # Correct GC with aggregate data
        #
        agg_data = library_cn_data.groupby(['chr', 'start'])['reads'].sum().reset_index()
        agg_data = agg_data.merge(cn_data[['chr', 'start', 'gc']].drop_duplicates())
    
        z = np.polyfit(agg_data['gc'].values, agg_data['reads'].values, 3)
        p = np.poly1d(z)
    
        library_cn_data['copy2'] = library_cn_data['reads'] / p(library_cn_data['gc'].values)
    
        library_corr_data = []
        for cell_id, cell_data in library_cn_data.groupby('cell_id'):
            correlation, pvalue = scipy.stats.spearmanr(cell_data['gc'], cell_data['copy2'])
            library_corr_data.append(dict(
                correlation=correlation,
                pvalue=pvalue,
                cell_id=cell_id,
            ))
        library_corr_data = pd.DataFrame(library_corr_data)
        library_corr_data['library_id'] = library_id
    
        corr_data.append(library_corr_data)
    
    corr_data = pd.concat(corr_data, ignore_index=True)

    # For cells with low read counts, some features may
    # be null.  Mask these cells and set the s phase
    # prediction to False
    null_features = corr_data.isnull().any(axis=1)
    corr_data = corr_data.fillna(0)

    model = pickle.load(open(classifier_filename))
    
    features = [
        'correlation',
    ]
    
    X = corr_data[features].values
    y = model.predict(X)

    # Set cells with null features to False
    y[null_features.values] = False

    corr_data['is_s_phase'] = y
    
    return corr_data[['cell_id', 'is_s_phase']]

