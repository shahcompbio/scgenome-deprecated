import logging
import dbclients.tantalus
import scgenome.db.search
import pandas as pd


def load_image_feature_data(
        library_id,
        results_storage_name='singlecellresults',
):
    """ Load image features for a library
    """
    tantalus_api = dbclients.tantalus.TantalusApi()

    storage_client = tantalus_api.get_storage_client(results_storage_name)

    try:
        results = scgenome.db.search.search_image_results(tantalus_api, library_id)
    except dbclients.basicclient.NotFoundError:
        logging.exception('no image data for {}'.format(library_id))
        raise

    file_instances = list(tantalus_api.get_dataset_file_instances(
        results['id'], 'resultsdataset', results_storage_name,
        filters={'filename__endswith': 'catalog.csv'}))

    if len(file_instances) != 1:
        raise ValueError(f'found {len(file_instances)} catalog.csv files')

    file_instance = file_instances[0]
    f = storage_client.open_file(file_instance['file_resource']['filename'])
    data = pd.read_csv(f, index_col=0)
    data['library_id'] = library_id

    return data
