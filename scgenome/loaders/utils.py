import os
import yaml


def find_filenames(filenames, suffix):
    return [f for f in filenames if f.endswith(suffix)]


def _find_manifest_filenames(
        ticket_id,
        local_cache_directory,
    ):
    ticket_directory = os.path.join(local_cache_directory, ticket_id)

    for dirpath, dirnames, filenames in os.walk(ticket_directory):
        for filename in filenames:
            if not filename == 'metadata.yaml':
                continue

            manifest_filename = os.path.join(dirpath, filename)

            yield manifest_filename


def find_results_directories(
        ticket_id,
        local_cache_directory,
    ):
    results_directories = {}

    for manifest_filename in _find_manifest_filenames(ticket_id, local_cache_directory):
        manifest = yaml.load(open(manifest_filename))

        results_type = manifest['meta']['type']

        # KLUDGE: alignment -> align
        if results_type == 'alignment':
            results_type = 'align'

        # KLUDGE: 0.3.1 -> v0.3.1
        if not manifest['meta']['version'].startswith('v'):
            manifest['meta']['version'] = 'v' + manifest['meta']['version']

        if results_type in results_directories:
            raise Exception('found {} and {} with results type {}'.format(
                results_directories[results_type],
                os.path.dirname(manifest_filename),
                results_type))

        results_directories[results_type] = os.path.dirname(manifest_filename)

    return results_directories


def get_version(results_dir):
    """ Get version for a given results.
    
    Args:
        results_dir (str): pseudobulk results directory
    
    Returns:
        str: results version
    """
    manifest_filename = os.path.join(results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))
    return manifest['meta']['version']


def get_pseudobulk_files(results_dir, suffix):
    """ Get files for libraries and samples by suffix
    
    Args:
        results_dir (str): pseudobulk results directory
        suffix (str): suffix of requested files
    
    Yields:
        (str, str, str): sample id, library id, filename
    """
    manifest_filename = os.path.join(results_dir, 'metadata.yaml')
    manifest = yaml.load(open(manifest_filename))

    tumour_samples = manifest['meta']['tumour_samples']
    filenames = manifest['filenames']

    for sample_info in tumour_samples:
        sample_id = sample_info['sample_id']
        library_id = sample_info['library_id']

        sample_lib_suffix = f'{sample_id}_{library_id}_{suffix}'
        sample_lib_filenames = list(filter(lambda a: a.endswith(sample_lib_suffix), filenames))

        if len(sample_lib_filenames) != 1:
            raise ValueError(f'found {len(sample_lib_filenames)} {suffix} files for {sample_id}, {library_id}, {results_dir}')

        sample_lib_filename = sample_lib_filenames[0]
        sample_lib_filepath = os.path.join(results_dir, sample_lib_filename)

        yield sample_id, library_id, sample_lib_filepath


