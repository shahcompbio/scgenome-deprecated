import collections
import os
import packaging.version
import yaml

import scgenome.utils


def find_filenames(filenames, suffix):
    return [f for f in filenames if f.endswith(suffix)]


def _find_manifest_filenames(results_dir):
    for dirpath, dirnames, filenames in os.walk(results_dir):
        for filename in filenames:
            if not filename == 'metadata.yaml':
                continue

            manifest_filename = os.path.join(dirpath, filename)

            yield manifest_filename


def find_results_directories(results_dir):
    results_directories = collections.defaultdict(list)

    for manifest_filename in _find_manifest_filenames(results_dir):
        manifest = yaml.safe_load(open(manifest_filename))

        # Kludge: mondrian v1 has name instead of type
        if 'type' not in manifest['meta']:
            manifest['meta']['type'] = manifest['meta']['name']

        results_type = manifest['meta']['type']

        # KLUDGE: 0.3.1 -> v0.3.1
        if not manifest['meta']['version'].startswith('v'):
            manifest['meta']['version'] = 'v' + manifest['meta']['version']

        results_directories[results_type].append(os.path.dirname(manifest_filename))

    return results_directories


def find_results_directory(results_dir, analysis_type):
    analysis_results_dirs = find_results_directories(results_dir)

    if analysis_type not in analysis_results_dirs:
        raise Exception(f'no analysis type {analysis_type} found at {results_dir}')

    if len(analysis_results_dirs[analysis_type]) > 1:
        raise Exception(f'found {len(analysis_results_dirs[analysis_type])} analyses of type {analysis_type} found at {results_dir}')

    return analysis_results_dirs[analysis_type][0]


def get_version(results_dir):
    """ Get version for a given results.
    
    Args:
        results_dir (str): pseudobulk results directory
    
    Returns:
        str: results version
    """
    manifest_filename = os.path.join(results_dir, 'metadata.yaml')
    manifest = yaml.safe_load(open(manifest_filename))
    return manifest['meta']['version']


def find_results_filepath(results_dir, filename_suffix, result_type=None, analysis_type=None):
    """ Get filepaths for libraries and samples by suffix
    
    Args:
        results_dir (str): pseudobulk results directory
        filename_suffix (str): suffix of requested files
        result_type (str): requested result type

    KwArgs:
        analysis_type (str): check analysis type

    Returns:
        str: filepath
    """

    if analysis_type is not None:
        results_dir = find_results_directory(results_dir, analysis_type)

    manifest_filename = os.path.join(results_dir, 'metadata.yaml')
    manifest = yaml.safe_load(open(manifest_filename))

    # Kludge: mondrian v1 has name instead of type
    if 'type' not in manifest['meta']:
        manifest['meta']['type'] = manifest['meta']['name']

    if analysis_type is not None and manifest['meta']['type'] != analysis_type:
        raise Exception(f"expected analysis {analysis_type} and found {manifest['meta']['type']}")

    # Two formats:
    #  - filenames for just the names of the files,
    #  - files for file information keyed by filename.
    if 'filenames' in manifest:
        filenames = list(filter(lambda a: a.endswith(filename_suffix), manifest['filenames']))

        if len(filenames) != 1:
            raise Exception(f'found {len(filenames)} {filename_suffix} files for {results_dir}: {filenames}')

        filename = filenames[0]
        filepath = os.path.join(results_dir, filename)

    elif 'files' in manifest:
        assert result_type is not None

        files = list(filter(lambda a: a[1]['result_type'] == result_type and not a[0].endswith('.yaml'), manifest['files'].items()))

        if len(files) != 1:
            raise Exception(f'found {len(files)} {filename_suffix} files for {results_dir}: {files}')

        filename = files[0][0]
        filepath = os.path.join(results_dir, filename)

    else:
        raise ValueError()

    return filepath


def get_pseudobulk_files(results_dir, suffix):
    """ Get files for libraries and samples by suffix
    
    Args:
        results_dir (str): pseudobulk results directory
        suffix (str): suffix of requested files
    
    Yields:
        (str, str, str): sample id, library id, filename
    """

    manifest_filename = os.path.join(results_dir, 'metadata.yaml')
    manifest = yaml.safe_load(open(manifest_filename))

    if packaging.version.parse(manifest['meta']['version']) < packaging.version.parse('v0.5.0'):
        for a in _get_pseudobulk_files_v_lt_050(results_dir, suffix):
            yield a
        return

    filenames = list(filter(lambda a: a.endswith(suffix), manifest['filenames']))

    if len(filenames) != 1:
        raise ValueError(f'found {len(filenames)} {suffix} files for {results_dir}: {filenames}')

    filename = filenames[0]
    filepath = os.path.join(results_dir, filename)

    yield None, None, filepath


def _get_pseudobulk_files_v_lt_050(results_dir, suffix):
    manifest_filename = os.path.join(results_dir, 'metadata.yaml')
    manifest = yaml.safe_load(open(manifest_filename))

    tumour_samples = manifest['meta']['tumour_samples']
    filenames = manifest['filenames']

    for sample_info in tumour_samples:
        sample_id = sample_info['sample_id']
        library_id = sample_info['library_id']

        sample_lib_suffix = f'{sample_id}_{library_id}_{suffix}'
        sample_lib_filenames = list(filter(lambda a: a.endswith(sample_lib_suffix), filenames))

        if len(sample_lib_filenames) != 1:
            raise ValueError(
                f'found {len(sample_lib_filenames)} {suffix} files for {sample_id}, {library_id}, {results_dir}')

        sample_lib_filename = sample_lib_filenames[0]
        sample_lib_filepath = os.path.join(results_dir, sample_lib_filename)

        yield sample_id, library_id, sample_lib_filepath


def concat_results(results):
    """ Merge a set of results tables

    Args:
        tables (list of dict of DataFrame): list of named results tables

    Returns:
        dict of DataFrame: named results tables with unified categories
    """

    names = set()
    for tables in results:
        for name in tables.keys():
            names.add(name)

    concatenated = {}
    for name in names:
        concatenated[name] = []
        for tables in results:
            if name in tables:
                concatenated[name].append(tables[name])

    for name in names:
        concatenated[name] = scgenome.utils.concat_with_categories(concatenated[name], ignore_index=True)

    return concatenated
