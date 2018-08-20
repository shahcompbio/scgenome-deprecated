from setuptools import setup

setup(name='scgenome',
      version='0.0.1',
      description='Code for analyzing single cell genomes',
      author='Shah Lab',
      url='https://www.shahlab.ca/',
      packages=['scgenome'],
      entry_points={
        'console_scripts': [
            'reduce_copynumber_dims = scgenome.reduce_copynumber_dims:main',
            'cluster_cell_embeddings = scgenome.cluster_cell_embeddings:main',
            'filter_copynumber_cells = scgenome.filter_copynumber_cells:main',
            'filter_copynumber_bins = scgenome.filter_copynumber_bins:main',
            'filter_copynumber_contiguous_duplicate_bins = '
            'scgenome.filter_copynumber_contiguous_duplicate_bins:main',
        ]
      })
