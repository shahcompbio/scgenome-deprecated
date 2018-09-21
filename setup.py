from setuptools import setup, find_packages

setup(name='scgenome',
      version='0.0.1',
      description='Code for analyzing single cell genomes',
      author='Shah Lab',
      url='https://www.shahlab.ca/',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
            'cluster_copynumber = scgenome.scripts.cluster_cn:main',
            'filter_copynumber_cells = scgenome.filter_copynumber_cells:main',
            'filter_copynumber_bins = scgenome.filter_copynumber_bins:main',
            'filter_copynumber_contiguous_duplicate_bins = scgenome.filter_copynumber_contiguous_duplicate_bins:main',
            'generate_copynumber_clone_tree_hac = scgenome.generate_copynumber_clone_tree_hac:main'
        ]
      })
