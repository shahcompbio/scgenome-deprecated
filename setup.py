from setuptools import setup, find_packages

setup(name='scgenome',
      version='0.0.1',
      description='Code for analyzing single cell genomes',
      author='Shah Lab',
      url='https://www.shahlab.ca/',
      packages=find_packages(),
      install_requires=[
          'click', 'numpy', 'matplotlib', 'pandas', 'adjusttext'
      ],
      entry_points={
        'console_scripts': [
            'cluster_cells = scgenome.scripts.cluster_cells:main',
            'filter_copynumber = scgenome.scripts.filter_copynumber:main',
            'generate_copynumber_clone_tree_hac = scgenome.generate_copynumber_clone_tree_hac:main'
        ]
      })
