from setuptools import setup, find_packages

setup(name='scgenome',
      version='0.0.1',
      description='Code for analyzing single cell genomes',
      author='Shah Lab',
      url='https://www.shahlab.ca/',
      packages=find_packages(),
      entry_points={
        'console_scripts': [
        ]
      },
      package_data={'':['data/mm10.fa.fai', 'data/GRCh37-lite.fa.fai', 'dtypes/*.yaml']},
    )
