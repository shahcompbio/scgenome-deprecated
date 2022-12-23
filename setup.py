from setuptools import setup, find_packages
import versioneer

setup(
    name='scgenome',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Code for analyzing single cell whole genomes',
    author='Shah Lab',
    url='https://github.com/shahcompbio/scgenome',
    packages=find_packages(),
    install_requires=[
        'adjustText',
        'anndata',
        'bamread',
        'biopython',
        'bokeh',
        'brewer2mpl',
        'Click',
        'csverve>=0.3.1',
        'hdbscan',
        'hmmlearn',
        'ipython!=8.7.0',
        'jupyter',
        'lda',
        'matplotlib',
        'nose',
        'numba',
        'numexpr',
        'numpy',
        'oauthlib',
        'pandas',
        'pyBigWig',
        'pyfaidx',
        'pyranges',
        'pysam',
        'PyYAML',
        'scikit-learn',
        'scipy',
        'seaborn',
        'statsmodels',
        'umap-learn',
        'wgs_analysis',
        'nose',
        'biopython',
        'pypeliner',
        'joblib==1.2.0',
        'sphinx==5.1.1'
    ],
    package_data={
        'scgenome': [
            'data/*',
            'dtypes/*.yaml',
            'datasets/data/*'
        ],
    },
)
