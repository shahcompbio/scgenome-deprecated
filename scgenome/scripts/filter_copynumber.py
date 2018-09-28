import click
import pandas as pd
from scgenome import cnfilter


@click.command()
@click.argument('cn_path', type=click.Path(exists=True))
@click.argument('filt_path', type=click.Path())
@click.option('--cell-scores', '-s', type=click.Path(exists=True))
@click.option('--qdnaseq-blacklist', '-q', type=click.Path(exists=True))
@click.option('--filter-contig-dup-bins', '-d', 'filt_dup', is_flag=True)
def main(cn_path, filt_path, cell_scores, qdnaseq_blacklist, filt_dup):
    """
    Filters cells and copynumber bins found in CN_PATH, outputs the results
    to FILT_PATH.
    """

    cn = pd.read_table(cn_path, dtype={'chr': str})

    if cell_scores is not None:
        scores = pd.read_csv(cell_scores)
        filtered_cn = cnfilter.filter_cells(cn, scores)
        click.echo('[+] cell score filter')
    else:
        click.echo('[ ] cell score filter')

    if qdnaseq_blacklist is not None:
        blacklist = pd.read_table(qdnaseq_blacklist, true_values=['TRUE'])
        filtered_cn = cnfilter.filter_qdnaseq_bins(filtered_cn, blacklist)
        click.echo('[+] QDNAseq blacklist filter')
    else:
        click.echo('[ ] QDNAseq blacklist filter')

    if filt_dup:
        filtered_cn = cnfilter.remove_contiguous_duplicate_bins(filtered_cn)
        click.echo('[+] duplicate contiguous bin filter')
    else:
        click.echo('[ ] duplicate contiguous bin filter')

    filtered_cn.to_csv(filt_path, sep='\t', index=False)
    return


if __name__ == '__main__':
    main()
