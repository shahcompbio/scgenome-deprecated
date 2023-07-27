import pandas as pd
import pyranges as pr


class Ranges(object):

    def _convert_to_pyranges(self, data):
        assert 'chr' in data.columns
        assert 'start' in data.columns
        assert 'end' in data.columns

        pr_data = pr.PyRanges(data.rename(columns={
            'chr': 'Chromosome',
            'start': 'Start',
            'end': 'End',
        }))

        return pr_data

    def _convert_to_dataframe(self, ranges):
        data = ranges.as_df().rename(columns={
            'Chromosome': 'chr',
            'Start': 'start',
            'End': 'end',
        })

        return data

    def tile_data(self, data, binsize):
        ranges = self._convert_to_pyranges(data)

        bins = pr.gf.tile_genome(ranges, binsize)

        bins_df = self._convert_to_dataframe(bins)

        return bins_df

    def intersect_two_regions(self, data1, data2):
        a = self._convert_to_pyranges(data1)
        b = self._convert_to_pyranges(data2)

        intersect_1 = a.intersect(b)
        intersect_2 = b.intersect(a)

        intersect = pd.merge(
            self._convert_to_dataframe(intersect_1),
            self._convert_to_dataframe(intersect_2),
            on=['chr', 'start', 'end'])

        return intersect
