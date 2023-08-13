"""
input: fields cn, copy state giemsa etc
       dict of colors hex
       number of levels
        cmap name
return:
        rgb color ref
        hex color ref
        colormap




"""

from collections import defaultdict

import matplotlib
import matplotlib.cm
import matplotlib.cm
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import to_rgb, to_hex
from numpy import ndarray


class Colors(object):
    def __init__(
            self,
            field_name=None,
            levels=None,
            cmap=None,
            hex_color_reference=None,
            vmin=None,
            vmax=None
    ):

        if field_name is not None:
            self.hex_color_reference = self.infer_color_reference_from_fieldname(field_name)
            self.cmap = self.get_cmap_from_reference(vmin=vmin, vmax=vmax)

        elif levels is not None:
            cmap = self.infer_colormap_from_levels(len(levels))
            self.cmap = self.load_preset_colormap(cmap)
            self.hex_color_reference = dict(zip(levels, list(self.cmap(np.linspace(0, 1, len(levels))))))
        elif cmap is not None:
            self.cmap = self.load_preset_colormap(cmap)
            self.hex_color_reference = self.get_reference_from_cmap(vmin=vmin, vmax=vmax)
        elif hex_color_reference is not None:
            self.hex_color_reference = hex_color_reference
            self.cmap = self.get_cmap_from_reference(vmin=vmin, vmax=vmax)

        self.rgb_color_reference = self.translate_hex_to_rgb(self.hex_color_reference)

    @property
    def cn_color_reference(self):
        color_reference = {
            0: '#3182BD',
            1: '#9ECAE1',
            2: '#CCCCCC',
            3: '#FDCC8A',
            4: '#FC8D59',
            5: '#E34A33',
            6: '#B30000',
            7: '#980043',
            8: '#DD1C77',
            9: '#DF65B0',
            10: '#C994C7',
            11: '#D4B9DA',
        }

        color_reference = defaultdict(lambda: '#D4B9DA', color_reference)

        return color_reference

    @property
    def cyto_band_giemsa_stain_hex_color_reference(self):
        """
         Adapted from: https://github.com/bernatgel/karyoploteR/blob/master/R/color.R
        """
        color_reference = {
            'gneg': '#FFFFFF',
            'gpos25': '#C8C8C8',
            # 'gpos33': '#D2D2D2',
            'gpos50': '#C8C8C8',
            # 'gpos66': '#A0A0A0',
            'gpos75': '#828282',
            'gpos100': '#000000',
            'gpos': '#000000',
            'stalk': '#647FA4',  # repetitive areas
            'acen': '#D92F27',  # centromeres
            'gvar': '#DCACAC',  # previously '#DCDCDC'
        }
        return color_reference

    def infer_color_reference_from_fieldname(self, field_name):

        if field_name.lower() in ('copy', 'state', 'cn'):
            return self.cn_color_reference
        elif field_name.lower() in ('giemsa', 'cyto_band_giemsa_stain'):
            return self.cyto_band_giemsa_stain_hex_color_reference
        else:
            raise Exception()

    def infer_colormap_from_levels(self, num_levels):
        if num_levels <= 10:
            return 'tab10'
        elif num_levels <= 20:
            return 'tab20'
        else:
            return 'hsv'

    def get_cmap_from_reference(self, vmin=None, vmax=None):

        vmin = min(self.hex_color_reference.keys()) if vmin is None else vmin
        vmax = min(self.hex_color_reference.keys()) if vmax is None else vmax

        color_list = [self.hex_color_reference[cn] for cn in range(vmin, vmax + 1)]

        cmap = ListedColormap(color_list)
        return cmap

    def get_reference_from_cmap(self, vmin=None, vmax=None):
        colors = {v: to_hex(self.cmap(v)) for v in range(vmin, vmax + 1)}
        return colors

    def load_preset_colormap(self, palette):
        matplotlib_cmaps = plt.colormaps()

        if palette in matplotlib_cmaps:

            palette = matplotlib.cm.get_cmap(palette)
        else:
            palette = sns.color_palette(palette, as_cmap=True)

        return palette

    def translate_hex_to_rgb(self, color_reference):

        new_reference = {}
        for k, v in color_reference.items():
            new_reference[k] = to_rgb(v)

        if isinstance(color_reference, defaultdict):
            default_value = to_rgb(color_reference['GETDEFAULTCOLOR'])
            new_reference = defaultdict(lambda: default_value, new_reference)

        return new_reference


def map_cn_colors_to_matrix(X: ndarray) -> ndarray:
    """ Create an array of colors from an array of copy number states

    Parameters
    ----------
    X : ndarray
        copy number states

    Returns
    -------
    ndarray
        colors with shape X.shape + (3,)
    """

    color_reference = Colors(field_name='cn').rgb_color_reference

    col = []

    for idx in range(len(X)):
        col.append(np.array([color_reference[v] for v in X[idx]]))

    col = np.array(col)

    return col


def map_colormap_to_levels(levels, colors=None):
    if colors is not None:
        try:
            color_ref = Colors(field_name=colors).rgb_color_reference
            return color_ref
        except:
            pass

    return Colors(levels=levels).rgb_color_reference

def cn_legend(ax, frameon=True, loc=2, bbox_to_anchor=(0., 1.), title='Copy Number'):
    """ Display a legend for copy number state colors

    Parameters
    ----------
    ax : Axes
        matplotlib Axes on which to show legend
    frameon : bool, optional
        show frame, by default True
    loc : int, optional
        location of the legend, by default 2
    bbox_to_anchor : tuple, optional
        bounding box to which to anchor legend location, by default (0., 1.)

    Returns
    -------
    Legend
        legend object
    """

    color_reference = Colors('cn').hex_color_reference

    states = []
    patches = []
    for s, h in color_reference.items():
        states.append(s)
        patches.append(Patch(facecolor=h, edgecolor=h))

    ncol = min(3, int(len(states) ** (1 / 2)))

    legend = ax.legend(patches, states, ncol=ncol,
                       frameon=frameon, loc=loc, bbox_to_anchor=bbox_to_anchor,
                       facecolor='white', edgecolor='white', fontsize='4',
                       title=title, title_fontsize='6')
    legend.set_zorder(level=200)

    return legend
