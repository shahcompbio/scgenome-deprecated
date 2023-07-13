from collections import defaultdict

import matplotlib
import numpy as np
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
from numpy import ndarray


class Colors(object):
    def __init__(self, category):

        category = category.lower()
        assert category in ['cn', 'giemsa', 'cyto_band_giemsa_stain']

        if category == 'cn':
            self.hex_color_reference = self.cn_color_reference
        elif category == 'giemsa' or category == 'cyto_band_giemsa_stain':
            self.hex_color_reference = self.cyto_band_giemsa_stain_hex_color_reference

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

    color_reference = Colors('cn').rgb_color_reference

    col = []

    for idx in range(len(X)):
        col.append(np.array([color_reference[v] for v in X[idx]]))

    col = np.array(col)

    return col


def infer_colormap_from_levels(num_levels):
    if num_levels <= 10:
        return 'tab10'
    elif num_levels <= 20:
        return 'tab20'
    else:
        return 'hsv'


def map_colormap_to_levels(levels, colors=None):
    if colors is None:
        cmap_name = infer_colormap_from_levels(len(levels))
        cmap = matplotlib.cm.get_cmap(cmap_name)
        level_colors = dict(zip(levels, list(cmap(np.linspace(0, 1, len(levels))))))
    else:
        level_colors = Colors(colors).rgb_color_reference

    return level_colors


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
