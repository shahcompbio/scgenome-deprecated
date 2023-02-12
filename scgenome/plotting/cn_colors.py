import numpy as np

from matplotlib.patches import Patch
from numpy import ndarray


color_reference = {
    0:'#3182BD',
    1:'#9ECAE1',
    2:'#CCCCCC',
    3:'#FDCC8A',
    4:'#FC8D59',
    5:'#E34A33',
    6:'#B30000',
    7:'#980043',
    8:'#DD1C77',
    9:'#DF65B0',
    10:'#C994C7',
    11:'#D4B9DA',
}


def hex_to_rgb(h):
    if h is None:
        return np.array((0, 0, 0), dtype=int)
    h = h.lstrip('#')
    return np.array(tuple(np.uint8(int(h[i:i+2], 16)) for i in (0, 2 ,4)), dtype=int)


def map_cn_colors(X: ndarray) -> ndarray:
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
    X_colors = np.zeros(X.shape + (3,), dtype=int)
    X_colors[X < 0, :] = 0
    X_colors[X > max(color_reference.keys()), :] = max(color_reference.keys())
    for state, hex in color_reference.items():
        X_colors[X == state, :] = hex_to_rgb(hex)
    return X_colors


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
    states = []
    patches = []
    for s, h in color_reference.items():
        states.append(s)
        patches.append(Patch(facecolor=h, edgecolor=h))

    ncol = min(3, int(len(states)**(1/2)))

    legend = ax.legend(patches, states, ncol=ncol,
        frameon=frameon, loc=loc, bbox_to_anchor=bbox_to_anchor,
        facecolor='white', edgecolor='white', fontsize='4',
        title=title, title_fontsize='6')
    legend.set_zorder(level=200)

    return legend

