"""
Provides plotting styles.

Anthony Brown Aug 2015 - Oct 2018
"""

import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler

from agabpylib.plotting.distinct_colours import get_distinct

def useagab(usetex=True, fontfam='serif', sroncolours=True, ncolors=4, axislinewidths=1, linewidths=2,
        lenticks=6):
    """
    Configure the plotting style to my liking.

    Parameters
    ----------

    None

    Keywords
    --------

    usetex : boolean
        Whether or not to use LaTeX text (default True).
    fontfam : boolean
        Font family to use (default 'serif')
    sroncolours : boolean
        If true use colour-blind proof distinct colours (https://personal.sron.nl/~pault/).
    ncolors : int
        Number of distinct colours to use (applies to SRON colours only, default 4)
    axislinewidths : float
        Width of lines used to draw axes (default 1)
    linewidths : float
        Width of lines used to draw plot elements (default 2)
    lenticks : float
        Length of major tickmarks in points (default 6, minor tick marks adjusted automatically)

    Returns
    -------

    Nothing
    """
    line_colours = get_distinct(ncolors)
    if (usetex):
        rc('text', usetex=True)
        rc('text.latex', preamble=r'\usepackage{amsmath}')
    rc('font', family=fontfam, size=18)
    rc('xtick.major', size=lenticks)
    rc('xtick.minor', size=lenticks*2/3)
    rc('ytick.major', size=lenticks)
    rc('ytick.minor', size=lenticks*2/3)
    rc('lines', linewidth=linewidths)
    rc('axes', linewidth=axislinewidths)
    if (sroncolours):
        rc('axes', prop_cycle=(cycler('color',line_colours)))
    else:
        rc('axes', prop_cycle=(cycler('color',plt.cm.tab10.colors)))
    rc('xtick', direction='out')
    rc('ytick', direction='out')
    rc('grid', color='cbcbcb')
    rc('grid', linestyle='-')
    rc('grid', linewidth=0.5)
    rc('grid', alpha=1.0)
    rc('figure', facecolor='ffffff')
    rc('figure', dpi=80)
    rc('figure.subplot', bottom=0.125)

def apply_tufte(ax, withgrid=False, minorticks=False, gridboth=False):
    """
    Apply the "Tufte" style to the plot axes contained in the input axis object. This mimics the sparse
    style advocated by Tufte in his book "The Visual Display of Quantitative Information".

    Parameters
    ----------

    ax - The axis object to configure.

    Keywords
    --------

    withgrid : boolean
        When true a grid is displayed in the plot background
    minorticks : boolean
        When true minor tickmarks are drawn.
    gridboth : boolean
        When true minor tickmarks are also used for the grid

    Returns
    -------

    Nothing.
    """

    # Move left and bottom spines outward by 5 points
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(ax.spines[axis].get_linewidth())
    ax.tick_params('both', width=ax.spines['bottom'].get_linewidth(), which='both')
    if withgrid:
        if (gridboth):
            ax.grid(which='both')
        else:
            ax.grid()
    if minorticks:
        ax.minorticks_on()
    ax.set_facecolor('w')
