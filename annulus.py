
from matplotlib.transforms import Affine2D

import mpl_toolkits.axisartist.floating_axes as floating_axes

import numpy as np
import  mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator, MaxNLocator, \
     DictFormatter

import matplotlib.pyplot as plt

def setup_axes2(fig, rect):
    """
    With custom locator and formatter.
    Note that the extreme values are swapped.
    """

    #tr_scale = Affine2D().scale(np.pi/180., 1.)

    tr = PolarAxes.PolarTransform()

    pi = np.pi
    angle_ticks = [(0, r"0"),
                   (.25*pi, r"45"),
                   (.5*pi, r"90"),(.75*pi, r"135")]
    grid_locator1 = FixedLocator([v for v, s in angle_ticks])
    tick_formatter1 = DictFormatter(dict(angle_ticks))

    grid_locator2 = MaxNLocator(20)

    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                        extremes=(2*pi, 0, 1.1, 0.90),
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1,
                                        tick_formatter2=None,
                                        )

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # create a parasite axes whose transData in RA, cz
    aux_ax = ax1.get_aux_axes(tr)

    aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                        # drawn twice, and possibly over some other
                        # artists. So, we decrease the zorder a bit to
                        # prevent this.

    return ax1, aux_ax


fig=plt.figure()

ax2, aux_ax2 = setup_axes2(fig, 111)

theta = np.random.rand(10)*.5*np.pi
radius = np.random.rand(10)+1.
aux_ax2.scatter(theta, radius)
plt.grid(True)
plt.show()
