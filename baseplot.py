import os

import matplotlib
import numpy

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


class SimplePlot(object):

    def __init__(self):
        self.prepare_matplotlib()
        self.fig = plt.figure()

    def prepare(self, **kwargs):
        """
        Before plotting:
        Add axes to Figure, etc
        """
        self.ax = self.fig.add_subplot(1, 1, 1)

    def produce(self):
        """
        Do the Plotting
        """
        pass

    def finalize(self, filepath='test.pdf'):
        """
        Apply final settings, autoscale etc
        Save the plot
        :param filepath:
        """

        self.autoscale(margin=0.1)

        directory = os.path.dirname(filepath)
        if not os.path.exists(directory):
            os.makedirs(directory)

        self.fig.savefig(filepath)

    def prepare_matplotlib(self):
        # matplotlib.use('agg')
        # matplotlib.use('pdf')

        matplotlib.rcParams['lines.linewidth'] = 2
        matplotlib.rcParams['font.family'] = 'sans-serif'
        matplotlib.rcParams['font.style'] = 'normal'
        matplotlib.rcParams['font.size'] = 22.
        matplotlib.rcParams['legend.fontsize'] = 14.
        matplotlib.rcParams['text.usetex'] = False
        # Axes
        matplotlib.rcParams['axes.linewidth'] = 2.0
        # Saving
        matplotlib.rcParams['savefig.bbox'] = 'tight'
        matplotlib.rcParams['savefig.dpi'] = 300
        matplotlib.rcParams['savefig.format'] = 'pdf'

    #
    # Helper functions
    #
    def set_preset_text(self, text, loc='topright', **kwargs):
        """
        Possible Positions : topleft, topright
        """

        if loc == 'topleft':
            kwargs.update({'x' : 0.0, 'y' : 1.01, 'va' : 'bottom',
                            'ha' : 'left'})
        elif loc == 'topright':
            kwargs.update({'x' : 1.0, 'y' : 1.01, 'va' : 'bottom',
                            'ha' : 'right'})
        else:
            raise Expection()
        print kwargs
        self.ax.text(s=text, transform=self.ax.transAxes, color='Black', **kwargs)

    def set_style(self, style, show_cme=True, **kwargs):
        """
        Some preset styles
        """
        if style == 'cmsprel':
            cmstext = "CMS Preliminary"
        else:
            cmstext = "CMS"

        self.set_preset_text(cmstext, loc='topleft')
        if show_cme:
            self.set_preset_text(r"$\sqrt{s} = 7\/ \mathrm{TeV}$", loc='topleft')

    def autoscale(self, xmargin=0.0, ymargin=0.0, margin=0.0):
        # User defined autoscale with margins
        x0, x1 = tuple(self.ax.dataLim.intervalx)
        if margin > 0:
            xmargin = margin
            ymargin = margin
        if xmargin > 0:
            if self.ax.get_xscale() == 'linear':
                delta = (x1 - x0) * xmargin
                x0 -= delta
                x1 += delta
            else:
                delta = (x1 / x0) ** xmargin
                x0 /= delta
                x1 *= delta
            self.ax.set_xlim(x0, x1)
        y0, y1 = tuple(self.ax.dataLim.intervaly)
        if ymargin > 0:
            if self.ax.get_yscale() == 'linear':
                delta = (y1 - y0) * ymargin
                y0 -= delta
                y1 += delta
            else:
                delta = (y1 / y0) ** ymargin
                y0 /= delta
                y1 *= delta
            self.ax.set_ylim(y0, y1)

    def log_locator_filter(self, x, pos):
        """Add minor tick labels in log plots at 2* and 5*
        """
        s = str(int(x))
        if len(s) == 4:
            return ''
        if s[0] in ('2', '5'):
            return s
        return ''

    def steppify_bin(self, arr, isx=False):
        """Produce stepped array of arr, also of x
        """
        if isx:
            newarr = numpy.array(zip(arr[0], arr[1])).ravel()
        else:
            newarr = numpy.array(zip(arr, arr)).ravel()
        return newarr
