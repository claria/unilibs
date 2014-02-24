import os
from abc import ABCMeta, abstractmethod
import numpy as np

import matplotlib
import matplotlib.pyplot as plt


class BasePlot(object):
    __metaclass__ = ABCMeta

    def __init__(self, output_fn='test', output_ext=('png',), style='none'):

        self.init_matplotlib()
        self.fig = plt.figure()

        self.output_fn = output_fn
        self.output_ext = output_ext
        self.style = style

    def do_plot(self):
        """
        Run all three plotting steps
        """
        self.prepare()
        self.produce()
        self.finalize()

    @abstractmethod
    def prepare(self, **kwargs):
        """
        Before plotting:
        Add axes to Figure, etc
        """
        pass

    @abstractmethod
    def produce(self):
        """
        Do the Plotting
        """
        pass

    @abstractmethod
    def finalize(self):
        """
        Apply final settings, autoscale etc
        Save the plot
        """
        self._save_fig()
        plt.close(self.fig)

    def _save_fig(self):
        """
        Save Fig to File and create directory structure
        if not yet existing.
        """
        #Check if directory exists and create if not
        directory = os.path.dirname(self.output_fn)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        for ext in self.output_ext:
            filename = "{}.{}".format(self.output_fn, ext)
            self.fig.savefig(filename, bbox_inches='tight')

    @staticmethod
    def init_matplotlib():
        """
        Initialize matplotlib with following rc
        """
        matplotlib.rcParams['lines.linewidth'] = 2
        matplotlib.rcParams['font.family'] = 'sans-serif'
        matplotlib.rcParams['font.style'] = 'normal'
        matplotlib.rcParams['font.size'] = 20.
        matplotlib.rcParams['legend.fontsize'] = 14.
        matplotlib.rcParams['text.usetex'] = False
        # Axes
        matplotlib.rcParams['axes.linewidth'] = 2.0
        # Saving
        matplotlib.rcParams['savefig.bbox'] = 'tight'
        matplotlib.rcParams['savefig.dpi'] = 90
        matplotlib.rcParams['savefig.format'] = 'png'

    #
    # Helper functions
    #
    @staticmethod
    def set_preset_text(ax, text, loc='topright', **kwargs):
        """
        Possible Positions : topleft, topright
        """
        if loc == 'topleft':
            kwargs.update({'x': 0.0, 'y': 1.01, 'va': 'bottom',
                           'ha': 'left'})
        elif loc == 'topright':
            kwargs.update({'x': 1.0, 'y': 1.01, 'va': 'bottom',
                           'ha': 'right'})
        else:
            raise Exception()

        ax.text(s=text, transform=ax.transAxes, color='Black', **kwargs)

    def set_style(self, ax, style, show_cme=True):
        """
        Some preset styles
        """
        if style == 'none':
            pass
        elif style == 'cmsprel':
            self.set_preset_text(ax, "CMS Preliminary", loc='topleft')
            if show_cme:
                self.set_preset_text(ax, r"$\sqrt{s} = 7\/ \mathrm{TeV}$",
                                     loc='topleft', )
        else:
            self.set_preset_text(ax, "CMS", loc='topleft')
            if show_cme:
                self.set_preset_text(ax, r"$\sqrt{s} = 7\/ \mathrm{TeV}$",
                                     loc='topleft', )

    @staticmethod
    def autoscale(ax, xmargin=0.0, ymargin=0.0, margin=0.0):
        # User defined autoscale with margins
        x0, x1 = tuple(ax.dataLim.intervalx)
        if margin > 0:
            xmargin = margin
            ymargin = margin
        if xmargin > 0:
            if ax.get_xscale() == 'linear':
                delta = (x1 - x0) * xmargin
                x0 -= delta
                x1 += delta
            else:
                delta = (x1 / x0) ** xmargin
                x0 /= delta
                x1 *= delta
            ax.set_xlim(x0, x1)
        y0, y1 = tuple(ax.dataLim.intervaly)
        if ymargin > 0:
            if ax.get_yscale() == 'linear':
                delta = (y1 - y0) * ymargin
                y0 -= delta
                y1 += delta
            else:
                delta = (y1 / y0) ** ymargin
                y0 /= delta
                y1 *= delta
            ax.set_ylim(y0, y1)

    @staticmethod
    def log_locator_filter(x, pos):
        """
        Add minor tick labels in log plots at 2* and 5*
        """
        s = str(int(x))
        if len(s) == 4:
            return ''
        if s[0] in ('2', '5'):
            return s
        return ''

    @staticmethod
    def steppify_bin(arr, isx=False):
        """
        Produce stepped array of arr, also of x
        """
        if isx:
            newarr = np.array(zip(arr[0], arr[1])).ravel()
        else:
            newarr = np.array(zip(arr, arr)).ravel()
        return newarr

    @staticmethod
    def set(obj, *args, **kwargs):
        """
        Apply Settings in kwargs, while defaults are set
        """
        funcvals = []
        for i in range(0, len(args) - 1, 2):
            funcvals.append((args[i], args[i + 1]))
        funcvals.extend(kwargs.items())
        for s, val in funcvals:
            attr = getattr(obj, s)
            if callable(attr):
                attr(val)
            else:
                setattr(obj, attr, val)


class GenericPlot(BasePlot):
    """
    Very simple generic plotting script
    A list of datasets has to be provided.
    A dataset is a dict with x and y keys and dx,dy
    """
    def __init__(self, datasets,
                 output_fn='test.png',
                 output_ext=('png', ),
                 props=None,
                 **kwargs):
        super(GenericPlot, self).__init__(output_fn=output_fn,
                                          output_ext=output_ext,
                                          **kwargs)

        self.output_fn = output_fn
        self.datasets = datasets

        self.ax = self.fig.add_subplot(111)
        self.props = props if props else {}
        self.props.update(kwargs.get('post_props', {}))
        self.pre_props = kwargs.get('pre_props', {})

    def prepare(self, **kwargs):

        for artist, props in self.props.items():
            obj = getattr(self, artist)
            self.set(obj, **props)

    def produce(self):
        for dataset in self.datasets:
            plot_type = dataset.get('plot_type', 'plot')
            if plot_type == 'plot':
                self.ax.plot(dataset['x'],
                             dataset['y'],
                             label=dataset.get('label', ''),
                             **dataset.get('props', {}))
            elif plot_type == 'errorbar':
                self.ax.errorbar(x=dataset['x'],
                                 xerr=dataset['dx'],
                                 y=dataset['y'],
                                 yerr=dataset['dy'],
                                 fmt='+',
                                 label=dataset.get('label', ''),
                                 **dataset.get('props', {}))
            elif plot_type == 'fill_between':
                self.ax.fill_between(x=dataset['x'],
                                     y1=dataset['y'] - dataset['dy'][0],
                                     y2=dataset['y'] + dataset['dy'],
                                     #label=dataset.get('label',''),
                                     **dataset.get('props', {}))
        self.ax.legend()

    def finalize(self):

        for artist, props in self.props.items():
            obj = getattr(self, artist)
            self.set(obj, **props)

        self.autoscale(self.ax, margin=0.1)
        self._save_fig()
        plt.close(self.fig)
