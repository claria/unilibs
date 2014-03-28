#! /usr/bin/env python2
import argparse
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from unilibs.pdf import PDF
from unilibs.baseplot import BasePlot
import helper

# PDF Plotting Style

pdf_fill_kwargs = [
    {
        'facecolor': '#19cce6',
        'alpha': 1.0,
        'edgecolor': 'Black',
        'linewidth': 1.0,
    },
    {
        'facecolor': 'None',
        'alpha': 1.0,
        'edgecolor': 'Blue',
        'linewidth': 1.0,
        'hatch': '//',
    },
    {
        'facecolor': 'none',
        'alpha': 1.0,
        'edgecolor': 'Green',
        'linewidth': 1.0,
        'hatch': '\\\\',
    },
    {
        'facecolor': 'none',
        'alpha': 1.0,
        'edgecolor': 'Gold',
        'linewidth': 1.0,
        'hatch': 'x',
    },


]
# Special HERAPDF Style
herapdf_fill_kwargs = {
    'exp': {'color': 'OrangeRed',
            'edgecolor': 'Black',
            'alpha': 1.0},
    'mod': {'color': 'Gold',
            'edgecolor': 'Black',
            'alpha': 1.0},
    'par': {'color': 'DarkGreen',
            'edgecolor': 'Black',
            'alpha': 1.0},
}


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('plot',
                        type=str,
                        choices=['pdf', 'dualpdf', 'ratio',
                                 'pdfratio', 'ratiooverview', 'pdfoverview'],
                        help='What to plot')

    parser.add_argument('pdfsets', type=str, nargs='+',
                        help='PDFSets')

    parser.add_argument('--flavors',
                        type=int,
                        nargs='+',
                        default=[0],
                        help='Flavors to plot')

    parser.add_argument("--q2", type=float,
                        default=1.9, help="Scale Q^2 of the PDF to evolve to")

    parser.add_argument("-o", "--output_folder",
                        default='./',
                        help="Output folder to save plots.")

    parser.add_argument("-t", "--output_type",
                        default=['pdf'],
                        type=str,
                        choices=['pdf', 'png', 'eps', 'svg', 'ps'],
                        nargs='+',
                        help="Plot Output format")

    parser.add_argument('--xscale', type=str, default='log',
                        choices=['linear', 'log'],
                        help='Xscale : log or linear')

    parser.add_argument('--yscale', type=str, default='linear',
                        choices=['linear', 'log'],
                        help='Yscale : log or linear')

    parser.add_argument('--uncertainty', type=str, default='default',
                        choices=['default', 'experimental', 'herapdf', 'none'],
                        help='Set wished uncertainty')

    parser.add_argument('--aroundunity', action='store_true',
                        help='Center uncertainties around unity')

    args = parser.parse_args()
    args = vars(args)

    lhgrid_filenames = args['pdfsets']

    # Create PDFs
    pdfs = []
    for lhgrid_filename in lhgrid_filenames:
        pdfs.append(PDF(lhgrid_filename,
                        flavors=args['flavors'],
                        q2=args['q2'],
                        x_range=np.logspace(-4., -0.0001, 501)))

    # if args['plot'] == 'pdf':
    #     for flavor in args['flavors']:
    #         kwargs = args
    #         output_fn = "{0}/{3}/{0}_{1}_{2}".format(pdfs[0].label,
    #                                                  flavor,
    #                                                  str(kwargs['q2']).replace('.', '_'),
    #                                                  kwargs['plot'])
    #         pdfplot = SimplePDFPlot(pdfs, flavor, kwargs['q2'],
    #                                 output_fn=output_fn,
    #                                 output_ext=kwargs['output_type'],
    #                                 )
    #         pdfplot.do_plot()
    #
    # elif args['plot'] == 'ratio':
    #     for flavor in args['flavors']:
    #         kwargs = args.copy()
    #         print 'ka;sdfj', kwargs['aroundunity']
    #         output_fn = "{0}/{3}/{0}_{1}_{2}".format(pdfs[0].label,
    #                                                  flavor,
    #                                                  str(kwargs['q2']).replace('.', '_'),
    #                                                  kwargs['plot'])
    #
    #         ratioplot = SimpleRatioPlot(pdfs, flavor, kwargs.pop('q2'),
    #                                     output_fn=output_fn,
    #                                     output_ext=kwargs.pop('output_type'),
    #                                     uncertainty=kwargs.pop('uncertainty'),
    #                                     **kwargs
    #                                     )
    #         ratioplot.do_plot()

    if args['plot'] == 'pdfratio':
        for flavor in args['flavors']:
            output_fn = os.path.join(args['output_folder'], 
                                 "{0}/{3}/{0}_{1}_{2}".format(pdfs[0].label,
                                                   flavor,
                                                   str(args['q2']).replace('.', '_'),
                                                   args['plot'])
                                 )
            pdfratioplot = SimplePDFRatioPlot(pdfs, flavor, args['q2'],
                                              output_fn=output_fn,
                                              output_ext=args['output_type'],
                                              xscale=args['xscale'],
                                              yscale=args['yscale'],
                                              uncertainty=args['uncertainty'])
            pdfratioplot.do_plot()

    elif args['plot'] == 'pdfoverview':
        output_fn = os.path.join(args['output_folder'],
                                 "{0}/{3}/{0}_{2}".format(pdfs[0].label,
                                                   'ov',
                                                   str(args['q2']).replace('.', '_'),
                                                   args['plot'])
                                 )
        pdfoverviewplot = SimplePDFOverviewPlot(
            pdfs, args['flavors'], args['q2'],
            output_fn=output_fn,
            output_ext=args['output_type'])
        pdfoverviewplot.do_plot()

    # elif args['plot'] == 'ratiooverview':
    #     output_fn = "{0}/{3}/{0}_{2}.{3}".format(pdfs[0].label,
    #                                              None,
    #                                              str(args['q2']).replace('.', '_'),
    #                                              args['output_type'],
    #                                              args['plot'])
    #     ratioplot = SimpleRatioOverviewPlot(pdfs, args['flavors'], args['q2'],
    #                                         output_fn=output_fn, )
    #     ratioplot.do_plot()
    #
    # elif args['plot'] == 'dualpdf':
    #     for flavor in args['flavors']:
    #         output_fn = "{0}/{3}/{0}_{1}_{2}.{3}".format(pdfs[0].label,
    #                                                      flavor,
    #                                                      str(args['q2']).replace('.', '_'),
    #                                                      args['output_type'],
    #                                                      args['plot'])
    #         dualpdfplot = SimpleDualPDFPlot(pdfs, flavor, args['q2'],
    #                                         output_fn=output_fn, )
    #         dualpdfplot.do_plot()


def plot_simple_pdf(ax,
                    pdfs,
                    flavor,
                    legend=False,
                    overview_scaling=False,
                    uncertainty='default',
                    **kwargs):
    if overview_scaling is True:
        plot_scale_factor = helper.get_plot_scalefactor(flavor)
    else:
        plot_scale_factor = 1.0

    for n, pdf in enumerate(pdfs):

        ax.plot(pdf.x, pdf.get_pdf_central(flavor) * plot_scale_factor,
                color=kwargs.get('color', pdf_fill_kwargs[n]['edgecolor']),
                linewidth=kwargs.get('linewidth', 0.5),
                linestyle=kwargs.get('linestyle', '--'))

        if uncertainty is 'none':
            if legend is True:
                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label=helper.get_pdflabel(
                                                     pdf.label),
                                                 fill=False,
                                                 linewidth=0,
                                                 edgecolor='none')
                ax.add_patch(p)

        # HERAPDF style plotting
        if uncertainty == 'herapdf':
            p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                             label=helper.get_pdflabel(
                                                 pdf.label),
                                             color='white', )
            ax.add_patch(p)
            # Add boxes with labels
            if legend is True:
                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label='Exp. Uncert.',
                                                 **herapdf_fill_kwargs[
                                                     'exp'])
                ax.add_patch(p)

                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label='Mod. Uncert.',
                                                 **herapdf_fill_kwargs[
                                                     'mod'])
                ax.add_patch(p)
                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label='Par. Uncert.',
                                                 **herapdf_fill_kwargs[
                                                     'par'])
                ax.add_patch(p)

            ax.fill_between(
                pdf.x,
                (pdf.get_pdf_central(flavor) -
                 np.sqrt(pdf.get_mod_uncert(flavor)[0] ** 2 +
                         pdf.get_pdf_uncert(flavor)[0] ** 2 +
                         pdf.get_par_uncert(flavor)[0] ** 2)) * plot_scale_factor,
                (pdf.get_pdf_central(flavor) +
                 np.sqrt(pdf.get_pdf_uncert(flavor)[1] ** 2 +
                         pdf.get_mod_uncert(flavor)[1] ** 2 +
                         pdf.get_par_uncert(flavor)[1] ** 2)) * plot_scale_factor,
                **herapdf_fill_kwargs['par'])

            ax.fill_between(
                pdf.x,
                (pdf.get_pdf_central(flavor) -
                 np.sqrt(pdf.get_mod_uncert(flavor)[0] ** 2 +
                         pdf.get_pdf_uncert(flavor)[0] ** 2)) * plot_scale_factor,
                (pdf.get_pdf_central(flavor) +
                 np.sqrt(pdf.get_pdf_uncert(flavor)[1] ** 2 +
                         pdf.get_mod_uncert(flavor)[1] ** 2)) * plot_scale_factor,
                **herapdf_fill_kwargs['mod'])

            ax.fill_between(
                pdf.x,
                (pdf.get_pdf_central(flavor) -
                 pdf.get_pdf_uncert(flavor)[0]) * plot_scale_factor,
                (pdf.get_pdf_central(flavor) +
                 pdf.get_pdf_uncert(flavor)[1]) * plot_scale_factor,
                **herapdf_fill_kwargs['exp'])

            # Common PDF Plotting
        if uncertainty == 'default' or uncertainty == 'experimental':
            if uncertainty == 'default':
                ax.fill_between(
                    pdf.x,
                    (pdf.get_pdf_central(flavor) -
                     pdf.get_tot_uncert(flavor)[0]) * plot_scale_factor,
                    (pdf.get_pdf_central(flavor) +
                     pdf.get_tot_uncert(flavor)[1]) * plot_scale_factor,
                    **pdf_fill_kwargs[n])
            if uncertainty == 'experimental':
                ax.fill_between(
                    pdf.x,
                    (pdf.get_pdf_central(flavor) -
                     pdf.get_pdf_uncert(flavor)[0]) * plot_scale_factor,
                    (pdf.get_pdf_central(flavor) +
                     pdf.get_pdf_uncert(flavor)[1]) * plot_scale_factor,
                    **pdf_fill_kwargs[n])

            if legend is True:
                fill = not 'hatch' in pdf_fill_kwargs[n]
                hatch = pdf_fill_kwargs[n].get('hatch', None)
                if fill:
                    color = pdf_fill_kwargs[n]['facecolor']
                else:
                    color = pdf_fill_kwargs[n]['edgecolor']

                alpha = pdf_fill_kwargs[n]['alpha']
                linewidth = pdf_fill_kwargs[n]['linewidth']
                edgecolor = pdf_fill_kwargs[n]['edgecolor']

                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label=helper.get_pdflabel(
                                                     pdf.label),
                                                 fill=fill,
                                                 facecolor=color,
                                                 alpha=alpha,
                                                 linewidth=linewidth,
                                                 edgecolor=edgecolor,
                                                 hatch=hatch)
                ax.add_patch(p)


def plot_simple_ratio(ax, pdfs, flavor, legend=False, trueratio=False,
                      uncertainty='default', y_shift=0.0):
    ref_pdf = pdfs[0]
    for n, pdf in enumerate(pdfs):

        #ax.plot(pdf.x, pdf.get_pdf_central(flavor) /
                       #ref_pdf.get_pdf_central(flavor) -1.,
                #color=pdf_fill_kwargs[n]['edgecolor'],
                #linewidth=1.5,
                #linestyle='--')

        # HERAPDF style plotting
        if uncertainty == 'herapdf':
            # Add boxes with labels
            if legend is True:
                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label='Exp. Uncert.',
                                                 **herapdf_fill_kwargs['exp'])
                ax.add_patch(p)
                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label='Mod. Uncert.',
                                                 **herapdf_fill_kwargs['mod'])
                ax.add_patch(p)
                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label='Par. Uncert.',
                                                 **herapdf_fill_kwargs['par'])
                ax.add_patch(p)

            ax.fill_between(
                pdf.x,
                (pdf.get_pdf_central(flavor) -
                 np.sqrt(pdf.get_mod_uncert(flavor)[0] ** 2 +
                         pdf.get_pdf_uncert(flavor)[0] ** 2 +
                         pdf.get_par_uncert(flavor)[0] ** 2)) /
                pdf.get_pdf_central(flavor) + y_shift,
                (pdf.get_pdf_central(flavor) +
                 np.sqrt(pdf.get_pdf_uncert(flavor)[1] ** 2 +
                         pdf.get_mod_uncert(flavor)[1] ** 2 +
                         pdf.get_par_uncert(flavor)[1] ** 2)) /
                pdf.get_pdf_central(flavor) + y_shift,
                **herapdf_fill_kwargs['par'])

            ax.fill_between(
                pdf.x,
                (pdf.get_pdf_central(flavor) -
                 np.sqrt(pdf.get_mod_uncert(flavor)[0] ** 2 +
                         pdf.get_pdf_uncert(flavor)[0] ** 2)) /
                pdf.get_pdf_central(flavor) + y_shift,
                (pdf.get_pdf_central(flavor) +
                 np.sqrt(pdf.get_pdf_uncert(flavor)[1] ** 2 +
                         pdf.get_mod_uncert(flavor)[1] ** 2)) /
                pdf.get_pdf_central(flavor) + y_shift,
                **herapdf_fill_kwargs['mod'])

            ax.fill_between(
                pdf.x,
                (pdf.get_pdf_central(flavor) -
                 pdf.get_pdf_uncert(flavor)[0]) /
                pdf.get_pdf_central(flavor) + y_shift,
                (pdf.get_pdf_central(flavor) +
                 pdf.get_pdf_uncert(flavor)[1]) /
                pdf.get_pdf_central(flavor) + y_shift,
                **herapdf_fill_kwargs['exp'])
            # Common PDF Plotting
        if uncertainty == 'default' or uncertainty == 'experimental':
            if uncertainty == 'default':
                if trueratio is True:
                    ax.fill_between(
                        pdf.x,
                        (pdf.get_pdf_central(flavor) -
                         pdf.get_tot_uncert(flavor)[0]) /
                        ref_pdf.get_pdf_central(flavor) + y_shift,
                        (pdf.get_pdf_central(flavor) +
                         pdf.get_tot_uncert(flavor)[1]) /
                        ref_pdf.get_pdf_central(flavor) + y_shift,
                        **pdf_fill_kwargs[n])
                else:
                    ax.fill_between(
                        pdf.x,
                        (pdf.get_pdf_central(flavor) -
                         pdf.get_tot_uncert(flavor)[0]) /
                        pdf.get_pdf_central(flavor) + y_shift,
                        (pdf.get_pdf_central(flavor) +
                         pdf.get_tot_uncert(flavor)[1]) /
                        pdf.get_pdf_central(flavor) + y_shift,
                        **pdf_fill_kwargs[n])
            if uncertainty == 'experimental':
                if trueratio is True:
                    ax.fill_between(
                        pdf.x,
                        (pdf.get_pdf_central(flavor) -
                         pdf.get_pdf_uncert(flavor)[0]) /
                        ref_pdf.get_pdf_central(flavor) + y_shift,
                        (pdf.get_pdf_central(flavor) +
                         pdf.get_pdf_uncert(flavor)[1]) /
                        ref_pdf.get_pdf_central(flavor) + y_shift,
                        **pdf_fill_kwargs[n])
                else:
                    ax.fill_between(
                        pdf.x,
                        (pdf.get_pdf_central(flavor) -
                         pdf.get_pdf_uncert(flavor)[0]) /
                        pdf.get_pdf_central(flavor) + y_shift,
                        (pdf.get_pdf_central(flavor) +
                         pdf.get_pdf_uncert(flavor)[1]) /
                        pdf.get_pdf_central(flavor) + y_shift,
                        **pdf_fill_kwargs[n])

            if legend is True:
                fill = not 'hatch' in pdf_fill_kwargs[n]
                hatch = pdf_fill_kwargs[n].get('hatch', None)
                if fill:
                    color = pdf_fill_kwargs[n]['facecolor']
                else:
                    color = pdf_fill_kwargs[n]['edgecolor']

                alpha = pdf_fill_kwargs[n]['alpha']
                linewidth = pdf_fill_kwargs[n]['linewidth']
                edgecolor = pdf_fill_kwargs[n]['edgecolor']

                p = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                                 label=helper.get_pdflabel(
                                                     pdf.label),
                                                 fill=fill,
                                                 facecolor=color,
                                                 alpha=alpha,
                                                 linewidth=linewidth,
                                                 edgecolor=edgecolor,
                                                 hatch=hatch)
                ax.add_patch(p)


class SimplePDFRatioPlot(BasePlot):
    def __init__(self, pdfs, flavor, q2, **kwargs):
        super(SimplePDFRatioPlot, self).__init__(output_fn=kwargs['output_fn'],
                                                 output_ext=kwargs['output_ext'],
                                                 style=kwargs.get('style', 'none'))
        self.pdfs = pdfs
        self.flavor = flavor
        self.q2 = q2
        self.props = kwargs
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1], hspace=0.05)
        self.ax1 = plt.subplot(gs[0])
        self.ax2 = plt.subplot(gs[1])

    def init_matplotlib(self):
        super(SimplePDFRatioPlot, self).init_matplotlib()
        matplotlib.rcParams['xtick.major.pad'] = 6
        matplotlib.rcParams['xtick.minor.pad'] = 8

    def prepare(self):

        self.set_style(style='cmsprel', show_cme=False, ax=self.ax1)
        self.set_preset_text(self.ax1,
                             r'{0}, $Q^2{{=}}{1}\mathrm{{GeV}}^2$'.format(
                                 helper.get_partonlabel(self.flavor),
                                 helper.get_q2label(self.q2)), )

    def produce(self):
        #Plot PDF in ax1

        plot_simple_pdf(self.ax1, self.pdfs, self.flavor, legend=True,
                        uncertainty=self.props.get('uncertainty', 'default'))
        #Plot Ratio in ax2
        plot_simple_ratio(self.ax2, self.pdfs, self.flavor,
                          uncertainty=self.props.get('uncertainty', 'default'),
                          #trueratio= not self.props.get('aroundunity', False),
                          trueratio = False,
                          y_shift=-1.0)
        self.ax1.set_xscale(self.props.get('xscale', 'log'))
        self.ax1.set_yscale(self.props.get('yscale', 'linear'), nonposy='clip')
        self.ax1.set_xlabel("$x$", ha='right')
        self.ax1.set_ylabel(r'$xf(x,Q^2)$', y = 1.0, size='large', ha='right')
        self.ax1.minorticks_on()
        self.ax1.set_xlim(1E-4, 0.95)
        self.ax1.set_ylim(ymin=0.)
        #if self.flavor == 8:
        #    self.ax1.set_ylim(ymax=0.6)

        self.ax1.legend(loc='best', prop={'size': 14})
        # minorLocator   = MultipleLocator(0.1)
        self.ax1.yaxis.grid(True, which='major')
        self.ax1.xaxis.grid(True, which='major')
        self.ax1.xaxis.set_ticklabels([])

        self.ax2.set_xscale(self.props.get('xscale', 'log'))
        self.ax2.set_xlabel("$x$", x=1.0, ha='right', size='x-large')
        self.ax2.set_ylabel(r'Fract. Uncert.')
        self.ax2.minorticks_on()
        self.ax2.autoscale(tight=True)
        self.ax2.set_xlim(1E-4, 0.95)
        self.ax2.set_ylim(-0.49, 0.49)

        self.ax2.legend(loc='best', prop={'size': 14})
        # minorLocator   = MultipleLocator(0.1)
        self.ax2.yaxis.grid(True, which='major')
        self.ax2.xaxis.grid(True, which='major')
        #self.autoscale(self.ax1, ymargin=0.1)

    def finalize(self):
        """
        Apply final settings, autoscale etc
        Save the plot
        """
        self._save_fig()
        plt.close(self.fig)


class SimplePDFPlot(BasePlot):
    def __init__(self, pdfs, flavor, q2, **kwargs):
        super(SimplePDFPlot, self).__init__(**kwargs)
        self.pdfs = pdfs
        self.q2 = q2
        self.flavor = flavor
        self.ax = self.fig.add_subplot(111)

    def prepare(self):
        self.set_style(style='cmsprel', show_cme=False, ax=self.ax)
        self.set_preset_text(self.ax,
                             r"{0}, $Q^2 = {1}\/\mathrm{{GeV}}^2$".format(
                                 helper.get_partonlabel(self.flavor),
                                 helper.get_q2label(self.q2)), )

    def produce(self):
        plot_simple_pdf(self.ax, self.pdfs, self.flavor, legend=True,
                        uncertainty=True)

        self.ax.set_xscale('log')
        self.ax.set_xlabel("$x$", ha='right')
        self.ax.set_ylabel(r'$xf(x,Q^2)$')
        self.ax.minorticks_on()
        self.ax.set_xlim(1E-3, 0.98)
        # self.ax.set_ylim(0.0, 1.)
        self.ax.autoscale(tight=True)

        self.ax.legend(loc='best', prop={'size': 14})
        # minorLocator   = MultipleLocator(0.1)
        self.ax.yaxis.grid(True, which='major')
        self.ax.xaxis.grid(True, which='major')

        # ax0.xaxis.set_ticklabels([])

    def finalize(self):
        """
        Apply final settings, autoscale etc
        Save the plot
        """
        self._save_fig()
        plt.close(self.fig)


# class SimpleRatioPlot(BasePlot):
#     def __init__(self, pdfs, flavor, q2, **kwargs):
#         super(SimpleRatioPlot, self).__init__(output_fn=kwargs['output_fn'],
#                                                  output_ext=kwargs['output_ext'],
#                                                  style=kwargs.get('style', 'none'))
#         self.pdfs = pdfs
#         self.q2 = q2
#         self.flavor = flavor
#         self.props = kwargs
#
#     def prepare(self):
#         self.ax = self.fig.add_subplot(111)
#         self.set_style(style='cmsprel', show_cme=False, ax=self.ax)
#         self.set_preset_text(self.ax,
#                              r"{0}, $Q^2={1}\/\mathrm{{GeV}}^2$".format(
#                                  helper.get_partonlabel(self.flavor),
#                                  helper.get_q2label(self.q2)), )
#
#     def produce(self):
#         plot_simple_ratio(self.ax, self.pdfs, self.flavor, legend=True,
#                           uncertainty=self.props.get('uncertainty', 'default'),
#                           trueratio= not self.props.get('aroundunity', False),
#                           )
#         print 'asdf', self.props['aroundunity']
#
#         self.ax.set_xscale('linear')
#         self.ax.set_xlabel("$x$")
#         self.ax.set_ylabel(r'Fractional uncertainties')
#         self.ax.minorticks_on()
#         self.ax.autoscale(tight=True)
#         self.ax.set_xlim(1E-4, 0.98)
#         self.ax.set_ylim(0.0, 2.0)
#
#         self.ax.legend(loc='upper left', prop={'size': 14})
#         # minorLocator   = MultipleLocator(0.1)
#         self.ax.yaxis.grid(True, which='major')
#         self.ax.xaxis.grid(True, which='major')
#
#     def finalize(self):
#         """
#         Apply final settings, autoscale etc
#         Save the plot
#         """
#         self._save_fig()
#         plt.close(self.fig)


class SimplePDFOverviewPlot(BasePlot):
    def __init__(self, pdfs, flavors, q2, **kwargs):
        super(SimplePDFOverviewPlot, self).__init__(**kwargs)
        self.pdfs = pdfs
        self.q2 = q2
        self.flavors = flavors

        self.p_label_pos = {0: {'x': 3E-3, 'y': 0.55, 'va': 'bottom', 'ha': 'left'},
                            7: {'x': 0.2, 'y': 0.34, 'va': 'bottom', 'ha': 'left'},
                            8: {'x': 0.07, 'y': 0.6, 'va': 'bottom', 'ha': 'left'},
                            9: {'x': 4E-4, 'y': 0.42, 'va': 'bottom', 'ha': 'left'},
        }

    def prepare(self):
        self.ax = self.fig.add_subplot(111)
        self.set_style(style='cmsprel', show_cme=False, ax=self.ax)
        self.set_preset_text(self.ax,
                             r" $Q^2={0}\/\mathrm{{GeV}}^2$".format(
                                 helper.get_q2label(self.q2)), )

        #Add text for the individual PDFs

    def produce(self):
        for flavor in self.flavors:
            plot_simple_pdf(self.ax, [self.pdfs[0], ], flavor,
                            legend=False,
                            overview_scaling=True,
                            uncertainty='default')
            plot_simple_pdf(self.ax, [self.pdfs[1], ],
                            flavor, legend=False,
                            overview_scaling=True,
                            linestyle='-',
                            linewidth=2,
                            uncertainty='none',
                            color='black')

            s = helper.get_partonlabel(flavor, short=True)
            self.ax.text(s=s, transform=self.ax.transData, color='Black',
                         **self.p_label_pos[flavor])

        self.ax.set_xscale('log')
        self.ax.set_xlabel("$x$", ha='right')
        self.ax.set_ylabel(r'$xf(x,Q^2)$')
        self.ax.minorticks_on()
        self.ax.autoscale(tight=True)
        self.ax.set_xlim(1E-4, 0.98)
        self.ax.set_ylim(ymin=0., ymax=1.0)
        #Get artists and labels for legend and chose which ones to display
        handles, labels = self.ax.get_legend_handles_labels()
        #Create custom artists
        pdf1 = matplotlib.patches.Rectangle((0, 0), 1, 1,
                                            **pdf_fill_kwargs[0])
        pdf2 = plt.Line2D((0, 0), (1, 1), color='black', linestyle='-')
        #Create legend from custom artist/label lists
        self.ax.legend(handles + [pdf1, pdf2],
                       labels + [helper.get_pdflabel(self.pdfs[0].label),
                                 helper.get_pdflabel(self.pdfs[1].label)],
                       loc='best', prop={'size': 14})

        # minorLocator   = MultipleLocator(0.1)
        self.ax.yaxis.grid(True, which='major')
        self.ax.xaxis.grid(True, which='major')

    def finalize(self):
        """
        Apply final settings, autoscale etc
        Save the plot
        """
        self._save_fig()
        plt.close(self.fig)


# class SimpleRatioOverviewPlot(BasePlot):
#     def __init__(self, pdfs, flavors, q2, **kwargs):
#         super(SimpleRatioOverviewPlot, self).__init__(**kwargs)
#         self.pdfs = pdfs
#         self.q2 = q2
#         self.flavors = flavors
#
#     def prepare(self):
#         self.ax = self.fig.add_subplot(111)
#         self.set_style(style='cmsprel', show_cme=False, ax=self.ax)
#         self.set_preset_text(self.ax,
#                              r" $Q^2={0}\/\mathrm{{GeV}}^2$".format(
#                                  helper.get_q2label(self.q2)), )
#
#     def produce(self):
#
#         for flavor in self.flavors:
#
#             ref_pdf = self.pdfs[0]
#             for n, pdf in enumerate(self.pdfs):
#                 if pdf is ref_pdf:
#                     continue
#                 self.ax.plot(pdf.x, pdf.get_pdf_central(flavor) /
#                                     ref_pdf.get_pdf_central(flavor),
#                              linewidth=0.5,
#                              linestyle='--',
#                              label=str(flavor))
#
#         self.ax.set_xscale('log')
#         self.ax.set_xlabel("$x$")
#         self.ax.set_ylabel(r'Fractional uncertainties')
#         self.ax.minorticks_on()
#         self.ax.autoscale(tight=True)
#         self.ax.set_xlim(3E-4, 0.98)
#         self.ax.set_ylim(0.5, 1.5)
#
#         self.ax.legend(loc='best', prop={'size': 14})
#         # minorLocator   = MultipleLocator(0.1)
#         self.ax.yaxis.grid(True, which='major')
#         self.ax.xaxis.grid(True, which='major')
#
#
#     def finalize(self):
#         """
#         Apply final settings, autoscale etc
#         Save the plot
#         """
#         self._save_fig()
#         plt.close(self.fig)


# class SimpleDualPDFPlot(BasePlot):
#     def __init__(self, pdfs, flavor, q2, **kwargs):
#         super(SimpleDualPDFPlot, self).__init__(**kwargs)
#         self.pdfs = pdfs
#         self.flavor = flavor
#         self.q2 = q2
#
#     def init_matplotlib(self):
#         super(SimpleDualPDFPlot, self).init_matplotlib()
#         matplotlib.rcParams['xtick.major.pad'] = 6
#         matplotlib.rcParams['xtick.minor.pad'] = 8
#
#     def prepare(self):
#         gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.05)
#         self.ax1 = plt.subplot(gs[0])
#         self.ax2 = plt.subplot(gs[1])
#
#         self.set_style(style='cmsprel', show_cme=False, ax=self.ax1)
#         self.set_preset_text(self.ax1,
#                              r'{0}, $Q^2={1}\mathrm{{GeV}}^2$'.format(
#                                  helper.get_partonlabel(self.flavor),
#                                  helper.get_q2label(self.q2)), )
#
#     def produce(self):
#         from matplotlib.ticker import MaxNLocator
#
#         plot_simple_ratio(self.ax1, [self.pdfs[0], ], self.flavor, legend=True)
#         plot_simple_ratio(self.ax2, [self.pdfs[
#                                          1], ], self.flavor, legend=False)
#
#         self.fig.text(0.0, 0.5, 'Fractional uncertaintie',
#                       ha='center', va='center',
#                       rotation='vertical')
#
#         self.ax1.set_xscale('log')
#         self.ax1.set_xlabel("$x$")
#         # self.ax1.set_ylabel(r'Fractional uncertainties')
#         self.ax1.minorticks_on()
#         self.ax1.autoscale(tight=True)
#         self.ax1.set_xlim(3E-4, 0.98)
#         self.ax1.set_ylim(0.801, 1.199)
#
#         self.ax1.legend(loc='upper left', prop={'size': 14})
#         # minorLocator   = MultipleLocator(0.1)
#         self.ax1.yaxis.set_major_locator(MaxNLocator(8, prune='both'))
#         self.ax1.yaxis.grid(True, which='major')
#         self.ax1.xaxis.grid(True, which='major')
#         self.ax1.xaxis.set_ticklabels([])
#
#         self.ax2.set_xscale('log')
#         self.ax2.set_xlabel("$x$")
#         # self.ax2.set_ylabel(r'Fractional uncertainties')
#         self.ax2.minorticks_on()
#         self.ax2.autoscale(tight=True)
#         self.ax2.set_xlim(3E-4, 0.98)
#         self.ax2.set_ylim(0.801, 1.199)
#
#         self.ax2.legend(loc='best', prop={'size': 14})
#         # minorLocator   = MultipleLocator(0.1)
#         self.ax2.yaxis.grid(True, which='major')
#         self.ax2.xaxis.grid(True, which='major')
#         self.ax2.yaxis.set_major_locator(MaxNLocator(8, prune='both'))
#
#     def finalize(self):
#         """
#         Apply final settings, autoscale etc
#         Save the plot
#         """
#         self._save_fig()
#         plt.close(self.fig)


if __name__ == '__main__':
    sys.exit(main())
