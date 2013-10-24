#! /usr/bin/env python2

import sys

import numpy

sys.path.append("/home/aem/uni/sw/lhapdf/lib/python2.7/site-packages")
#os.environ['LHAPATH'] = '/home/aem/uni/sw/lhapdf/share/lhapdf/PDFsets'
#os.environ['LHAPATH'] = '/home/aem/uni/ana/pdfs/'

import lhapdf


class PDF(object):

    def __init__(self,
                 lhgrid_filename,
                 flavors=range(-6, 7),
                 q2=1.9,
                 x_range=numpy.logspace(-4, -0.001, 201)):

        self._lhgrid_filename = lhgrid_filename
        self._flavors = flavors
        self._identify_pdfset()
        self._xrange = x_range
        self._q2 = q2
        self._q = numpy.sqrt(q2)

        self._pdf = {}
        self._varpdf = {}
        self._pdf_uncert = {}
        self._pdf_central = {}

        self._mod_uncert = {}
        self._par_uncert = {}

    def get_x(self):
        return self._xrange

    x = property(get_x)

    def get_label(self):
        return self._lhgrid_filename.split('.')[0]

    label = property(get_label)

    def has_var(self):
        return self._has_var

    has_var = property(has_var)

    def _identify_pdfset(self):

        if self._lhgrid_filename.startswith('CT10'):
            self._pdf_type = 'EV'
            self._has_var = False
        elif self._lhgrid_filename.startswith('MSTW'):
            self._pdf_type = 'EV'
            self._has_var = False
        elif self._lhgrid_filename.startswith('NNPDF'):
            self._pdf_type = 'MC'
            self._has_var = False
        elif self._lhgrid_filename.startswith('HERAMC'):
            self._pdf_type = 'MC'
            self._has_var = False
        elif self._lhgrid_filename.startswith('HERA'):
            self._pdf_type = 'EV'
            self._has_var = True
        else:
            raise Exception('No PDF type identified: %s', self._lhgrid_filename)

    def _read_lhapdf(self, lhgrid_filename):
        pdf = {}
        for flavor in self._flavors:
            pdf[flavor] = self._get_lhapdf_flavor(flavor, lhgrid_filename)
        return pdf

    def _get_lhapdf_flavor(self, flavor, lhgrid_filename):

        lhapdf.initPDFSetByName(lhgrid_filename)
        npdfs = lhapdf.numberPDF()
        pdf = numpy.zeros((npdfs + 1, self._xrange.size))
        for member in range(0, npdfs + 1):
            lhapdf.initPDF(member)
            for (i, xi) in enumerate(self._xrange):
                if flavor < 7:
                    pdf[member][i] = lhapdf.xfx(xi, self._q, flavor)
                elif flavor == 7:
                    #DVAL 1-(-1)
                    pdf[member][i] = lhapdf.xfx(xi, self._q, 1) - \
                                     lhapdf.xfx(xi, self._q, -1)
                elif flavor == 8:
                    #UVAL 2-(-2)
                    pdf[member][i] = lhapdf.xfx(xi, self._q, 2) - \
                                     lhapdf.xfx(xi, self._q, -2)
                elif flavor == 9:
                    #Light sea: xS=2(xubar + xdbar + xsbar)
                    pdf[member][i] = 2*(lhapdf.xfx(xi, self._q, -1) + \
                                     lhapdf.xfx(xi, self._q, -2) + \
                                     lhapdf.xfx(xi, self._q, -3))
                else:
                    raise Exception('Flavor not defined')
        return pdf

    def get_pdf_central(self, flavor):
        if not flavor in self._pdf_central:
            self._calc_pdf_central(flavor)
        return self._pdf_central[flavor]

    def _calc_pdf_central(self, flavor):
        if self._pdf_type == 'MC':
            self._calc_pdf_mean(flavor)
        elif self._pdf_type == 'EV':
            self._calc_pdf_zeroth(flavor)

    def _calc_pdf_mean(self, flavor):
        if not flavor in  self._pdf:
            self._pdf = self._read_lhapdf(self._lhgrid_filename)
        self._pdf_central[flavor] = numpy.mean(
            self._pdf[flavor][1:], axis=0)

    def _calc_pdf_zeroth(self, flavor):
        if not flavor in self._pdf:
            self._pdf = self._read_lhapdf(self._lhgrid_filename)
        self._pdf_central[flavor] = self._pdf[flavor][0]

    def get_pdf_uncert(self, flavor):
        if not flavor in self._pdf_uncert:
            self._calc_pdf_uncert(flavor)
        return self._pdf_uncert[flavor]

    def get_mod_uncert(self, flavor):
        if not self._has_var:
            return None
        if not flavor in self._mod_uncert:
            self._calc_mod_uncert(flavor)
        return self._mod_uncert[flavor]

    def get_par_uncert(self, flavor):
        if not self._has_var:
            return None
        if not flavor in self._par_uncert:
            self._calc_par_uncert(flavor)
        return self._par_uncert[flavor]

    def get_tot_uncert(self, flavor):
        if not self._has_var:
            return self.get_pdf_uncert(flavor)
        else:
            tot_uncert_down = numpy.sqrt(self.get_pdf_uncert(flavor)[0] ** 2 +
                                         self.get_mod_uncert(flavor)[0] ** 2 +
                                         self.get_par_uncert(flavor)[0] ** 2)
            tot_uncert_up = numpy.sqrt(self.get_pdf_uncert(flavor)[1] ** 2 +
                                       self.get_mod_uncert(flavor)[1] ** 2 +
                                       self.get_par_uncert(flavor)[1] ** 2)
            return tot_uncert_down, tot_uncert_up

    def _calc_pdf_uncert(self, flavor):

        if not self._pdf:
            self._pdf = self._read_lhapdf(self._lhgrid_filename)

        if self._pdf_type == 'MC':
            self._calc_pdf_mc_uncert(flavor)
        elif self._pdf_type == 'EV':
            self._calc_pdf_ev_uncert(flavor)

    def _calc_mod_uncert(self, flavor):
        if not self._varpdf:
            self._varpdf = self._read_lhapdf(
                self._lhgrid_filename.replace('EIG', 'VAR'))

        mod_uncert_down = numpy.zeros(self._varpdf[flavor].shape[1])
        mod_uncert_up = numpy.zeros(self._varpdf[flavor].shape[1])

        #print "Assuming member 1 to 8 are model uncertainties"

        for member in range(1, 9):
            mod_uncert_up += numpy.square(numpy.maximum(
                self._varpdf[flavor][member] -
                self._varpdf[flavor][0], 0.0))
            mod_uncert_down += numpy.square(numpy.minimum(
                self._varpdf[flavor][member] -
                self._varpdf[flavor][0], 0.0))

        self._mod_uncert[flavor] = (numpy.sqrt(mod_uncert_down),
                                    numpy.sqrt(mod_uncert_up))

    def _calc_par_uncert(self, flavor):
        if not self._varpdf:
            self._varpdf = self._read_lhapdf(
                self._lhgrid_filename.replace('EIG', 'VAR'))

        par_uncert_down = numpy.zeros(self._varpdf[flavor].shape[1])
        par_uncert_up = numpy.zeros(self._varpdf[flavor].shape[1])

        #print "Assuming member 9 to x are par variations"

        for member in range(9, self._varpdf[flavor].shape[0]):
        #if member in [11,12]:
        #continue
        #for member in range(9,):
            par_uncert_up = numpy.maximum(par_uncert_up,
                                          self._varpdf[flavor][member] -
                                          self._varpdf[flavor][0])
            par_uncert_down = numpy.maximum(par_uncert_down,
                                            self._varpdf[flavor][0] -
                                            self._varpdf[flavor][member])

        self._par_uncert[flavor] = (par_uncert_down, par_uncert_up)

    def _calc_pdf_mc_uncert(self, flavor):
        if not self._pdf:
            self._pdf = self._read_lhapdf(self._lhgrid_filename)
        #npdfs, nbins = pdf.shape
        std = numpy.std(self._pdf[flavor][1:], axis=0)
        self._pdf_uncert[flavor] = (std, std)

    def _calc_pdf_ev_uncert(self, flavor):

        #npdfs,nbins = pdf.shape
        pdf_uncert_up = numpy.zeros(self._pdf[flavor].shape[1])
        pdf_uncert_down = numpy.zeros(self._pdf[flavor].shape[1])
        for i in range(1, (self._pdf[flavor].shape[0]) / 2 + 1):
            pdf_uncert_up += numpy.square(numpy.maximum(
                numpy.maximum(
                    self._pdf[flavor][2 * i - 1] - self._pdf[flavor][0],
                    self._pdf[flavor][2 * i] - self._pdf[flavor][0]),
                numpy.zeros((self._pdf[flavor].shape[1],))))
            pdf_uncert_down += numpy.square(numpy.minimum(
                numpy.minimum(
                    self._pdf[flavor][2 * i - 1] - self._pdf[flavor][0],
                    self._pdf[flavor][2 * i] - self._pdf[flavor][0]),
                numpy.zeros(self._pdf[flavor].shape[1])))

        pdf_uncert_up = numpy.sqrt(pdf_uncert_up)
        pdf_uncert_down = numpy.sqrt(pdf_uncert_down)

        self._pdf_uncert[flavor] = (pdf_uncert_down, pdf_uncert_up)


    #def get_symmetric_error(pdf):
    #"""Calculate Percentage Errors (experimental uncertainty) out of eigenvector
    #up/down variations"""
    #npdfs,nbins = pdf.shape
    #pdf_uncert = numpy.zeros(nbins)

    #for i in range(1,len(pdf)/2 + 1):
    #pdf_uncert += numpy.square(pdf[2*i-1] - pdf[2*i])

    #pdf_uncert = 0.5 * numpy.sqrt(pdf_uncert)

    #return (pdf_uncert, pdf_uncert)

    #if __name__ == '__main__':
    #sys.exit(main(sys.argv[1:]))

    #def calc_var_error(pdf):
    #npdfs,nbins =pdf.shape
    #uncert_up = numpy.zeros(nbins)
    #uncert_down = numpy.zeros(nbins)
    #
    #print "Assuming member 1 to 8 are model uncertainties"
    #for i in range(1,8):
    #uncert_up += numpy.square(numpy.maximum(pdf[i]-pdf[0],0.0))
    #uncert_down += numpy.square(numpy.minimum(pdf[i]-pdf[0],0.0))
    #
    #return numpy.sqrt(uncert_down), numpy.sqrt(uncert_up)
    #
    #def calc_par_error(pdf):
    #npdfs,nbins = pdf.shape
    #uncert_up = numpy.zeros(nbins)
    #uncert_down = numpy.zeros(nbins)
    #print "Assuming member 9 to x ar par variations"
    #
    #for i in range(9,npdfs):
    #uncert_up = numpy.maximum(uncert_up,pdf[i]-pdf[0])
    #uncert_down = numpy.maximum(uncert_down,pdf[0]-pdf[i])
    #
    #return uncert_down, uncert_up


