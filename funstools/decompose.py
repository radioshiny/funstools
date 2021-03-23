from .cube2map import Cube2map
from .smooth import smooth3d, radsmo2d, boxsmo1d
from .makemaps import get_rms
import numpy as np
from scipy.signal import find_peaks
from astropy.modeling.models import Gaussian1D
from astropy.table import Table
from astropy.io.ascii import write
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


def predict_peaks(data, rms):
    if not isinstance(data, np.ndarray):
        raise TypeError("'data' should be given as numpy.ndarray.")
    if not isinstance(rms, float):
        raise TypeError("'rms' should be given as float.")
    if not len(data.shape) == 1:
        raise TypeError("'data' should be given as 1D-spectrum.")
    grad = -np.gradient(np.gradient(data))
    peaks, _ = find_peaks(grad, prominence=rms/3)
    return peaks[data[peaks] > rms]


def find_comps(data, rms, mask=None, ss=0., vs=0.):
    if not isinstance(data, np.ndarray):
        raise TypeError("'data' should be given as numpy.ndarray.")
    if not len(data.shape) == 3:
        raise TypeError("'data' should be given as 3D-cube.")
    if isinstance(rms, np.ndarray):
        if not data[0].shape == rms.shape:
            raise TypeError("'rms' shape is not match to data.")
    elif isinstance(rms, float):
        rms = np.full_like(data[0], rms)
    else:
        raise TypeError("'rms' should be given as float or 2D-image.")
    if mask is None:
        mask = np.full_like(data[0], True)
    ny, nx = data[0].shape
    comps = np.zeros_like(data)
    for y in range(ny):
        for x in range(nx):
            if not mask[y, x]:
                comps[:, y, x] = np.nan
                continue
            pl = predict_peaks(data[:, y, x], rms[y, x])
            comps[:, y, x][pl] = 1.
    ss = int(np.ceil(float(ss)/2))
    vs = int(np.ceil(vs))
    if ss == 0 and vs == 0:
        return comps, np.sum(comps, 0)
    elif vs == 0:
        comps_smo = boxsmo1d(radsmo2d(comps, ss), vs)
    else:
        comps_smo = boxsmo1d(radsmo2d(comps, ss), vs)*vs
    # elif vs == 0:
    #     comps_smo = smooth1d(radsmo2d(comps, ss), vs)
    # else:
    #     comps_smo = smooth1d(radsmo2d(comps, ss), vs)*vs
    comps_cont = np.zeros_like(data).astype('int')
    for y in range(ny):
        for x in range(nx):
            if not mask[y, x]:
                continue
            pl, _ = find_peaks(comps_smo[:, y, x], height=0.3)
            comps_cont[:, y, x][pl] = 1
    return comps_cont, np.sum(comps_cont, 0)


class Decompose(Cube2map):
    """
    The Decompose is a class intended to decompose multiple Gaussian components
    from molecular line profiles in the FITS cube data.
    """

    def __init__(self, cube=None, ext=None, getrms='both', rmssize=None, max_rms=None,
                 snr=3., spasmo_find=3, velsmo_find=3, spasmo_fit=1, velsmo_fit=1,
                 mlim=0.1, wmin=0.1, wmax=1.):
        """
        Construct a 'Decompose' object.

        Parameters:
            cube, ext, getrms, rmssize, max_rms, snr : see help(Cube2map)
            velsmo_find : float
                Velocity smoothing channel size for finding components.
            velsmo_fit: : float
                Velocity smoothing channel size for data used decomposing.
            spasmo_find : float
                Spatial smoothing pixel size for finding components.
            spasmo_fit: : float
                Spatial smoothing pixel size for data used decomposing.
            mlim : float
                Mean velocity bound limit for Gaussian fit.
                The bounds is from Vpeak-mlim to Vpeak+mlim.
            wmin : float
                Lower limit of FWHM (km/s) for Gaussian fit.
            wmax : float
                Upper limit of FWHM (km/s) for Gaussian fit.
        """
        Cube2map.__init__(self, cube, ext, getrms, rmssize, max_rms, snr, velsmo_fit, spasmo_fit, True, True)
        self._velsmo_find = float(velsmo_find)
        self._spasmo_find = float(spasmo_find)
        self._mlim = float(mlim)
        self._smin = float(wmin)/np.sqrt(8*np.log(2))
        self._smax = float(wmax)/np.sqrt(8*np.log(2))

    _datafc = None
    _mdatafc = None
    _rmsfc = None
    _comps = None
    _nc = None
    _det_comp = None
    _fitres1 = None
    _fitres2 = None
    _fitres3 = None
    _fitres_after_fof = None

    @property
    def data_for_finding(self):
        """
        Smoothed data for finding component.
        """
        if self._datafc is None:
            self._datafc = smooth3d(self.data, self._spasmo_find, self._velsmo_find)
        return self._datafc

    @property
    def mdata_for_finding(self):
        """
        Smoothed data with masking noise channels
        """
        if self._mdatafc is None:
            self._mdatafc = smooth3d(self.mdata, self._spasmo_find, self._velsmo_find)
        return self._mdatafc

    @property
    def rms_for_finding(self):
        """
        RMS noise level data for smoothed data.
        """
        if self._rmsfc is None:
            self._rmsfc = get_rms(self.data_for_finding, self._getrms, self._rmssize)
        return self._rmsfc

    @property
    def comps(self):
        """
        Result of finding component.
        An 3d-cube array marked the peak position to 1 and the rest to 0.
        """
        if self._comps is None:
            self._comps, self._nc = find_comps(self.mdata_for_finding, self.rms_for_finding,
                                               self.det, self._spasmo_find, self._velsmo_find)
        return self._comps

    @property
    def nc(self):
        """
        Result of finding component.
        An 2d-image array with number of components for each pixel.
        """
        if self._nc is None:
            self._comps, self._nc = find_comps(self.mdata_for_finding, self.rms_for_finding,
                                              self.det, self._spasmo_find, self._velsmo_find)
        return self._nc

    @property
    def det_comp(self):
        if self._det_comp is None:
            dc = self.nc > 0
            self._det_comp = self.det & dc
        return self._det_comp

    def _px(self, d, r):
        return np.where(self.comps[:, d, r] == 1)[0]

    @staticmethod
    def _gauss(cn):
        if cn == 1:
            def _gauss1(x, a1, m1, s1):
                return a1*np.exp(-(x-m1)**2/2/s1**2)
            return _gauss1
        elif cn == 2:
            def _gauss2(x, a1, m1, s1, a2, m2, s2):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2)
                return f
            return _gauss2
        elif cn == 3:
            def _gauss3(x, a1, m1, s1, a2, m2, s2, a3, m3, s3):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2)
                return f
            return _gauss3
        elif cn == 4:
            def _gauss4(x, a1, m1, s1, a2, m2, s2, a3, m3, s3, a4, m4, s4):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2) \
                    + a4*np.exp(-(x-m4)**2/2/s4**2)
                return f
            return _gauss4
        elif cn == 5:
            def _gauss5(x, a1, m1, s1, a2, m2, s2, a3, m3, s3, a4, m4, s4, a5, m5, s5):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2) \
                    + a4*np.exp(-(x-m4)**2/2/s4**2) + a5*np.exp(-(x-m5)**2/2/s5**2)
                return f
            return _gauss5
        elif cn == 6:
            def _gauss6(x, a1, m1, s1, a2, m2, s2, a3, m3, s3, a4, m4, s4, a5, m5, s5, a6, m6, s6):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2) \
                    + a4*np.exp(-(x-m4)**2/2/s4**2) + a5*np.exp(-(x-m5)**2/2/s5**2) + a6*np.exp(-(x-m6)**2/2/s6**2)
                return f
            return _gauss6
        elif cn == 7:
            def _gauss7(x, a1, m1, s1, a2, m2, s2, a3, m3, s3, a4, m4, s4, a5, m5, s5, a6, m6, s6, a7, m7, s7):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2) \
                    + a4*np.exp(-(x-m4)**2/2/s4**2) + a5*np.exp(-(x-m5)**2/2/s5**2) + a6*np.exp(-(x-m6)**2/2/s6**2) \
                    + a7*np.exp(-(x-m7)**2/2/s7**2)
                return f
            return _gauss7
        elif cn == 8:
            def _gauss8(x, a1, m1, s1, a2, m2, s2, a3, m3, s3, a4, m4, s4, a5, m5, s5, a6, m6, s6, a7, m7, s7,
                        a8, m8, s8):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2) \
                    + a4*np.exp(-(x-m4)**2/2/s4**2) + a5*np.exp(-(x-m5)**2/2/s5**2) + a6*np.exp(-(x-m6)**2/2/s6**2) \
                    + a7*np.exp(-(x-m7)**2/2/s7**2) + a8*np.exp(-(x-m8)**2/2/s8**2)
                return f
            return _gauss8
        elif cn == 9:
            def _gauss9(x, a1, m1, s1, a2, m2, s2, a3, m3, s3, a4, m4, s4, a5, m5, s5, a6, m6, s6, a7, m7, s7,
                        a8, m8, s8, a9, m9, s9):
                f = a1*np.exp(-(x-m1)**2/2/s1**2) + a2*np.exp(-(x-m2)**2/2/s2**2) + a3*np.exp(-(x-m3)**2/2/s3**2) \
                    + a4*np.exp(-(x-m4)**2/2/s4**2) + a5*np.exp(-(x-m5)**2/2/s5**2) + a6*np.exp(-(x-m6)**2/2/s6**2) \
                    + a7*np.exp(-(x-m7)**2/2/s7**2) + a8*np.exp(-(x-m8)**2/2/s8**2) + a9*np.exp(-(x-m9)**2/2/s9**2)
                return f
            return _gauss9
        else:
            raise IOError("cn = {}, 'cn' should be in 1 to 9.".format(cn))

    def _fitter(self, y, rms, refit=None, guess=None):
        if refit is None:
            if guess is None:
                raise IOError("'refit' or 'guess' are required.")
            else:
                if len(guess) > 9:
                    guess = np.array(guess)
                    guess = np.sort(guess[np.argsort(y[guess])[len(guess)-9:]])
                cn = len(guess)
                ig = np.zeros((cn, 3))
                bmin = np.zeros((cn, 3))
                bmax = np.zeros((cn, 3))
                # amplitude
                ig[:, 0] = y[guess]
                bmax[:, 0] = y[guess]+rms
                # mean
                ig[:, 1] = self.x[guess]
                bmin[:, 1] = self.x[guess]-self._mlim
                bmax[:, 1] = self.x[guess]+self._mlim
                # stddev
                ig[:, 2] = 0.2/np.sqrt(8*np.log(2))
                bmin[:, 2] = self._smin
                bmax[:, 2] = self._smax
                return curve_fit(self._gauss(cn), self.x, y, ig.flatten(), bounds=(bmin.flatten(), bmax.flatten()),
                                 maxfev=1000*(3*cn+1))
        else:
            refit = refit.reshape(-1, 3)
            cn = len(refit)
            ig = np.zeros((cn, 3))
            bmin = np.zeros((cn, 3))
            bmax = np.zeros((cn, 3))
            iy = interp1d(self.x, y)
            # amplitude
            ig[:, 0] = refit[:, 0]
            bmax[:, 0] = np.max(np.array([iy(refit[:, 1]), refit[:, 0]]), axis=0)+rms
            # mean
            ig[:, 1] = refit[:, 1]
            bmin[:, 1] = refit[:, 1]-self._mlim
            bmax[:, 1] = refit[:, 1]+self._mlim
            # stddev
            ig[:, 2] = np.min(np.array([refit[:, 2], np.full(cn, self._smax)]), axis=0)
            bmin[:, 2] = self._smin
            bmax[:, 2] = self._smax
            try:
                return curve_fit(self._gauss(cn), self.x, y, ig.flatten(), bounds=(bmin.flatten(), bmax.flatten()),
                             maxfev=1000*(3*cn+1))
            except:
                print(ig)
                print(bmin)
                print(bmax)
                raise ValueError("'p0'(initial guess) is out of 'bounds'.")

    def _fitter_after_fof(self, y, rms, refit=None):
        if refit is None:
            raise IOError("'refit' or 'guess' are required.")
        else:
            refit = refit.reshape(-1, 3)
            cn = len(refit)
            ig = np.zeros((cn, 3))
            bmin = np.zeros((cn, 3))
            bmax = np.zeros((cn, 3))
            iy = interp1d(self.x, y)
            # amplitude
            ig[:, 0] = refit[:, 0]
            bmin[:, 0] = refit[:, 0]*0.9-rms
            bmax[:, 0] = refit[:, 0]*1.1+rms
            # mean
            ig[:, 1] = refit[:, 1]
            bmin[:, 1] = refit[:, 1]-self._mlim
            bmax[:, 1] = refit[:, 1]+self._mlim
            # stddev
            ig[:, 2] = refit[:, 2]
            bmin[:, 2] = np.max(np.array([refit[:, 2]*0.5, np.full(cn, self._smin)]), axis=0)
            bmax[:, 2] = np.min(np.array([refit[:, 2]*2, np.full(cn, self._smax)]), axis=0)
            try:
                return curve_fit(self._gauss(cn), self.x[self._rmssize:-self._rmssize], y[self._rmssize:-self._rmssize],
                                 ig.flatten(), bounds=(bmin.flatten(), bmax.flatten()), maxfev=1000*(3*cn+1))
            except:
                print(ig)
                print(bmin)
                print(bmax)
                raise ValueError("'p0'(initial guess) is out of 'bounds'.")

    def initial_fit(self):
        fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd'),
                   dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8'))
        pp = 0
        pi = 0
        ppp = np.nansum(self.det_comp)
        print('\n[ Initial fitting for making initial guess ]')
        print('Progress 0...|10..|20..|30..|40..|50..|60..|70..|80..|90..|100%')
        print('         ', end='')
        for d in range(self.nd):
            for r in range(self.nr):
                if not self.det_comp[d, r]:
                    continue
                fp, fc = self._fitter(self.mdata_for_finding[:, d, r], self.rms_for_finding[d, r], guess=self._px(d, r))
                fp = fp.reshape((-1, 3))
                for i in range(len(fp)):
                    fr.add_row((r, d, len(fp), i, fp[i, 0], fp[i, 1], fp[i, 2]))
                pp += 1
                if int(pp/ppp*50) > pi:
                    pi += 1
                    print('#', end='')
        fr['dv'] = fr['sd']*np.sqrt(8.*np.log(2))
        fr['area'] = fr['tp']*np.sqrt(2.*np.pi*fr['sd']**2.)
        self._fitres1 = fr
        print(' complete!')
        return self._fitres1

    @property
    def initial_fit_result(self):
        if self._fitres1 is None:
            self._fitres1 = self.initial_fit()
        return self._fitres1

    def second_fit(self):
        fit1 = self.initial_fit_result
        fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd'),
                   dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8'))
        pp = 0
        pi = 0
        ppp = np.nansum(self.det_comp)
        print('\n[ Second fitting with initial guess ]')
        print('Progress 0...|10..|20..|30..|40..|50..|60..|70..|80..|90..|100%')
        print('         ', end='')
        for d in range(self.nd):
            for r in range(self.nr):
                if not self.det_comp[d, r]:
                    continue
                ct = fit1[(fit1['rp'] == r) & (fit1['dp'] == d)]
                ig = []
                for i in range(len(ct)):
                    for c in ['tp', 'vp', 'sd']:
                        ig.append(ct[c][i])
                fp, fc = self._fitter(self.y[:, d, r], self.srms[d, r], refit=np.array(ig))
                fp = fp.reshape((-1, 3))
                for i in range(len(fp)):
                    fr.add_row((r, d, len(fp), i, fp[i, 0], fp[i, 1], fp[i, 2]))
                pp += 1
                if int(pp/ppp*50) > pi:
                    pi += 1
                    print('#', end='')
        fr['dv'] = fr['sd']*np.sqrt(8.*np.log(2))
        fr['area'] = fr['tp']*np.sqrt(2.*np.pi*fr['sd']**2.)
        self._fitres2 = fr
        print(' complete!')
        return self._fitres2

    @property
    def second_fit_result(self):
        if self._fitres2 is None:
            self._fitres2 = self.second_fit()
        return self._fitres2

    def _get_continuity(self, fit, d, r, ss):
        ss = int(np.ceil(float(ss)/2))
        cen = fit[(fit['rp'] == r) & (fit['dp'] == d)]
        sur = fit[(fit['rp'] >= r-ss) & (fit['rp'] <= r+ss)]
        sur = sur[(sur['dp'] >= d-ss) & (sur['dp'] <= d+ss)]
        cn = len(cen)
        ig = np.zeros((cn, 3))
        for c in range(cn):
            sc = sur[(sur['vp'] >= cen['vp'][c]-self._mlim) & (sur['vp'] <= cen['vp'][c]+self._mlim)]
            if len(sc) == 0:
                for i, n in enumerate(['tp', 'vp', 'sd']):
                    ig[c, i] = cen[n][c]
            else:
                for i, n in enumerate(['tp', 'vp', 'sd']):
                    ig[c, i] = np.median(sc[n])
        return ig.flatten()

    def final_fit(self):
        fit2 = self.second_fit_result
        fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8'))
        pp = 0
        pi = 0
        ppp = np.nansum(self.det_comp)
        print('\n[ Final fitting for continuity ]')
        print('Progress 0...|10..|20..|30..|40..|50..|60..|70..|80..|90..|100%')
        print('         ', end='')
        for d in range(self.nd):
            for r in range(self.nr):
                if not self.det_comp[d, r]:
                    continue
                ig = self._get_continuity(fit2, d, r, self._spasmo_find)
                fp, fc = self._fitter(self.y[:, d, r], self.srms[d, r], refit=ig)
                fp = fp.reshape((-1, 3))
                for i in range(len(fp)):
                    fr.add_row((r, d, len(fp), i, fp[i, 0], fp[i, 1], fp[i, 2], fp[i, 0]/self.rms[d, r]))
                pp += 1
                if int(pp/ppp*50) > pi:
                    pi += 1
                    print('#', end='')
        fr['dv'] = fr['sd']*np.sqrt(8.*np.log(2))
        fr['area'] = fr['tp']*np.sqrt(2.*np.pi*fr['sd']**2.)
        self._fitres3 = fr
        print(' complete!')
        return self._fitres3

    @property
    def final_fit_result(self):
        if self._fitres3 is None:
            self._fitres3 = self.final_fit()
        return self._fitres3

    def fit_after_fof(self, fit):
        fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd'),
                   dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8'))
        pp = 0
        pi = 0
        ppp = np.nansum(self.det_comp)
        print('\n[ Second fitting with initial guess ]')
        print('Progress 0...|10..|20..|30..|40..|50..|60..|70..|80..|90..|100%')
        print('         ', end='')
        for d in range(self.nd):
            for r in range(self.nr):
                if not self.det_comp[d, r]:
                    continue
                ct = fit[(fit['rp'] == r) & (fit['dp'] == d)]
                if len(ct) == 0:
                    continue
                ig = []
                for i in range(len(ct)):
                    for c in ['tp', 'vp', 'sd']:
                        ig.append(ct[c][i])
                fp, fc = self._fitter(self.y[:, d, r], self.srms[d, r], refit=np.array(ig))
                fp = fp.reshape((-1, 3))
                for i in range(len(fp)):
                    fr.add_row((r, d, len(fp), i, fp[i, 0], fp[i, 1], fp[i, 2]))
                pp += 1
                if int(pp/ppp*50) > pi:
                    pi += 1
                    print('#', end='')
        fr['dv'] = fr['sd']*np.sqrt(8.*np.log(2))
        fr['area'] = fr['tp']*np.sqrt(2.*np.pi*fr['sd']**2.)
        self._fitres_after_fof = fr
        print(' complete!')
        return self._fitres_after_fof

    @property
    def fit_result_after_fof(self):
        if self._fitres_after_fof is None:
            print('Run Decompose.fit_after_fof() first.')
            return
        else:
            return self._fitres_after_fof

    def run_decompose(self, save=None):
        """
        Gaussian decomposing for cube data.
        This method include three following steps of multi-Gaussian fitting.
            Decompose.initial_fit()
            Decompose.second_fit()
            Decompose.final_fit()
        The result of decomposing will be given by astropy table.

        Parameters:
            save: path or filename
                Filename of table data.

        Returns:
             astropy.table.Table
        """
        self._fitres3 = self.final_fit()
        if save is not None:
            if isinstance(save, str):
                write(self._fitres3, save, overwrite=False)
        return self._fitres3

    @property
    def decompose_result(self):
        """
        Return the Gaussian decomposing result.
        The columns of table are ...
            'rp' : R.A. pixel number (from 0, python style)
            'dp' : Dec. pixel number (from 0, python style)
            'tn' : Total number of components of this pixel.
            'cn' : Component number (order, starting at low speed) in this pixel.
            'tp' : Peak temperature (amplitude of Gaussian func.)
            'vp' : Velocity position (mean value of Gaussian func.)
            'sd' : Sigma (standard deviation of Gaussian func.)
            'dv' : Line width (FWHM of Gaussian func.)

        Returns:
             astropy.table.Table
        """
        if self._fitres3 is None:
            self._fitres3 = self.final_fit()
        return self._fitres3

    def plot_fit(self, fit, vr=None, yr=None):
        """
        Interactive scan for fitting results.
        Plot fitting results over line profile on click position.

        Parameters:
            fit : astropy.table.Table
                Result table of decomposing.
                    Decompose.initial_fit_result
                    Decompose.second_fit_result
                    Decompose.final_fit_result
                    Decompose.decompose_result
            vr : list or tuple, optional
                Range of velocity axis for line plot.
            yr : list or tuple, optional
                Range of temperature axis for line plot.

        Returns:
            tuple (plt.figure, plt.subplot, plt.subplot)
        """
        plt.ion()
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        x = self.x
        if not isinstance(fit, Table):
            raise TypeError("'fit' is not astropy Table.".format(fit))
        if vr is None:
            vr = [x[self._rmssize], x[-self._rmssize]]
        elif isinstance(vr, list) or isinstance(vr, tuple):
            pass
        else:
            raise TypeError("'vr' is not list or tuple.".format(vr))
        if yr is None:
            yr = [-4*np.nanmedian(self.rms), np.nanmax(self.data)+np.nanmedian(self.rms)]
        elif isinstance(yr, list) or isinstance(yr, tuple):
            pass
        else:
            raise TypeError("'yr' is not list or tuple".format(yr))
        # make map
        det = np.full(self.data.shape, 0.)
        det[np.nan_to_num(self.data) > 3*np.nanmedian(self.rms)] = 1.
        dsum = np.nansum(det, axis=(1, 2))
        imap = np.sum(self.data[dsum > 2*np.mean(dsum)], axis=0)
        # map drawing
        lsfig = plt.figure(figsize=(16, 9))
        lsmap = lsfig.add_axes([0.06, 0.07, 0.36, 0.9], projection=self.wcs2d)
        lsmap.set_xlabel('R.A.')
        lsmap.set_ylabel('Dec.')
        lsmap.imshow(imap, origin='lower', vmin=0, cmap='inferno')
        self._pp = lsmap.scatter(-100, -100, s=60, color='lime')
        lsmap.set_xlim(-0.5, self.nr-0.5)
        lsmap.set_ylim(-0.5, self.nd-0.5)
        # line plot
        lbl = lsfig.add_axes([0.47, 0.07, 0.17, 0.3])
        lcl = lsfig.add_axes([0.47, 0.37, 0.17, 0.3])
        ltl = lsfig.add_axes([0.47, 0.67, 0.17, 0.3])
        lbc = lsfig.add_axes([0.64, 0.07, 0.17, 0.3])
        lslin = lsfig.add_axes([0.64, 0.37, 0.17, 0.3])
        ltc = lsfig.add_axes([0.64, 0.67, 0.17, 0.3])
        lbr = lsfig.add_axes([0.81, 0.07, 0.17, 0.3])
        lcr = lsfig.add_axes([0.81, 0.37, 0.17, 0.3])
        ltr = lsfig.add_axes([0.81, 0.67, 0.17, 0.3])

        def _onclick(event):
            if event.inaxes == lsmap.axes:
                self._pr = int(round(event.xdata))
                self._pd = int(round(event.ydata))
            elif event.inaxes == lbl.axes:
                self._pr -= 1
                self._pd -= 1
            elif event.inaxes == lcl.axes:
                self._pr -= 1
            elif event.inaxes == ltl.axes:
                self._pr -= 1
                self._pd += 1
            elif event.inaxes == lbc.axes:
                self._pd -= 1
            elif event.inaxes == ltc.axes:
                self._pd += 1
            elif event.inaxes == lbr.axes:
                self._pr += 1
                self._pd -= 1
            elif event.inaxes == lcr.axes:
                self._pr += 1
            elif event.inaxes == ltr.axes:
                self._pr += 1
                self._pd += 1
            else:
                return
            r = self._pr
            d = self._pd
            ri = [-1, -1, -1, 0, 0, 0, +1, +1, +1]
            di = [-1, 0, +1, -1, 0, +1, -1, 0, +1]
            for i, ax in enumerate([lbl, lcl, ltl, lbc, lslin, ltc, lbr, lcr, ltr]):
                ax.clear()
                ax.set_xlim(*vr)
                ax.set_ylim(*yr)
                if not i in [0, 1, 2]:
                    ax.set_yticklabels([])
                if not i in [0, 3, 6]:
                    ax.set_xticklabels([])
                rr = r+ri[i]
                dd = d+di[i]
                self._pp.set_offsets([rr-1, dd-1])
                if 0 <= rr < self.nr and 0 <= dd < self.nd:
                    srms = self._snr*self.rms_for_finding[dd, rr]
                    rms = self._snr*self.srms[dd, rr]
                    ax.plot([x[0], x[-1]], [0, 0], lw=1, alpha=0.5, color='black')
                    ax.plot([x[0], x[-1]], [srms, srms], ls='dotted', lw=1, alpha=0.3, color='blue')
                    ax.plot([x[0], x[-1]], [rms, rms], ls='dotted', lw=1, alpha=0.5, color='red')
                    ax.step(x, self.mdata_for_finding[:, dd, rr], color='blue', lw=1, alpha=0.6)
                    ax.step(x, self.sdata[:, dd, rr], color='black', lw=1)
                    ct = fit[(fit['rp'] == rr) & (fit['dp'] == dd)]
                    if len(ct) == 0:
                        continue
                    cp = self._px(dd, rr)
                    for c in cp:
                        xx = [x[c], x[c]]
                        yy = [yr[0]/2, 0]
                        ax.plot(xx, yy, color='red', lw=1, alpha=0.5)
                    cm = Gaussian1D(ct['tp'][0], ct['vp'][0], ct['sd'][0])
                    for c in range(1, len(ct)):
                        cm += Gaussian1D(ct['tp'][c], ct['vp'][c], ct['sd'][c])
                    if len(ct) == 1:
                        ax.plot(x, cm(x), color='green', lw=1, alpha=0.5)
                    else:
                        for c in range(len(ct)):
                            ax.plot(x, cm[c](x), color='green', lw=1, alpha=0.5)
                    ax.plot(x, cm(x), color='red', lw=2)

            lbc.set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
            lcl.set_ylabel(r'$T_\mathrm{A}^* [K]$')
            c = self.wcs2d.pixel_to_world(r, d)
            c = c.to_string('hmsdms').split(' ')
            lslin.annotate('R.A.: {} ({})'.format(r, c[0]), (0.03, 0.92), xycoords='axes fraction')
            lslin.annotate('Dec.: {} ({})'.format(d, c[1]), (0.03, 0.86), xycoords='axes fraction')
            lslin.figure.canvas.draw()

        lsfig.canvas.mpl_connect('button_press_event', _onclick)
        return lsfig, lsmap, lslin


