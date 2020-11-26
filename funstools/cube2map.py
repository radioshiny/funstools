from .checkshape import check_all, make_velo
from .makemaps import get_rms, get_mask, get_det, wcs2d, header2d
from .smooth import smooth3d
import numpy as np
from warnings import warn
from matplotlib import pyplot as plt


class Cube2map:
    """
    The Cube2map is a class intended to make the RMS error, 1st or 2nd moment, or channel maps
    from the FITS cube data with the desired level of noise masking and smoothing.
    """

    def __init__(self, cube=None, ext=None, getrms='both', rmssize=None, max_rms=None,
                 snr=3., velocity_smo=1, spatial_smo=1, masking=True, smoothing=True):
        """
        Construct a 'Cube2map' object.

        Parameters:
            cube : file path, hdulist, or hdu
                Name or location of the input cube data.
            ext : int, optional
                Extension number in HDUlist with multiple HDU.
            getrms : {'both', 'left', 'right'}
                Location of channels to use for RMS error measurement. Default is 'both'.
            rmssize : int, optional
                Number of channels to use for RMS error measurement.
            max_rms : float, optional
                Make maps excluding pixels that have rms error higher than the 'max_rms'.
                If mas_rms is None, use a value that is twice the detectable level.
            snr : float
                Masking data having S/N ratio less than 'snr' value. Default is 3.
            velocity_smo: float
                FWHM channel size of the 1D Gaussian kernel for velocity smoothing. Default is 1.
            spatial_smo: float
                FWHM pixel size of the 2D Gaussian kernel for spatial smoothing. Default is 1.
            masking : bool
                Make maps with noise masking. Default is True
            smoothing : bool
                Make maps with smoothing. Default is True
        """
        self.cube = check_all(cube, 3, ext, verbose=True)
        self._header = self.cube.header
        self._data = self.cube.data
        self._nr = self.cube.header['NAXIS1']
        self._nd = self.cube.header['NAXIS2']
        self._nx = self.cube.header['NAXIS3']
        self._cw = self.cube.header['CDELT3']
        self._getrms = str(getrms)
        if rmssize is None:
            self._rmssize = int(self._nx/3)
        else:
            self._rmssize = rmssize
        if self._getrms == 'both':
            self._rmsch = np.arange(self._nx)[:self._rmssize]
            self._rmsch = np.append(self._rmsch, np.arange(self._nx)[-self._rmssize:])
        elif self._getrms == 'left':
            self._rmsch = np.arange(self._nx)[:self._rmssize]
        elif self._getrms == 'right':
            self._rmsch = np.arange(self._nx)[-self._rmssize:]
        else:
            raise ValueError("'{}' is not recognized. Possible options: 'left', 'right', or 'both'.".format(self._getrms))
        self._maxrms = max_rms
        self._masking = masking
        self._smoothing = smoothing
        self._snr = float(snr)
        self._velsmo = float(velocity_smo)
        self._spasmo = float(spatial_smo)


    _rms = None
    _mask = None
    _det = None
    _mdata = None
    _sdata = None
    _srms = None
    _smdata = None
    _x = None
    _y = None
    _m0 = None
    _m1 = None
    _m2 = None
    _wcs2d = None
    _header2d = None
    _pr = 0
    _pd = 0

    @property
    def x(self):
        """
        Return velocity array of channels. X-data for plotting a line profile.
        """
        if self._x is None:
            self._x = make_velo(self.cube)
        return self._x

    @property
    def header(self):
        """
        Return FITS header removed unused axes and histories.
        """
        return self._header

    @property
    def data(self):
        """
        Return original FITS data without masking and smoothing.
        """
        return self._data

    @property
    def nr(self):
        """
        Number of pixels of R.A. axis.
        """
        return self._nr

    @property
    def nd(self):
        """
        Number of pixels of Dec. axis.
        """
        return self._nd

    @property
    def nx(self):
        """
        Number of channels of velocity axis.
        """
        return self._nx

    @property
    def cw(self):
        """
        Velocity channel width.
        """
        return self._cw

    @property
    def rms(self):
        """
        Return a RMS error map of cube data based on initial parameters, 'getrms' and 'rmssize'.
        """
        if self._rms is None:
            self._rms = get_rms(self._data, self._getrms, self._rmssize)
        return self._rms

    @property
    def mask(self):
        """
        Test channels for detectable and return result as a 3d-array.
        Noise channels: 0. and detectable channels: 1.

        Returns:
            3d-array (cube) filled with 0. or 1.
        """
        if self._mask is None:
            self._mask = get_mask(self._data, self._snr, self.rms, self._maxrms, self._velsmo, verbose=True)
        return self._mask

    @property
    def det(self):
        """
        Test pixels for detectable and return result as a boolean 2d-array.
        Noisy or empty pixels: False and detectable pixels: True.

        Returns:
            boolean 2d-array (image)
        """
        if self._det is None:
            self._det = get_det(self.mask)
        return self._det

    @property
    def detmask(self):
        """
        Return Cube2map.det to mask.
        """
        detmask = np.full(self.det.shape, np.nan)
        detmask[self.det] = 1.
        return detmask

    @property
    def mdata(self):
        """
        Masked cube data. mdata = mask * data
        """
        if self._mdata is None:
            self._mdata = self.mask*self._data
        return self._mdata

    @property
    def sdata(self):
        """
        Smoothed cube data based on initial parameters, 'velocity_smo' and 'spatial_smo'.
        """
        if self._sdata is None:
            self._sdata = smooth3d(self._data, self._spasmo, self._velsmo)
        return self._sdata

    @property
    def srms(self):
        """
        RMS error of smoothed data.
        """
        if self._srms is None:
            self._srms = get_rms(self.sdata, self._getrms, self._rmssize)
        return self._srms

    @property
    def smdata(self):
        """
        Data cube with masking and smoothing.
        """
        if self._smdata is None:
            self._smdata = smooth3d(self.mdata, self._spasmo, self._velsmo)
        return self._smdata

    @property
    def y(self):
        """
        Data cube to be used when making maps based on initial options, 'masking' and 'smoothing'.

        Masking   Smoothing   Return
          True      True        smdata
          True      False       mdata
          False     True        sdata
          False     False       data
        """
        if self._y is None:
            if self._smoothing:
                if self._masking:
                    self._y = self.smdata
                else:
                    self._y = self.sdata
            else:
                self._y = self.mdata
        return self._y

    @property
    def m0(self):
        """
        Return the result of the most recently computed moment0 method.

        If there is no pre-computed moment 0 map,
        return the result of the moment0 method using the entire velocity range,
        excluding the channel used to measure the rms error.
        """
        if self._m0 is None:
            self._m0 = self.moment0(cr=(self._rmssize, -self._rmssize), verbose=False)
        return self._m0

    @property
    def m1(self):
        """
        Return the result of the most recently computed moment1 method.

        If there is no pre-computed moment 1 map,
        return the result of the moment1 method using the entire velocity range,
        excluding the channel used to measure the rms error.
        """
        if self._m1 is None:
            self._m1 = self.moment1(cr=(self._rmssize, -self._rmssize))
        return self._m1

    @property
    def m2(self):
        """
        Return the result of the most recently computed moment2 method.

        If there is no pre-computed moment 2 map,
        return the result of the moment2 method using the entire velocity range,
        excluding the channel used to measure the rms error.
        """
        if self._m2 is None:
            self._m2 = self.moment2(cr=(self._rmssize, -self._rmssize))
        return self._m2

    @property
    def wcs2d(self):
        """
        Return the 2D wcs information for making plot with WCSAxes in matplotlib.
        More information: http://docs.astropy.org/en/stable/visualization/wcsaxes/
        """
        if self._wcs2d is None:
            self._wcs2d = wcs2d(self.cube)
        return self._wcs2d

    @property
    def header2d(self):
        """
        Return the 2D fits header.
        """
        if self._header2d is None:
            self._header2d = header2d(self.cube)
        return self._header2d

    def v2ch(self, vr):
        """
        Convert velocity range to channel range.

        Parameters:
            vr : float, list or tuple
                Velocity value or range, [start, end] or (start, end).

        Returns:
            cr : list
                Channel value or range, [start, end+1] ('+1' follows the python indexing rule).
        """
        if vr is None:
            if self._rmssize is None:
                size = int(len(self.x)/3)
            else:
                size = int(self._rmssize)
            if self._getrms == 'left':
                return [size, self.nx]
            elif self._getrms == 'right':
                return [0, self.nx-size]
            else:
                return [size, self.nx-size]
        elif isinstance(vr, list) or isinstance(vr, tuple):
            if self.x[0] > vr[0] or self.x[-1] < vr[1]:
                warn('Given velocity range exceeds the velocity axis of cube data.')
            return [np.argmin(np.abs(self.x-float(vr[0]))), np.argmin(np.abs(self.x-float(vr[1])))+1]
        elif isinstance(vr, float) or isinstance(vr, int):
            if self.x[0] > vr or self.x[-1] < vr:
                warn('Given velocity range exceeds the velocity axis of cube data.')
            return np.argmin(np.abs(self.x-float(vr)))
        else:
            raise TypeError("'vr' (velocity range) is not a list or tuple.")

    def moment0(self, vr=None, cr=None, verbose=False):
        """
        Make moment 0 map for given velocity or channel range.

        Parameters:
            vr : list or tuple
                Velocity range [v1, v2]. Default is None.
            cr : list or tuple
                Channel number range [c1, c2], Default is None.
            verbose : bool
                Return sigma of moment0, Default is False.

            If both 'vr' and 'cr' is None, 'vr' is set entire velocity range
                except for the channel used to measure the RMS error.
            If both 'vr' and 'cr' are given, 'cr' is used first.
            If only one range is given without a key name, it is used as velocity.

        Returns:
            moment0 : 2d-array (image)

            (optional)
            sigma_moment0 : float
        """
        if cr is None:
            cr = self.v2ch(vr)
        if cr[1] < 0:
            cr = (cr[0], np.arange(self._nx)[cr[1]])
        if not (isinstance(cr, list) or isinstance(cr, tuple)):
            raise TypeError("'cr' (channel range) is not a list or tuple.")
        self._m0 = np.sum(self.y[cr[0]:cr[1]]*self.cw, axis=0)*self.detmask
        if verbose:
            dch = np.nansum(self.mask, axis=(1, 2))
            rms_dch = np.nanmedian(dch[self._rmsch])
            nch = np.sum(dch[cr[0]:cr[1]] > self._snr*rms_dch)
            if self._smoothing:
                mrms = np.nanmedian(self.srms*self.detmask)
            else:
                mrms = np.nanmedian(self.rms*self.detmask)
            m0rms = np.sqrt(nch)*mrms*self.cw
            print('\n[ Making Moment 0 (integrated intensity) map ]')
            print('Channel range      = {:d} ~ {:d}'.format(cr[0], cr[1]))
            print('Velocity range     = {:.2f} ~ {:.2f}'.format(self._x[cr[0]], self._x[cr[1]]))
            print('N total channel    = {:d}'.format(cr[1]-cr[0]))
            print('N detected channel = {:d}'.format(nch))
            print('Smoothing          = {}'.format(self._smoothing))
            print('Median RMS_line    = {:.3f} K'.format(mrms))
            print('-----\nRMS_moment0        = {:.3f} K km/s'.format(m0rms))
            return self._m0, m0rms
        else:
            return self._m0

    def moment1(self, vr=None, cr=None):
        """
        Make moment 1 map for given velocity or channel range.

        Parameters:
            vr : list or tuple
                Velocity range [v1, v2]. Default is None.
            cr : list or tuple
                Channel number range [c1, c2], Default is None.

            If both 'vr' and 'cr' is None, 'vr' is set entire velocity range
                except for the channel used to measure the RMS error.
            If both 'vr' and 'cr' are given, 'cr' is used first.
            If only one range is given without a key name, it is used as velocity.

        Returns:
            moment1 : 2d-array (image)
        """
        vel = np.ones_like(self.data)*self.x[:, np.newaxis, np.newaxis]
        if cr is None:
            cr = self.v2ch(vr)
        if isinstance(cr, list) or isinstance(cr, tuple):
            m0 = self.moment0(vr, cr)
            m0[m0 == 0.] = np.nan
            m1 = np.sum(vel[cr[0]:cr[1]]*self.y[cr[0]:cr[1]]*self.cw, axis=0)/m0*self.detmask
            m1[np.nan_to_num(m1) < self.x[cr[0]]] = np.nan
            m1[np.nan_to_num(m1) > self.x[cr[1]]] = np.nan
            self._m1 = m1
            # m1[np.nan_to_num(m1) > np.nanpercentile(m1, 99.5)] = np.nan
            # m1[np.nan_to_num(m1) < np.nanpercentile(m1, 0.5)] = np.nan
            return self._m1
        else:
            raise TypeError("'cr' (channel range) is not a list or tuple.")

    def moment2(self, vr=None, cr=None):
        """
        Make moment 2 map for given velocity or channel range.

        Parameters:
            vr : list or tuple
                Velocity range [v1, v2]. Default is None.
            cr : list or tuple
                Channel number range [c1, c2], Default is None.

            If both 'vr' and 'cr' is None, 'vr' is set entire velocity range
                except for the channel used to measure the RMS error.
            If both 'vr' and 'cr' are given, 'cr' is used first.
            If only one range is given without a key name, it is used as velocity.

        Returns:
            moment2 : 2d-array (image)
        """
        vel = np.ones_like(self.data)*self.x[:, np.newaxis, np.newaxis]
        if cr is None:
            cr = self.v2ch(vr)
        m0 = self.moment0(vr, cr)
        m1 = self.moment1(vr, cr)
        m0[m0 == 0] = np.nan
        self._m2 = np.sqrt(np.sum(self.y[cr[0]:cr[1]]*self.cw*(vel[cr[0]:cr[1]]-m1)**2., axis=0)/m0)*self.detmask
        return self._m2

    def tpeak(self, vr=None, cr=None):
        if cr is None:
            cr = self.v2ch(vr)
        if isinstance(cr, list) or isinstance(cr, tuple):
            return np.max(self.y[cr[0]:cr[1]], axis=0)
        else:
            raise TypeError("'cr' (channel range) is not a list or tuple.")

    def chmap(self, mn=9, vr=None, cr=None):
        """
        Make channel maps and return array of maps and list of labels.

        Parameters:
            mn : int
                Number of channel maps to generate. Default is 9. Must be non-negative.
            vr, cr: same as in the moment method.

        Returns:
            maps, labels : tuple, (3d-array, list of str)
                maps.shape = (mn, nd, nr)
                maps[i] = 2d-array (image) of moment 0 map of a 'i'st channel section.
                labels[i] = string indicating velocity range of a 'i'st channel map.
        """
        if cr is None:
            cr = self.v2ch(vr)
        if cr[1]-cr[0] < mn:
            raise ValueError("Channel range is too small")
        ci = np.linspace(cr[0], cr[1], mn+1, endpoint=True, dtype=int)
        cm = np.zeros((mn, self.nd, self.nr))
        for i in range(mn):
            cm[i] = self.moment0(cr=[ci[i], ci[i+1]])
        # normalizing intensities because number of channel values are different.
        cn = ci[1:]-ci[:-1]
        weight = float(ci[-1]-ci[0])/float(mn)/cn.astype(float)
        cm = cm*weight[:, np.newaxis, np.newaxis]
        vv = self.x[ci]
        return cm, [r'${} \leq v_\mathrm{{ch}} < {}$'.format(vv[i], vv[i+1]) for i in range(mn)]

    def line_scan(self, vr=None, yr=None):
        """
        Interactive scan for line profiles.
        Plot line profile on click position.

        Parameters:
            vr : list or tuple, optional
                Range of velocity axis for line plot.
            yr : list or tuple, optional
                Range of temperature axis for line plot.
            cut : list [ra_start, dec_start, ra_size, dec_size], optional
                Position and size of map cut in axis fraction.

        Returns:
            tuple (plt.figure, plt.subplot, plt.subplot)
        """
        plt.ion()
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        x = self.x
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
                if 0 <= rr < self.nr and 0 <= dd < self.nd:
                    self._pp.set_offsets([rr-1, dd-1])
                    srms = self._snr*self.srms[dd, rr]
                    rms = self._snr*self.rms[dd, rr]
                    ax.plot([x[0], x[-1]], [0, 0], lw=1, alpha=0.5, color='black')
                    ax.plot([x[0], x[-1]], [srms, srms], ls='dotted', lw=1, alpha=0.3, color='blue')
                    ax.plot([x[0], x[-1]], [rms, rms], ls='dotted', lw=1, alpha=0.5, color='red')
                    ax.step(x, self.sdata[:, dd, rr], color='blue', lw=1, alpha=0.6)
                    ax.step(x, self.data[:, dd, rr], color='red', lw=1)

            lbc.set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
            lcl.set_ylabel(r'$T_\mathrm{A}^* [K]$')
            c = self.wcs2d.pixel_to_world(r, d)
            lslin.annotate('Coord.: {}'.format(c.to_string('hmsdms')), (0.03, 0.92), xycoords='axes fraction')
            lslin.annotate('Pixel: [{}, {}]'.format(d, r), (0.03, 0.86), xycoords='axes fraction')
            lslin.figure.canvas.draw()

        lsfig.canvas.mpl_connect('button_press_event', _onclick)
        return lsfig, lsmap, lslin

