from .checkshape import get_fits, isfits, check_axis
from .smooth import smooth1d
from .userinput import yesno
from warnings import warn
import numpy as np
from os.path import isfile
from itertools import groupby
from operator import itemgetter
from astropy.io import fits
from astropy.wcs import WCS


def get_rms(data, where='both', size=None, ext=None):
    """
    Measure rms noise level from spectrum and return.

    Parameters:
        data : file path, object, or ndarray
            Name or location of data to measure rms.
        where : {'both', 'left', 'right'}
            RMS is measured from the channels in the selected direction.
            Default is 'both'
        size : int
            Size of channel to use for rms measurement.
            Default is 1/3 of channels.
        ext : int, optional
            Extension number when input FITS is multiple HDUList.

    Returns:
        ndarray
    """
    temp = get_fits(data, 'data', ext)
    if size is None:
        size = int(len(temp)/3)
    elif (size > len(temp)/2 and where == 'both') or (size > len(temp)):
        raise ValueError('RMS channel size is too large.')
    else:
        size = int(size)
    if where == 'left':
        return np.std(temp[:size], axis=0)
    elif where == 'right':
        return np.std(temp[-size:], axis=0)
    elif where == 'both':
        return np.std(np.append(temp[:size], temp[-size:], axis=0), axis=0)
    else:
        raise ValueError("'{}' is not recognized. Possible options: 'left', 'right', or 'both'.".format(where))


# def get_mask(data, snr=3., rms=None, max_rms=None, velocity_smo=1, ext=None, ext_rms=None, verbose=True):
def get_mask(data, snr=3., rms=None, max_rms=None, ext=None, ext_rms=None, getrms='both', rmssize=None, verbose=True, debug=False):
    """
    Check the detection of emission lines for each channel in the cube data
    and return a mask that makes the undetected channels zero.

    Parameters:
        data : file path, object, or ndarray
            Name or location of data to check the detection.
        snr : float
            Signal-to-noise ratio for detection.
        rms : file path, object, or ndarray, optional
            Name or location of rms noise level.
            If no data is given, estimated from the left and right 1/3 channels.
        max_rms : float
            Maximum limit of rms. Pixels with worse rms are excluded.
            Default is '2 x snr x median(rms)'
        # velocity_smo : float
        #     Channel smoothing size for detection.
        ext : int, optional
            Extension number when input FITS is multiple HDUList.
        ext_rms : int, optional
            Extension number when input rms FITS is multiple HDUList.
        verbose : bool
            Print summary of detection.

    Returns:
        ndarray
    """
    temp = get_fits(data, 'data', ext)
    if len(temp.shape) != 3:
        raise TypeError('Input data is not 3D-cube.')
    if rms is None:
        if rmssize is None:
            warn('No RMS data. RMS is estimated from the left and right 1/3 channels.')
            rms = get_rms(temp)
        else:
            rms = get_rms(temp, where=getrms, size=rmssize)
    else:
        if isinstance(rms, float):
            rms = np.full(temp.shape, 1.)*rms
        else:
            rms = get_fits(rms, 'data', ext_rms)
            if not temp[0].shape == rms.shape:
                raise TypeError('Input RMS does not match the shape of input data.')
    snr = float(snr)
    # if max_rms is None:
    #     max_rms = np.nanmedian(rms)*snr*2.
    # else:
    #     max_rms = float(max_rms)
    if max_rms is not None:
        max_rms = float(max_rms)
    mask = np.zeros_like(temp)
    # smodata = smooth1d(temp, velocity_smo)
    # remove smoothing in make_mask, just use input data
    empty, noisy = 0, 0
    for d in range(rms.shape[0]):
        for r in range(rms.shape[1]):
            if np.isnan(rms[d, r]):
                empty += 1
                mask[:, d, r] = np.nan
                continue
            if max_rms is not None:
                if rms[d, r] > max_rms:
                    noisy += 1
                    mask[:, d, r] = np.nan
                    continue
            # y = smodata[:, d, r]
            y = temp[:, d, r]
            pc = np.where(y > 0.)[0]
            for k, g in groupby(enumerate(pc), lambda x: x[0]-x[1]):
                g = list(map(itemgetter(1), g))
                if len(g) > 3 and any(y[g] > snr*rms[d, r]):
                    mask[g, d, r] = 1.
    obs = np.isfinite(temp[0])

    npx = np.nansum(mask, axis=(1, 2))
    if rmssize is None:
        size = int(len(temp)/3)
        medch = np.mean(np.append(npx[:size], npx[-size:]))
        stdch = np.std(np.append(npx[:size], npx[-size:]))
    else:
        medch = np.mean(np.append(npx[:rmssize], npx[-rmssize:]))
        stdch = np.std(np.append(npx[:rmssize], npx[-rmssize:]))
    detch = npx-medch > stdch*snr
    mask[~detch] = 0.
    mask[:, ~obs] = np.nan

    nch = np.nansum(mask, axis=0)
    mask[:, nch < 4] = 0.

    nch = np.nansum(mask, axis=0)
    sumpx = np.nansum(temp*mask, axis=0)
    rmspx = rms*np.sqrt(nch)
    detpx = sumpx > rmspx*snr
    mask[:, ~detpx] = 0.
    mask[:, ~obs] = np.nan

    if debug:
        print('\n===== [ DEBUG - START ] =====')
        # np.save('debug_get_mask_nch.temp', nch)
        _nodet = sumpx < rmspx*snr
        print('tot / det / no_det = {} / {} / {}'.format(np.sum(np.isfinite(sumpx)), np.sum(detpx), np.sum(_nodet)))
        # print('minimum_nch_det = {}'.format(min_nch_det))
        # _check = np.logical_and(det, nch < 4)
        # print('maximum_nch_no_det = {}'.format(max_nch_nod))
        from matplotlib import pyplot as plt
        _fig, _ax = plt.subplots(1, 3)
        _ax[0].hist(nch.flatten(), bins=np.arange(1, np.nanmax(nch)+1))
        # _ax[0].axvline(min_nch_det, color='red')
        _ax[1].imshow(sumpx/rmspx, origin='lower', interpolation='nearest')
        # _ax[2].imshow(_check, origin='lower', interpolation='nearest')
        _ax[2].step(np.arange(len(npx)), npx, where='mid')
        _ax[2].axhline(medch+stdch*snr, color='red')
        _ax[2].scatter(np.arange(len(npx))[detch], npx[detch], c='black')
        print('===== [ DEBUG - END ] =====\n')

    if verbose:
        print('\n[ Masking noise channels ]')
        print('Total pixels        = {:d}'.format(np.prod(rms.shape)))
        print('Empty pixels        = {:d} ; NaN pixel'.format(empty))
        print('Noise pixels        = {:d} ; > max_rms limit'.format(noisy))
        print('Undetected pixels   = {:d} ; < s/n ratio cut'.format(np.prod(rms.shape)-empty-noisy-detpx.sum()))
        print('Available pixels    = {:d}'.format(detpx.sum()))
        tc, dc = np.nansum(temp*0.+1.), np.nansum(mask)
        print('-----\nMasked channels     = {:d}'.format(int(tc-dc)))
        print('Detectable channels = {:d} / {:d} ( {:.1f} % )'.format(int(dc), int(tc), dc/tc*100.))

    return mask


def get_det(mask):
    """
    Return the result of detection as a boolean array.

    Parameters:
        mask : ndarray
            Detection mask returned from 'get_mask'.

    Returns:
        ndarray (image, boolean type)
    """
    det = np.nansum(mask, axis=0)
    # det[det < 4] = np.nan
    # det = det > np.nanpercentile(det, 1)
    return det > 3


def save_fits(file, data, header=None, overwrite=None):
    """
    Save data and header to FITS file.

    Parameters:
        file : file path and name
            File name to save.
        data : ndarray
            Cube or image data to save.
        header : fits.header
            Header of data
        overwrite : bool, optional
            Overwrite?
    """
    done = {False: '{}: File was saved.'.format(file), True: '{}: File was overwritten.'.format(file)}
    if isfile(file):
        if overwrite is None:
            overwrite = yesno('{}: File exists. Overwrite?'.format(file), 'no')
        if not overwrite:
            print('{}: File exists.'.format(file))
            return None
    else:
        overwrite = False
    if isfits(data):
        data.writeto(file, overwrite)
        return done[overwrite]
    else:
        if header is None:
            raise IOError('No header.')
        else:
            if isinstance(header, fits.Header):
                fits.PrimaryHDU(data, header).writeto(file, overwrite=overwrite)
                return done[overwrite]
            else:
                raise TypeError('{}: can not find fits.Header.'.format(header))


def wcs2d(hdu):
    """
    Return 2D wcs information from 2D or 3D FITS header.
    This is a required input when drawing a map with matplotlib
    using the wcs coordinate system.

    Parameters:
        hdu : fits.HDUList, PrimaryHDU, or ImageHDU

    Returns:
        wcs.WCS (with 2 axis)
    """
    if isinstance(hdu, fits.Header):
        temp = np.zeros(tuple(hdu['NAXIS'+str(i)] for i in range(hdu['NAXIS'], 0, -1)))
        hdu = fits.PrimaryHDU(temp, hdu)
    temp = check_axis(hdu, 2)
    return WCS(get_fits(temp, 'header'))


def header2d(hdu):
    """
    Return 2D fits header from 2D or 3D FITS header.

    Parameters:
        hdu : fits.HDUList, PrimaryHDU, or ImageHDU

    Returns:
        fits.Header (with only 2axis)
    """
    if isinstance(hdu, fits.Header):
        temp = np.zeros(tuple(hdu['NAXIS'+str(i)] for i in range(hdu['NAXIS'], 0, -1)))
        hdu = fits.PrimaryHDU(temp, hdu)
    temp = check_axis(hdu, 2)
    return get_fits(temp, 'header')

