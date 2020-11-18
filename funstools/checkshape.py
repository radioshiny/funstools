from astropy.io import fits
from astropy import units as u
from astropy import constants as co
import numpy as np
from warnings import warn


def isfits(data):
    """
    Return True if data is a regular FITS image file (HDUList, PrimaryHDU, or ImageHDU).

    Parameters:
        data (file path or object): Name or location to check if data is a FITS image.

    Returns:
        (bool)
    """
    if isinstance(data, str):
        try:
            fits.open(data)
            return True
        except:
            return False
    elif isinstance(data, fits.PrimaryHDU) or isinstance(data, fits.ImageHDU) or isinstance(data, fits.HDUList):
        return True
    else:
        return False


def get_fits(data, form, ext=None, verbose=False):
    """
    Return FITS data to the assigned form.

    Parameters:
        data : file path, object, or ndarray
            Name or location of data to get in the assigned form.
        form : {'hdulist', 'hdu', 'data', 'header'}
            The assigned form of data.
        ext : int, optional
            Extension number in HDUlist with multiple HDU.
        verbose : bool
            Print information.

    Returns:
        hdulist, hdu, data, or header
    """
    if form not in ['hdulist', 'hdu', 'data', 'header']:
        raise TypeError("Possible forms: 'hdulist', 'hdu', 'data', or 'header'.")
    if form is 'data' and isinstance(data, np.ndarray):
        return data*1.
    elif form is 'hdu' and (isinstance(data, fits.PrimaryHDU) or isinstance(data, fits.ImageHDU)):
        return fits.PrimaryHDU(data.data, data.header)
    elif form is 'hdulist' and isinstance(data, fits.HDUList):
        return fits.HDUList([data[i] for i in range(len(data))])
    if isinstance(data, str):
        hdulist = fits.open(data)
        df = 'fits.hdulist'
    elif isinstance(data, fits.HDUList):
        hdulist = data
        df = 'fits.hdulist'
    else:
        hdulist = None
    if hdulist is None:
        if isinstance(data, fits.PrimaryHDU) or isinstance(data, fits.ImageHDU):
            hdu = data
            df = 'fits.image'
        else:
            hdu = None
    else:
        if verbose:
            hdulist.info()
        if ext is None:
            ext = 0
            if len(hdulist) > 1:
                while hdulist[ext].header['NAXIS'] < 2:
                    ext += 1
                warn("'ext' is not assigned. Extension {} is opened.".format(ext))
        hdu = hdulist[ext]
    if hdu is None:
        if isinstance(data, np.ndarray):
            df = 'numpy.ndarray'
        else:
            data = None
    else:
        try:
            remove_history = hdu.header['HISTORY']
            for i in range(len(remove_history)):
                hdu.header.remove('HISTORY')
        except:
            pass
        data = hdu.data
    if form == 'hdulist' and hdulist is not None:
        return hdulist
    elif form == 'hdu' and hdu is not None:
        return hdu
    elif form == 'data' and data is not None:
        return data
    elif form == 'header' and hdu is not None:
        return hdu.header
    else:
        raise IOError("'{}' does not contain '{}'".format(df, form))


def check_all(file, naxis=3, ext=None, verbose=False):
    return check_velo(check_axis(file, naxis, ext, verbose))


def check_axis(file, naxis, ext=None, verbose=False):
    """
    Check FITS structure and return HDU data with the structure and header you want.

    Parameters:
        file : FITS file path or object (HDUList or HDU)
            Name or location of data to check and open in assigned structure.
        naxis : int
            Number of axes (dimension of shape) that you want.
        ext : int
            Extension number in HDUlist with multiple HDU.
        verbose : bool
            If True, print information.

    Returns:
        HDU (image data and header)
    """
    hdu = get_fits(file, 'hdu', ext, verbose)
    if verbose:
        print('\nOriginal shape is {}'.format(hdu.data.shape), end='')
    if len(hdu.data.shape) == naxis and hdu.header['NAXIS'] == naxis:
        if verbose:
            print(' that is matched the required shape.')
        return hdu
    elif len(hdu.data.shape) < naxis or hdu.header['NAXIS'] < naxis:
        if verbose:
            print(' that is NOT matched the required shape.')
        raise TypeError('FITS data does not have as many axes as required')
    else:
        for ra in range(naxis+1, 5):
            remove_fields = hdu.header['*'+str(ra)]
            remove_fields += hdu.header['PC'+str(ra)+'*']
            for rf in remove_fields:
                try:
                    hdu.header.remove(rf)
                except:
                    pass
        while len(hdu.data.shape) > naxis:
            hdu.data = hdu.data[0]
        if verbose:
            print(".\nReshaped to {}.".format(str(hdu.data.shape)))
        return hdu


def check_velo(hdu, ext=None):
    """
    Check 3rd axis of cube FITS,
    if type of 3rd axis is not 'VRAD' (radial velocity), reform to 'VRAD',
    and return HDU data and header.

    Parameters:
        hdu : file path, object (HDUList or HDU)
             Name or location of cube data to check and reform to 'VRAD' type
        ext : int
            Extension number in HDUlist with multiple HDU.

    Returns:
        HDU (image data and header)
    """
    hdu = get_fits(hdu, 'hdu', ext)
    h = hdu.header
    if not h['CTYPE3'] in ['VRAD', 'FREQ']:
        raise TypeError('{} is not supported type.'.format(h['CTYPE3']))
    ch = (np.arange(h['NAXIS3'])+1.-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
    cp = h['CRPIX3']
    if round(cp, 0) == round(cp, 5):
        cp = int(cp)
    else:
        cp = int(cp)
        h['CRPIX3'] = cp
        h['CRVAL3'] = ch[cp]
        ch = (np.arange(h['NAXIS3'])+1.-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
    if h['CTYPE3'] == 'FREQ':
        rf = h['CRVAL3']
        try:
            fu = u.Unit(h['CUNIT3'])
        except:
            fu = u.Hz
        ch = (ch*fu).to(u.m, equivalencies=u.spectral())
        rf = (rf*fu).to(u.m, equivalencies=u.spectral())
        ch = ((ch-rf)/rf*co.c).to(u.km/u.s)
        h['CTYPE3'] = 'VRAD'
    else:
        try:
            ch = ch*u.Unit(h['CUNIT3'])
        except:
            ch = ch*u.m/u.s
        ch = ch.to(u.km/u.s)
    if ch[1]-ch[0] < 0:
        ch = np.flip(ch, axis=0)
        hdu.data = np.flip(hdu.data, axis=0)
        cp = h['NAXIS3']-cp+1
    h['CRVAL3'] = ch[cp].value
    h['CDELT3'] = (ch[cp]-ch[cp-1]).value
    h['CRPIX3'] = cp
    h['CUNIT3'] = 'km/s'
    return fits.PrimaryHDU(hdu.data, h)


def make_velo(hdu, ext=None):
    """
    Make and return velocity array of cube fits.

    Parameters:
        hdu : file path, object (HDUList or HDU), or header
             Name or location of cube data to get header information.
        ext : int
            Extension number in HDUlist with multiple HDU.

    Returns:
        array converted in km/s unit.
    """
    if isinstance(hdu, fits.Header):
        h = hdu
    else:
        h = get_fits(hdu, 'header', ext)
    if not h['CTYPE3'] == 'VRAD':
        raise TypeError('{} is not supported type.'.format(h['CYTPE3']))
    ch = (np.arange(h['NAXIS3'])+1.-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
    try:
        ch = ch*u.Unit(h['CUNIT3'])
    except:
        ch = ch*u.m/u.s
    ch = ch.to(u.km/u.s).value.round(5)
    return ch
