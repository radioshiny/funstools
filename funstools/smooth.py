from astropy.convolution import convolve
from astropy.convolution.kernels import Gaussian1DKernel, Gaussian2DKernel, Tophat2DKernel, Box1DKernel
from astropy import units as u
import numpy as np


def smooth1d(data, fwhm):
    """
    Convolve data with 1d Gaussian kernel and return.
    For cube data, the convolution applies only to velocity axis.

    Parameters:
        data : ndarray (1d or 3d array)
        fwhm : float
            FWHM size of 1d Gaussian kernel.

    Returns:
        ndarray (smoothed data)
    """
    if float(fwhm) == 0.:
        return data
    kernel = Gaussian1DKernel(fwhm/np.sqrt(8.*np.log(2.)))
    kl = len(kernel.array)
    data = np.r_[data[kl-1:0:-1], data, data[-2:-kl-1:-1]]
    if len(data.shape) == 1:
        return convolve(data, kernel, boundary=None)[kl-1:-kl+1]
    elif len(data.shape) == 2:
        raise TypeError('1D-smoothing can not be applied to 2D-array.')
    elif len(data.shape) == 3:
        smodata = np.full(data.shape, np.nan)
        for i in range(data.shape[1]):
            for j in range(data.shape[2]):
                smodata[:, i, j] = convolve(data[:, i, j], kernel, boundary=None)
        return smodata[kl-1:-kl+1]


def smooth2d(data, fwhm, pa=None):
    """
    Convolve data with 2d Gaussian kernel and return.
    For cube data, the convolution applies only to RA-Dec plane.

    Parameters:
        data : ndarray (2d or 3d array)
        fwhm : float, list or tuple
            FWHM size of 2d Gaussian kernel in x and y before rotating.
            if list or tuple, fwhm[0] = x_size and fwhm[1] = y_size.
        pa : float
            Position angle (counterclockwise) in degree.

    Returns:
        ndarray (smoothed data)
    """
    if float(fwhm) == 0.:
        return data
    if isinstance(fwhm, int) or isinstance(fwhm, float):
        kernel = Gaussian2DKernel(fwhm/np.sqrt(8.*np.log(2.)))
    elif len(fwhm) == 2:
        if pa is None:
            kernel = Gaussian2DKernel(fwhm[0]/np.sqrt(8.*np.log(2.)), fwhm[1]/np.sqrt(8.*np.log(2.)))
        else:
            theta = pa*u.deg.to(u.rad)
            kernel = Gaussian2DKernel(fwhm[0]/np.sqrt(8.*np.log(2.)), fwhm[1]/np.sqrt(8.*np.log(2.)), theta)
    else:
        raise ValueError('2D-smoothing size: float (fwhm) or tuple (x_fwhm, y_fwhm)')
    if len(data.shape) == 1:
        raise TypeError('2D-smoothing can not be applied to 1D-array.')
    elif len(data.shape) == 2:
        return convolve(data, kernel, boundary='extend', preserve_nan=True)
    elif len(data.shape) == 3:
        smodata = np.full(data.shape, np.nan)
        for i in range(len(data)):
            smodata[i] = convolve(data[i], kernel, boundary='extend', preserve_nan=True)
        return smodata


def smooth3d(data, fwhm_spatial, fwhm_velocity, pa=None):
    """
    Convolve cube data with 1d and 2d Gaussian kernel and return.

    Parameters:
        data : ndarray (3d array)
        fwhm_spatial : float, list or tuple
            FWHM size of 2d Gaussian kernel in x and y before rotating.
            if list or tuple, fwhm[0] = x_size and fwhm[1] = y_size.
        fwhm_velocity : float
            FWHM size of 1d Gaussian kernel.
        pa : float
            Position angle (counterclockwise) in degree.

    Returns:
        ndarray (smoothed data)
    """
    if len(data.shape) == 3:
        return smooth1d(smooth2d(data, fwhm_spatial, pa), fwhm_velocity)
    else:
        raise TypeError('3D-smoothing can not be applied to 1D or 2D-array.')


def radsmo2d(data, radius):
    """
    Convolve data with 2d Tophat kernel and return.
    Tophat is an isotropic radial smoothing filter.
    For cube data, the convolution applies only to RA-Dec plane.

    Parameters:
        data : ndarray (2d or 3d array)
        radius : int
            Radius of Tophat kernel.

    Returns:
        ndarray (smoothed data)
    """
    if not isinstance(radius, int):
        raise TypeError("'radius' should be given as int.")
    if int(radius) == 0:
        return data
    kernel = Tophat2DKernel(radius)
    if isinstance(data, np.ndarray):
        if len(data.shape) == 2:
            return convolve(data, kernel, boundary='extend', preserve_nan=True)
        elif len(data.shape) == 3:
            smodata = np.full_like(data, np.nan)
            for i in range(len(data)):
                smodata[i] = convolve(data[i], kernel, boundary='extend', preserve_nan=True)
            return smodata
        else:
            raise TypeError("'data' should be given as 2D-image or 3D-cube.")
    else:
        raise TypeError("'data' should be given as numpy.ndarray")


def boxsmo1d(data, width):
    """
    Convolve data with 1d Box kernel and return.
    For cube data, the convolution applies only to RA-Dec plane.

    Parameters:
        data : ndarray (2d or 3d array)
        width : float
            Width of Box kernel.

    Returns:
        ndarray (smoothed data)
    """
    if not (isinstance(width, int) or isinstance(width, float)):
        raise TypeError("'width' should be given as int or float.")
    if float(width) == 0.:
        return data
    kernel = Box1DKernel(width)
    if isinstance(data, np.ndarray):
        if len(data.shape) == 1:
            return convolve(data, kernel, boundary='extend', preserve_nan=True)
        elif len(data.shape) == 3:
            smodata = np.full_like(data, np.nan)
            for i in range(data.shape[1]):
                for j in range(data.shape[2]):
                    smodata[:, i, j] = convolve(data[:, i, j], kernel, boundary='extend', preserve_nan=True)
            return smodata
        else:
            raise TypeError("'data' should be given as 1D-spectrum or 3D-cube.")
    else:
        raise TypeError("'data' should be given as numpy.ndarray")
