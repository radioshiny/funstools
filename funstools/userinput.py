from astropy.io import fits
from .checkshape import make_velo
from warnings import warn
import numpy as np


def yesno(question, default=None):
    valid = {'y':True, 'ye':True, 'yes':True, 'n':False, 'no': False}
    if default is None:
        prompt = ' [y/n]'
    elif default == 'yes':
        prompt = ' [Y/n], [Enter] is Yes.'
    elif default == 'no':
        prompt = ' [y/N], [Enter] is No.'
    else:
        raise ValueError('Invalid default answer.')
    while True:
        answer = input(question+prompt).lower()
        if default is not None and answer == '':
            return valid[default]
        elif answer in valid:
            return valid[answer]
        else:
            print("Please respond with 'yes' or 'no' ('y' or 'n').")


def v2ch(header, vr):
    """
    Convert velocity range to channel range.

    Parameters:
        vr : list or tuple
            Velocity range, [start, end] or (start, end).

    Returns:
        cr : list
            Channel range, [start, end+1] ('+1' follows the python indexing rule).
    """
    if vr is None:
        raise ValueError("'vr', velocity range is required.")
    if not (isinstance(vr, list) or isinstance(vr, tuple)):
        raise TypeError('{} is not list or tuple.'.format(vr))
    if isinstance(header, fits.Header):
        x = make_velo(header)
    else:
        raise TypeError('{} is not fits.Header.'.format(header))

    if x[0] > vr[0] or x[-1] < vr[1]:
        warn('Given velocity range exceeds the velocity axis of cube data.')
    return [np.argmin(np.abs(x-float(vr[0]))), np.argmin(np.abs(x-float(vr[1])))+1]
