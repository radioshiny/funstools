import numpy as np
from astropy import units as u


_freq = {'13CO': 110.2015430*u.GHz,
         'C18O': 109.7821734*u.GHz,
         'N2HP': 93.1737637*u.GHz,
         'HCOP': 89.1885247*u.GHz,
         'CS': 97.9809533*u.GHz,
         'SO': 99.2998700*u.GHz,
         'NH2D': 85.9262780*u.GHz,
         'H13COP': 86.7542884*u.GHz,
         70: (70*u.um).to(u.GHz, equivalencies=u.spectral()),
         100: (100*u.um).to(u.GHz, equivalencies=u.spectral()),
         160: (160*u.um).to(u.GHz, equivalencies=u.spectral()),
         250: (250*u.um).to(u.GHz, equivalencies=u.spectral()),
         350: (350*u.um).to(u.GHz, equivalencies=u.spectral()),
         500: (500*u.um).to(u.GHz, equivalencies=u.spectral())}

_wlen = dict(zip(_freq.keys(), [_freq[i].to(u.mm, equivalencies=u.spectral()) for i in _freq.keys()]))


_hgbs_beam = {70: np.sqrt(8.80*9.60)*u.arcsec,
              100: np.sqrt(6.66*6.89)*u.arcsec,
              160: np.sqrt(10.55*12.08)*u.arcsec,
              250: np.sqrt(18.9*18.0)*u.arcsec,
              350: np.sqrt(25.8*24.6)*u.arcsec,
              500: np.sqrt(38.3*35.2)*u.arcsec}


def trao_beam(line):
    a = 0.03333836
    b = 6.0357315
    c = 0.64283393
    w = _wlen[line].to(u.m).value
    return (a*(w/b)**c*u.rad).to(u.arcsec).round(1)


def trao_aeff(line):
    a = 45.43378521
    b = 2.31232982
    c = 10.23642098
    w = _wlen[line].to(u.mm).value
    return ((a*np.exp(-(b/w)**c))/100).round(3)


def trao_beff(line):
    a = 52.40088923
    b = 0.14793845
    c = -14.82565966
    f = _freq[line].to(u.GHz).value
    return ((a-np.exp(b*f+c))/100).round(3)


def hgbs_beam(wl):
    if wl in [70, 100, 160, 250, 350, 500]:
        return _hgbs_beam[wl]
    else:
        raise ValueError('Undefined wave length.')
