# FUNStools

`funstools` is python toolkit for FUNS project.

`funstools` is a toolkit that contains the useful functions and core tools 
needed for performing FUNS project with python. It is based on the `numpy` and 
`astropy` packages, which have been validated by astronomers for a long time.

## Install

```bash
pip install git+https://github.com/radioshiny/funstools.git
```

## Requirement

`astropy`
`matplotlib`
`scipy`

## User Manual

### RMS noise map

#### Using `get_rms`
`get_rms` returns the RMS noise map of input cube data.
```python
from funstools import get_rms
rms = get_rms(cube, where='both', size=200)
```

#### Using `Cube2map`
`Cube2map` is more useful to make maps including RMS noise map from cube data.
```python
from funstools import Cube2map
a = Cube2map('sample_cube.fits', getrms='both', rmssize=200, velocity_smo=2, spatial_smo=2)
rms = a.rms
rms_smoothed = a.srms
```
You can save maps you made to fits using `save_fits`.  
`Cube2map.header2d` is FITS header for 2d-image data.
```python
from funstools import save_fits
save_fits('rms.fits', a.rms, a.header2d, overwrite=True)
save_fits('rms_smoothed.fits', a.srms, a.header2d, overwrite=True)
```
`Cube2map.wcs2d` is 2d (RA-Dec) wcs information for matplotlib projection.
```python
from matplotlib import pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection=a.wcs2d)
ax.imshow(a.rms, origin='lower')
fig.savefig('rms.pdf')
```

### Moments maps
`Cube2map.moment0()` and `Cube2map.m0` return the integrated intensity map of 
input cube data.
```python
# Using velocity range
a.moment0(vr=[v1, v2])

# Using channel range
a.moment0(cr=[ch1, ch2])

# Automatic
a.m0

# Save to fits
save_fits('moment0.fits', a.m0, a.header2d, overwrite=True)

# Make figure
fig = plt.figure()
ax = fig.add_subplot(111, projection=a.wcs2d)
ax.imshow(a.m0, origin='lower', vmin=0, vmax=np.nanpercentile(a.m0, 99))
fig.savefig('moment0.pdf')
```

### Channel maps

### Line scanning

### Gaussian decomposing

### Finding velocity-coherent structure

### Calculating physical properties