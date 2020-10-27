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
rms = cube.rms
rms_smoothed = cube.srms
```
You can save maps you made to fits using `save_fits`.  
`Cube2map.header2d` is FITS header for 2d-image data.
```python
from funstools import save_fits
save_fits('rms.fits', cube.rms, cube.header2d, overwrite=True)
save_fits('rms_smoothed.fits', cube.srms, cube.header2d, overwrite=True)
```
`Cube2map.wcs2d` is 2d (RA-Dec) wcs information for matplotlib projection.
```python
from matplotlib import pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection=cube.wcs2d)
ax.imshow(cube.rms, origin='lower')
fig.savefig('rms.pdf')
```

### Moments maps

#### Moment 0 (integrated intensity)

`Cube2map.moment0()` and `Cube2map.m0` return the integrated intensity map of input cube data.

![equation](https://latex.codecogs.com/svg.latex?M_0=\Delta&space;v\sum&space;I_i)

```python
# Using velocity range
cube.moment0(vr=[v1, v2])

# Using channel range
cube.moment0(cr=[ch1, ch2])

# Return the most recently computed moment 0 map
# If there is not pre-computed moment 0 map,
#    return the moment 0 map with cr=[rmssize, -rmssize]
cube.m0

# Return with RMS_moment0
m0, m0rms = cube.moment0(verbose=True)

# Save to fits
save_fits('moment0.fits', cube.m0, cube.header2d, overwrite=True)

# Make figure
fig = plt.figure()
ax = fig.add_subplot(111, projection=cube.wcs2d)
ax.imshow(cube.m0, origin='lower', vmin=0, vmax=np.nanpercentile(cube.m0, 99.9))
fig.savefig('moment0.pdf')
```

#### Moment 1 (intensity weighted mean velocity)

`Cube2map.moment1()` and `Cube2map.m1` return the intensity weighted mean velocity map of input cube data.

![equation](https://latex.codecogs.com/svg.latex?M_1=\frac{\sum&space;I_i&space;v_i}{M_0})

```python
# Using channel range
cube.moment1(cr=[ch1, ch2])

# Return the most recently computed moment 1 map
cube.m1
```

#### Moment 2 (velocity dispersion)

`Cube2map.moment2()` and `Cube2map.m2` return the intensity weighted mean velocity map of input cube data.

![equation](https://latex.codecogs.com/svg.latex?M_2=\sqrt{\frac{\sum&space;I_i\left(v_i-M_1\right)^2}{M_0}})

```python
# Using channel range
cube.moment2(cr=[ch1, ch2])

# Return the most recently computed moment 2 map
cube.m2
```

#### Example script for `Cube2map.moment#()`

```python
from funstools import Cube2map
from matplotlib import pyplot as plt
from multiaxes import Multiaxes

# set baseline channel size
rs = 200

# load cube data
cube = Cube2map('/.../TRAO-FUNS/SerB/release/SerB_C18O_v10_match_cube.fits', 
                rmssize=rs, velocity_smo=2, spatial_smo=2)

# set map x-y ratio
xyr = float(cube.nr/cube.nd)

# draw figure using Multiaxes
mx = Multiaxes(2, 3, 1, xyr, 0.36, 0.56, cb=0.1, clab=0.34, scale=0.7, proj=cube.wcs2d)
mx.shareaxes((False, True), 0.1)
fig, ax, cax = mx.drawfig()

# set map title, colormap, vmin, vmax
title = [r'Integrated inteisity (K km s$^{-1}$)',
         r'Intensity weighted mean velocity (km s$^{-1}$)',
         r'Velocity dispersion (km s$^{-1}$)']
cmaps = ['inferno', 'jet', 'coolwarm']
vmin = [0., 7.3, 0.3]
vmax = [4.5, 9.8, 1.3]

# make moment maps and draw figure
map = []
for i in range(3):
    map.append(getattr(cube, 'moment{}'.format(i))())
    cmap = plt.get_cmap(cmaps[i])
    cmap.set_bad('grey')
    cmaps[i] = cmap
    cs = ax[i].imshow(map[i], cmap=cmaps[i], vmin=vmin[i], vmax=vmax[i])
    plt.colorbar(cs, cax[i], orientation='horizontal')
    cax[i].set_xlabel(title[i])
    mx.topcolorbar(cax[i])
    if i == 0:
        ax[i].coords[0].set_axislabel('R.A. (J2000)')
        ax[i].coords[1].set_axislabel('Dec. (J2000)')
    else:
        ax[i].coords[0].set_axislabel(' ')
        ax[i].coords[1].set_axislabel(' ')

# save figure
fig.savefig('images/SerB_C18O_moment_maps.png')
```

![example](./images/SerB_C18O_moment_maps.png) 

### Channel maps

`Cube2map.chmap()` return the channel maps and labels for their velocity ranges.

```python
# Using velocity range
maps, labels = cube.chmap(mn=9, vr=(vr1, vr2))

# Using channel range
maps, labels = cube.chmap(mn=9, cr=(cr1, cr2))
```

#### Example for `Cube2map.chmap()`

```python
from funstools import Cube2map
from matplotlib import pyplot as plt
from multiaxes import Multiaxes

# set baseline channel size
rs = 200

# load cube data
cube = Cube2map('/.../TRAO-FUNS/SerB/release/SerB_C18O_v10_match_cube.fits', 
                rmssize=rs, velocity_smo=2, spatial_smo=2)

# draw figure using Multiaxes
mx = Multiaxes(2, 3, 2, 1, 0.36, 0.56, margin=(0.02, 0.02, 0.71, 0.02), scale=0.7, proj=cube.wcs2d)
mx.shareaxes((True, True), 0.1)
fig, ax, _ = mx.drawfig()
cax = mx.sharecolorbar('right', 0.15, 0.1)

# set map value_range, colormap
vr = [0., 1.55]
cmap = plt.get_cmap('nipy_spectral')
cmap.set_bad('grey')

# make moment maps and draw figure
maps, labels = cube.chmap(6, (6.3, 9.9))
for i in range(6):
    iax = ax[i//3, i%3]
    cs = iax.imshow(maps[i], cmap=cmap, vmin=vr[0], vmax=vr[1])
    iax.set_xlim(26.5, 66.5)
    iax.set_ylim(163.5, 203.5)
    iax.annotate(labels[i]+' km/s', xy=(0.55, 0.05), xycoords='axes fraction', 
                 ha='center', va='baseline', color='white', fontsize='large')
    if i == 3:
        iax.coords[0].set_axislabel('R.A. (J2000)')
        iax.coords[1].set_axislabel('Dec. (J2000)')
    else:
        iax.coords[0].set_axislabel(' ')
        iax.coords[1].set_axislabel(' ')
plt.colorbar(cs, cax)
cax.set_ylabel(r'Integrated intensity (K km s$^{-1}$)')

# save figure
fig.savefig('images/SerB_C18O_channel_maps.png')
```

![example](./images/SerB_C18O_channel_maps.png) 

### Line scanning

`Cube2map.line_scan()` is the viewer of line profiles for 3 by 3 pixels.

Since the `onclick` event is activated, you can move a location by clicking a
pixel on the left intensity map or clicking the direction you want to go on the
right line profile map.

```python
# set the velocity range (x-axis range of plot)
cube.line_scan(vr=(5.5, 11.))
```

![example](./images/SerB_C18O_line_scan.png)

### Full line set scanning

`full_line_scan()` is the viewer of line profiles for full line sets of TRAO-FUNS.

```python
from funstools import full_line_scan

loc = '/.../TRAO-FUNS/SerB/release/'
full_line_scan(loc, vr=)

``` 

### Gaussian decomposing



### Finding velocity-coherent structure

### Calculating physical properties