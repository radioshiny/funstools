from .cube2map import Cube2map
import os
from glob import glob
from matplotlib import pyplot as plt
import numpy as np


def full_line_scan(loc='', vr=None, yi=1.5, cut=None, ver='otfpro'):
    """
    Interactive line scan for full data set of FUNS project.
    This recognizes only the 'fullcube' fits files in stored in the 'loc' directory.
    The C18O file is required and the remaining lines are optional.

    Parameters:
        loc : path
            Directory path stored fullcube fits files.
        vr : list or tuple, optional
            Velocity (x-axis) range of the plot for line profile.
        yi : float
            Temperature (y-axis) interval to plot multi-lines.
        cut : list [ra_start, dec_start, ra_size, dec_size], optional
            Position and size of map cut in axis fraction.

    Returns:
        matplotlib.figure, matplotlib.axes, matplotlib.axes
            Tuple with 3 components.
            You can modify the detail elements to save and print the figure.
    """
    if ver == 'otfpro':
        release = False
    elif ver == 'release':
        release = True
    else:
        raise ValueError("ver = {}, full_line_scan only work for 'otfpro' or 'release' version.")
    if not len(loc) == 0:
        if not loc[-1] == '/':
            loc = loc+'/'
    if release:
        base = glob(loc+'*C18O*_match_cube.fits')
    else:
        base = glob(loc+'*C18O*_fullcube.fits')
    if len(base) == 0:
        raise IOError('{}: No such file.'.format(base))
    else:
        args = (base[0].split('/')[-1]).split('_')
    files = []
    # lines = ['13CO', 'C18O', 'CS', 'HCOP', 'N2HP', 'SO', 'NH2D', 'H13COP']
    dcw = [1, 0, 1, 0, 1, 1, 0, 0]
    acw = [0, 1, 0, 1, 0, 0, 1, 1]
    ors = [200, 330]
    lines = ['H13COP', 'NH2D', 'N2HP', 'SO', 'HCOP', 'CS', 'C18O', '13CO']
    tlin = [r'$\mathrm{H^{13}CO^+}$', r'$\mathrm{NH_2D}$', r'$\mathrm{N_2H^+}$', 'SO', r'$\mathrm{HCO^+}$', 'CS', r'$\mathrm{C^{18}O}$', r'$\mathrm{^{13}CO}$']
    lnam = []
    tnam = []
    rs = []
    for ln, lin in enumerate(lines):
        if release:
            nam = [loc+args[0]+'_'+lin+'_v10_match_cube.fits', loc+args[0]+'_'+lin+'_match_cube.fits']
        else:
            nam = [loc+args[0]+'_'+lin+'_all_'+args[3]+'_v10_fullcube.fits', loc+args[0]+'_'+lin+'_all_'+args[3]+'_fullcube.fits']
        if os.path.exists(nam[dcw[ln]]):
            files.append(nam[dcw[ln]])
            rs.append(ors[dcw[ln]])
            lnam.append(lin)
            tnam.append(tlin[ln])
        elif os.path.exists(nam[acw[ln]]):
            files.append(nam[acw[ln]])
            rs.append(ors[acw[ln]])
            lnam.append(lin)
            tnam.append(tlin[ln])
        else:
            print('{}\n{}\n: No such files.'.format(nam[0], nam[1]))
    if 'C18O' in lnam:
        li = lnam.index('C18O')
        l0 = Cube2map(files[li], rmssize=rs[li], velocity_smo=3, spatial_smo=2)
    else:
        l0 = Cube2map(files[-1], rmssize=rs[-1], velocity_smo=3, spatial_smo=2)
    lset = []
    for ln, lin in enumerate(files):
        lset.append(Cube2map(lin, rmssize=rs[ln], smoothing=False, masking=False))

    plt.ion()
    x = l0.x
    if vr is None:
        vr = [x[rs[-1]], x[-rs[-1]]]
    elif isinstance(vr, list) or isinstance(vr, tuple):
        pass
    else:
        raise TypeError("'vr' is not list or tuple.".format(vr))
    if cut is None:
        cut = [0., 0., 1., 1.]
    elif isinstance(cut, list) and len(cut) == 4:
        pass
    else:
        raise TypeError("'cut' is not readable form of [xs, ys, xl, yl].")
    # cut data
    rcut = [round(l0.nr*cut[0]), round(l0.nr*(cut[0]+cut[2]))]
    dcut = [round(l0.nd*cut[1]), round(l0.nd*(cut[1]+cut[3]))]
    cdat = l0.data[:, dcut[0]:dcut[1], rcut[0]:rcut[1]]
    crms = l0.rms[dcut[0]:dcut[1], rcut[0]:rcut[1]]
    # make map
    cdet = np.full(cdat.shape, 0.)
    cdet[np.nan_to_num(cdat) > 3*np.nanmedian(crms)] = 1.
    csum = np.nansum(cdet, axis=(1, 2))
    imap = np.sum(cdat[csum > 2*np.mean(csum)], axis=0)
    # map drawing
    lsfig = plt.figure(figsize=(16, 9))
    # lsfig.subplots_adjust(0.1, 0.1, 0.95, 0.95, 0.2, 0.2)
    lsmap = lsfig.add_axes([0.06, 0.07, 0.55, 0.9], projection=l0.wcs2d)
    lsmap.set_xlabel('R.A.')
    lsmap.set_ylabel('Dec.')
    lsmap.imshow(imap, origin='lower', vmin=0, cmap='inferno')
    cp = lsmap.scatter(-100, -100, s=60, color='lime')
    lsmap.set_xlim(rcut[0]-0.5, rcut[1]-0.5)
    lsmap.set_ylim(dcut[0]-0.5, dcut[1]-0.5)
    # line plot
    # lslin = lsfig.add_subplot(122)
    lslin = lsfig.add_axes([0.67, 0.07, 0.28, 0.9])
    if lnam[-1] == '13CO':
        ymax = 0.5*np.nanmax(lset[-1].data)
    else:
        ymax = max([np.nanmax(i.data) for i in lset])
    yr = [-4*np.nanmedian(crms), yi*(len(lset)-1)+ymax+np.nanmedian(crms)]
    print(ymax, yr)
    wset = [i.wcs2d for i in lset]

    def _onclick(event):
        if event.inaxes != lsmap.axes: return
        ri = int(round(event.xdata))
        di = int(round(event.ydata))
        if not release:
            rd = l0.wcs2d.wcs_pix2world(ri, di, 0)
            px = [i.wcs_world2pix(rd[0], rd[1], 0) for i in wset]
        lslin.clear()
        lslin.set_xlabel(r'$V_\mathrm{LSR}$ [km/s]')
        lslin.set_ylabel(r'$T_\mathrm{A}^* [K]$')
        lslin.set_xlim(*vr)
        lslin.set_ylim(*yr)
        lslin.set_yticks(np.arange(0, yr[1], yi))
        for li, ld in enumerate(lset):
            if release:
                r, d = ri*1, di*1
            else:
                r = int(round(px[li][0].item()))
                d = int(round(px[li][1].item()))
            lslin.plot([x[0], x[-1]], [li*yi, li*yi], lw=1, alpha=0.5, color='black')
            # lslin.plot([x[0], x[-1]], [3*l0.srms[d, r], 3*l0.srms[d, r]], ls='dotted', lw=1, alpha=0.3, color='blue')
            # lslin.plot([x[0], x[-1]], [3*ld.rms[d, r]+li*yi, 3*ld.rms[d, r]+li*yi], ls='dotted', lw=1, alpha=0.5, color='blue')
            # lslin.step(x, l0.sdata[:, d, r], color='blue', lw=1, alpha=0.6)
            if (0 <= r < ld.nr) and (0 <= d < ld.nd):
                cp.set_offsets([r, d])
                if lnam[li] == 'N2HP':
                    lslin.step(ld.x+8., ld.data[:, d, r]+li*yi, color='red', lw=1)
                elif lnam[li] == '13CO':
                    lslin.step(ld.x, 0.5*ld.data[:, d, r]+li*yi, color='red', lw=1)
                    lslin.annotate(r'$(\times 0.5)$', (0.97, ((li+0.3)*yi-yr[0])/(yr[1]-yr[0])), xycoords='axes fraction', ha='right')
                else:
                    lslin.step(ld.x, ld.data[:, d, r]+li*yi, color='red', lw=1)
            else:
                pass
            lslin.annotate('[{}, {}]'.format(d, r), (0.03, ((li+0.3)*yi-yr[0])/(yr[1]-yr[0])), xycoords='axes fraction')
            lslin.annotate(tnam[li], (1.01, (li*yi-yr[0])/(yr[1]-yr[0])), xycoords='axes fraction')
        c = l0.wcs2d.pixel_to_world(ri, di)
        lslin.annotate('Coord.: {}'.format(c.to_string('hmsdms')), (0.03, 0.95), xycoords='axes fraction')
        lslin.figure.canvas.draw()

    lsfig.canvas.mpl_connect('button_press_event', _onclick)
    return lsfig, lsmap, lslin
