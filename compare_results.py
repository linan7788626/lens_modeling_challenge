import sie_lensing_simulation as sls
import numpy as np
import pylab as pl


def single_run_test(ind, ysc1, ysc2, q, vd, pha, zl, zs):
    nnn = 300  # Image dimension
    bsz = 9.0  # arcsecs
    dsx = bsz / nnn         # pixel size of SDSS detector.

    xi1, xi2 = sls.make_r_coor(nnn, dsx)

    xc1 = 0.0
    xc2 = 0.0
    rc = 0.0  # Core size of lens (in units of Einstein radius).
    re = sls.re_sv(vd, zl, zs)  # Einstein radius of lens.
    lpar = np.asarray([xc1, xc2, q, rc, re, pha])
    ai1, ai2, kappa_out, shear1, shear2, mua = sls.lensing_signals_sie(
        xi1, xi2, lpar)

    # a11, a12 = np.gradient(ai1, dsx)
    # a21, a22 = np.gradient(ai2, dsx)
    # kappa_out = 0.5 * (a11 + a22)

    yi1 = xi1 - ai1
    yi2 = xi2 - ai2

    return xi1, xi2, yi1, yi2, kappa_out, mua

if __name__ == '__main__':
    num_imgs = 1
    sourcpos = 0.0

    ysc1 = [0.4]
    ysc2 = [-0.3]
    zl = 0.298  # zl is the redshift of the lens galaxy.
    zs = 1.0
    vd = [320]  # Velocity Dispersion.
    q = [0.5]
    pha = [-45.0]

    i = 0
    xf1, xf2, yf1, yf2, kappa_in, mua_in = single_run_test(
        i, ysc1[i], ysc2[i], q[i], vd[i], pha[i], zl, zs)
    kappa_in = kappa_in

    levels = [-1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0]

    kappa_out = np.loadtxt(
        "./fits_outputs/data_from_Quinn/bestfit_q0.5.logkappa")

    pl.figure(figsize=(8, 8))
    pl.contour(xf1, xf2, np.log10(kappa_in), levels, colors=('k',))
    pl.contour(xf1, xf2, kappa_out, levels, colors=('r',))

    mua_out = np.loadtxt("./fits_outputs/data_from_Quinn/bestfit_q0.5.maglog")

    pl.figure(figsize=(8, 8))
    pl.contour(xf1, xf2, np.log10(np.abs(mua_in)), colors=('k',))
    pl.contour(xf1, xf2, mua_out, colors=('r',))

    pl.figure(figsize=(8, 8))
    pl.contour(yf1, yf2, np.log10(np.abs(mua_in)), colors=('k',))
    pl.contour(yf1, yf2, 10.0**mua_out, colors=('r',))

    vpix = np.loadtxt("./fits_outputs/data_from_Quinn/src_pixel.dat")

    pl.show()
