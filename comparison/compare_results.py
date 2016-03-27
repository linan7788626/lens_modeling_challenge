#!/usr/bin/env python
import numpy as np
import libv4_cv as lv4
import mycosmology as mm
import scipy.signal as ss
import astropy.io.fits as pyfits
from astropy.cosmology import Planck13
import pylab as pl
import scipy.interpolate as sci
import pixcos2pixsdss as p2p
import congrid
import matplotlib.cm as cm


def rebin_psf(input_psf, new_shape):
    nxo, nyo = np.shape(input_psf)
    nxn, nyn = new_shape

    xo = np.linspace(0, nxo - 1.0, nxo) + 0.5
    yo = np.linspace(0, nyo - 1.0, nyo) + 0.5
    xo, yo = np.meshgrid(xo, yo)
    xo = xo.reshape((nxo * nyo))
    yo = yo.reshape((nxo * nyo))
    zo = input_psf.reshape((nxo * nyo))

    xn = np.linspace(0, nxo - 1.0, nxn) + 0.5
    yn = np.linspace(0, nyo - 1.0, nyn) + 0.5
    xn, yn = np.meshgrid(xn, yn)

    # print np.max(xo),np.min(xo)
    # print np.max(xn),np.min(xn)

    res = sci.griddata(np.array([xo, yo]).T, zo, (xn, yn), method='linear')
    return res


# def re0_sigma(sigma):
#    cv = 3e5
#    Dds = 1.0
#    Ds = 2.0
#    res = 4.0*np.pi*(sigma/cv)**2.0*Dds/Ds
#    return res

nMgyCount_r = 0.004760406   # nanomaggies per count for SDSS detector.
sky_r = 5.98          # SDSS typical r band sky
softbias = 1000.0        # SDSS softbias
Mgy2nanoMgy = 10e+9         # nanoMaggy to Maggy
aa_r = -24.149
kk = 0.156347
skycount = sky_r / (nMgyCount_r)
expsdss = 53.9
gain = 4.7
airmass = 1.201824
factor = 10.0**(0.4 * (aa_r + kk * airmass))


def psf_gaussian_norm(x1, x2, mu, sigma):
    r = np.sqrt(x1 * x1 + x2 * x2)
    res = 1.0 / (sigma * np.sqrt(2.0 * np.pi)) * \
        np.exp(-(r - mu)**2.0 / 2.0 * sigma**2.0)
    return res


def noise_map(nx1, nx2, nstd, NoiseType):
    if NoiseType == 'Poisson':
        noise = np.random.poisson(nstd, (nx1, nx2)) - nstd
    if NoiseType == 'Gaussian':
        noise = nstd * np.random.normal(0.0, 1.0, (nx1, nx2))
    return noise
#--------------------------------------------------------------------


def make_r_coor(nc, dsx):

    bsz = nc * dsx
    x1 = np.linspace(0, bsz - dsx, nc) - bsz / 2.0 + dsx / 2.0
    x2 = np.linspace(0, bsz - dsx, nc) - bsz / 2.0 + dsx / 2.0

    x2, x1 = np.meshgrid(x1, x2)
    return x1, x2


def make_c_coor(nc, dsx):

    bsz = nc * dsx
    x1, x2 = np.mgrid[0:(bsz - dsx):nc * 1j, 0:(bsz - dsx)
                         :nc * 1j] - bsz / 2.0 + dsx / 2.0
    return x1, x2

#--------------------------------------------------------------------


def lens_equation_sie(x1, x2, lpar):
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = lpar[0]
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]
    q = lpar[2]  # Ellipticity of lens.
    rc = lpar[3]  # Core size of lens (in units of Einstein radius).
    re = lpar[4]  # Einstein radius of lens.
    pha = lpar[5]  # Orintation of lens.

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1 - xc1) * cosa + (x2 - xc2) * sina
    xt2 = (x2 - xc2) * cosa - (x1 - xc1) * sina

    phi = np.sqrt(xt2 * xt2 + xt1 * q * xt1 * q + rc * rc)
    sq = np.sqrt(1.0 - q * q)
    pd1 = phi + rc / q
    pd2 = phi + rc * q
    fx1 = sq * xt1 / pd1
    fx2 = sq * xt2 / pd2
    qs = np.sqrt(q)

    a1 = qs / sq * np.arctan(fx1)
    a2 = qs / sq * np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 = -sina

    fx11 = xt11 / pd1 - xt1 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd1 * pd1)
    fx22 = xt22 / pd2 - xt2 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd2 * pd2)
    fx12 = xt12 / pd1 - xt1 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd1 * pd1)
    fx21 = xt21 / pd2 - xt2 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd2 * pd2)

    a11 = qs / (1.0 + fx1 * fx1) * fx11
    a22 = qs / (1.0 - fx2 * fx2) * fx22
    a12 = qs / (1.0 + fx1 * fx1) * fx12
    a21 = qs / (1.0 - fx2 * fx2) * fx21

    rea11 = (a11 * cosa - a21 * sina) * re
    rea22 = (a22 * cosa + a12 * sina) * re
    rea12 = (a12 * cosa - a22 * sina) * re
    rea21 = (a21 * cosa + a11 * sina) * re

    y11 = 1.0 - rea11
    y22 = 1.0 - rea22
    y12 = 0.0 - rea12
    y21 = 0.0 - rea21

    jacobian = y11 * y22 - y12 * y21
    mu = 1.0 / jacobian

    res1 = (a1 * cosa - a2 * sina) * re
    res2 = (a2 * cosa + a1 * sina) * re
    return res1, res2, mu
#--------------------------------------------------------------------


def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
    return (xnew, ynew)


def gauss_2d(x, y, par):
    (xnew, ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt(((xnew**2) * par[4] + (ynew**2) / par[4])) / np.abs(par[1])
    res = par[0] * np.exp(-res0**2.0)
    return res


def re_sv(sv, z1, z2):
    res = 4.0 * np.pi * (sv**2.0 / mm.vc**2.0) * \
        mm.Da2(z1, z2) / mm.Da(z2) * mm.apr
    return res

#----Fundamental Plain------------------------------


def Brightness(Re, Vd):
    a = 1.49
    b = 0.2
    c = -8.778
    mag_e = ((np.log10(Re) - a * np.log10(Vd) - c) / b) + \
        20.09  # Bernardi et al 2003
    nanoMgy = Mgy2nanoMgy * 10.0**(-(mag_e - 22.5) / 2.5)
    counts = nanoMgy / nMgyCount_r

    return counts


def de_vaucouleurs_2d(x, y, par):
    #[I0, Re, xc1,xc2,q,pha]
    # print "I0",par[0]
    # print "Re",par[1]
    (xnew, ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt((xnew**2) * par[4] + (ynew**2) / par[4]) / par[1]
    #res = par[0]*np.exp(-par[1]*res0**0.25)
    res = par[0] * np.exp(-7.669 * (res0**0.25 - 1.0))
    soften = par[0] * np.exp(-7.669 * ((0.2)**0.25 - 1.0))
    res[res > soften] = soften
    return res

# ----de Vaucouleurs profile-------------------------
# def deVaucouleurs(x,y,xc,yc,counts,R,e,phi):
    #theta   =phi*np.pi/180.

    #xx      =x-xc
    #yy      =y-yc
    #rx      =xx*np.cos(theta)+yy*np.sin(theta)
    #ry      =-xx*np.sin(theta)+yy*np.cos(theta)
    #rr      =np.sqrt(rx*rx/(1.0-e)+ry*ry*(1.0-e))
    #image   =counts*np.exp(-7.669*((rr/R)**0.25-1.0))
    #soften  =counts*np.exp(-7.669*((0.02)**0.25-1.0))
    #ix      =np.where(image>=soften)
    # image[ix]=soften

    # return image
#--------------------------------------------------------------------


def single_run_test(ind, ysc1, ysc2, q, vd, pha, zl, zs):
    dsx_sdss = 0.396         # pixel size of SDSS detector.
    R = 3.0000     #
    nnn = 300  # Image dimension
    bsz = 9.0  # arcsecs
    dsx = bsz / nnn         # pixel size of SDSS detector.
    nstd = 59  # ^2

    xx01 = np.linspace(-bsz / 2.0, bsz / 2.0 - dsx, nnn) + 0.5 * dsx
    xx02 = np.linspace(-bsz / 2.0, bsz / 2.0 - dsx, nnn) + 0.5 * dsx
    xi2, xi1 = np.meshgrid(xx01, xx02)
    #----------------------------------------------------------------------
    dsi = 0.03
    g_source = pyfits.getdata("./439.0_149.482739_1.889989_processed.fits")
    g_source = np.array(g_source, dtype="<d") * 10.0
    g_source[g_source <= 0.0001] = 1e-6
    #----------------------------------------------------------------------
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = 0.0
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0
    # q   = 0.7       #Ellipticity of lens.
    rc = 0.0  # Core size of lens (in units of Einstein radius).
    re = re_sv(vd, zl, zs)  # Einstein radius of lens.
    # pha = 45.0      #Orintation of lens.
    lpar = np.asarray([xc1, xc2, q, rc, re, pha])
    print lpar
    #----------------------------------------------------------------------
    ai1, ai2, mua = lens_equation_sie(xi1, xi2, lpar)

    a11, a12 = np.gradient(ai1, dsx)
    a21, a22 = np.gradient(ai2, dsx)

    kappa = 0.5 * (a11 + a22)

    yi1 = xi1 - ai1
    yi2 = xi2 - ai2

    g_limage = lv4.call_ray_tracing(g_source, yi1, yi2, ysc1, ysc2, dsi)
    g_limage[g_limage <= 0.0001] = 1e-6
    g_limage = p2p.cosccd2mag(g_limage)
    g_limage = p2p.mag2sdssccd(g_limage)

    # pl.figure()
    # pl.contourf(g_limage)
    # pl.colorbar()

    #-------------------------------------------------------------
    # Need to be Caliborate the mags
    dA = Planck13.comoving_distance(zl).value * 1000. / (1 + zl)
    Re = dA * np.sin(R * np.pi / 180. / 3600.)
    counts = Brightness(Re, vd)
    vpar = np.asarray([counts, R, xc1, xc2, q, pha])
    #g_lens = deVaucouleurs(xi1,xi2,xc1,xc2,counts,R,1.0-q,pha)
    g_lens = de_vaucouleurs_2d(xi1, xi2, vpar)

    # pl.figure()
    # pl.contourf(g_lens)
    # pl.colorbar()

    g_clean_ccd = g_lens * 0.0 + g_limage
    output_filename = "./output_fits/clean_lensed_imgs.fits"
    pyfits.writeto(output_filename, g_clean_ccd, clobber=True)
    from scipy.ndimage.filters import gaussian_filter
    #-------------------------------------------------------------
    g_images_psf = gaussian_filter(g_clean_ccd, 2.0)
    #g_images_psf = ss.convolve(g_clean_ccd,g_psf,mode="same")
    #g_images_psf = g_clean_ccd
    #-------------------------------------------------------------
    # Need to be Caliborate the mags
    g_noise = noise_map(nnn, nnn, np.sqrt(nstd), "Gaussian")
    output_filename = "./output_fits/noise_map.fits"
    pyfits.writeto(output_filename, g_noise, clobber=True)
    g_final = g_images_psf + g_noise
    #-------------------------------------------------------------
    output_filename = "./output_fits/lensed_imgs_only.fits"
    pyfits.writeto(output_filename, g_final, clobber=True)

    return xi1, xi2, yi1, yi2, kappa, mua

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

    kappa_out = np.loadtxt("./data_from_Quinn/kap_map.dat")

    pl.figure(figsize=(8, 8))
    pl.contour(xf1, xf2, np.log10(kappa_in), levels, colors=('k',))
    pl.contour(xf1, xf2, np.flipud(np.log10(kappa_out)), levels, colors=('r',))

    mua_out = np.loadtxt("./data_from_Quinn/mag_map.dat")

    pl.figure(figsize=(8, 8))
    pl.contour(xf1, xf2, mua_in, colors=('k',))
    pl.contour(xf1, xf2, np.flipud(mua_out), colors=('r',))

    pl.figure(figsize=(8, 8))
    pl.contour(yf1, yf2, mua_in, colors=('k',))
    pl.contour(yf1, yf2, np.flipud(mua_out), colors=('r',))

    #xpix1 = np.loadtxt("./data_from_Quinn/src_pixel.x")
    #xpix2 = np.loadtxt("./data_from_Quinn/src_pixel.y")

    vpix = np.loadtxt("./data_from_Quinn/src_pixel.dat")

    print np.shape(vpix)

    # pl.figure(figsize=(8,8))
    # pl.contourf(xpix1,xpix2,vpix)
    # pl.colorbar()

    pl.show()
