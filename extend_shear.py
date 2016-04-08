def deflect_SIE(lens, x, y):  # SIE lens model
    tr = pi * (lens.th / 180.0) + pi / 2.0
    sx = x − lens.x0
    sy = y − lens.y0
    cs = cos(tr)
    sn = sin(tr)
    sx r = sx∗cs + sy∗sn
    sy r = −sx∗sn + sy∗cs
    psi = sqrt(lens.fl∗∗2.0 ∗ (lens.rc∗∗2.0 + sx r∗∗2.0) + sy r∗∗2.0)
    dx tmp = (lens.bl∗sqrt(lens.fl) / sqrt(1.0−lens.fl∗∗2.0))∗ arctan(sqrt(1.0−lens.fl∗∗2.0)∗sx r / (psi + lens.rc))
    dy tmp = (lens.bl∗sqrt(lens.fl) / sqrt(1.0−lens.fl∗∗2.0))∗arctanh(sqrt(1.0−lens.fl∗∗2.0)∗sy r / (psi + lens.rc∗lens.fl∗∗2.0))
    dx = dx tmp∗cs − dy tmp∗sn dy = dx tmp∗sn + dy tmp∗cs
    # external shear
    tr2 = cs2 = sn2 =
    dx2 = dy2 =
    pi∗(lens.sa / 180.0) cos(2.0∗tr2) sin(2.0∗tr2)
    lens.ss∗(cs2∗sx + sn2∗sy)
    lens.ss∗(sn2∗sx−cs2∗sy)
    return array([dx + dx2, dy + dy2])
#############################################################################
# Convergence for SIE + external shear #
#############################################################################


def convergence_SIE(lens, x, y):
    tr = pi ∗(lens . th / 180.0) + pi / 2.0 72
    sx = x−lens.x0
    sy = y−lens.y0
    cs = cos(tr)
    sn = sin(tr)
    sx r = sx∗cs + sy∗sn sy r = −sx∗sn + sy∗cs
    psi = sqrt(lens.fl∗∗2.0 ∗ (lens.rc∗∗2.0 + sx r∗∗2.0) + sy r∗∗2.0) kappa tmp = (0.5∗lens.bl∗sqrt(lens.fl) / psi)
    return kappa tmp
#
# Potential for SIE + external shear #
#############################################################################


def potential(lens, x, y):  # SIE lens model
    tr = pi ∗(lens . th / 180.0) + pi / 2.0 sx = x−lens.x0
    sy = y−lens.y0 cs = cos(tr)
    sn = sin(tr)
    sx r = sx∗cs + sy∗sn
    sy r = −sx∗sn + sy∗cs
    psi = sqrt(lens.fl∗∗2.0 ∗ (lens.rc∗∗2.0 + sx r∗∗2.0) + sy r∗∗2.0)
    dx tmp = (lens.bl∗sqrt(lens.fl) / sqrt(1.0−lens.fl∗∗2.0))∗ arctan(sqrt(1.0−lens.fl∗∗2.0)∗sx r / (psi + lens.rc))
    dy tmp = (lens.bl∗sqrt(lens.fl) / sqrt(1.0−lens.fl∗∗2.0))∗arctanh(sqrt(1.0−lens.fl∗∗2.0)∗sy r / (psi + lens.rc∗lens.fl∗∗2.0))
    pot SIE = sx r∗dx tmp + sy r∗dy tmp − 0.5∗lens.bl∗sqrt(lens.fl)∗ lens.rc∗log((psi + lens.rc)∗∗2.0 + (1.0−(lens.fl∗∗2.0))∗(sx r∗∗2.0))

    # external shear
    tr2 = pi∗(lens.sa / 180.0) cs2 = cos(2.0∗ tr2)
    sn2 = sin(2.0∗tr2)
    pot exts = lens . ss∗(sn2∗sx∗sy + 0.5∗cs2∗(sx∗∗2.0−sy∗∗2.0))
    return pot_SIE + pot_exts


def phi_shears(x1, x2):
    res = shear1 / 2.0 * (x1 * x1 - x2 * x2) + shear2 * x1 * x2
    return res


def phi_kappa(x1, x2):
    res = kappa / 2.0 * (x1 * x1 + x2 * x2)
    return res


def potential_shear_kappa(x1, x2):
    phi_all = phi_sie(x1, x2) + phi_shear(x1, x2) + phi_kappa(x1, x2)
    return phi_all
