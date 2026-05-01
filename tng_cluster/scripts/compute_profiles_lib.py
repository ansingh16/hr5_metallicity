"""
Shared library functions for computing radial metallicity profiles.

Used by both 04_compute_profiles.py and 05_redshift_evolution.py.
"""

import os
import sys
import numpy as np
import h5py
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

cosmo = FlatLambdaCDM(H0=cfg.H0, Om0=cfg.OMEGA_M, Ob0=cfg.OMEGA_B)


def compute_abundance_ratio(metals_array, elem_idx, ref_idx, solar_elem, solar_ref):
    x = metals_array[:, elem_idx]
    y = metals_array[:, ref_idx]
    ratio = np.full(len(x), np.nan)
    valid = (x > 0) & (y > 0)
    ratio[valid] = np.log10((x[valid] / y[valid]) / (solar_elem / solar_ref))
    return ratio


def compute_feh(metals_array):
    return compute_abundance_ratio(
        metals_array, cfg.ELEM_FE, cfg.ELEM_H,
        cfg.FE_MASS_FRAC_SUN, cfg.H_MASS_FRAC_SUN
    )


def compute_ofe(metals_array):
    return compute_abundance_ratio(
        metals_array, cfg.ELEM_O, cfg.ELEM_FE,
        cfg.O_MASS_FRAC_SUN, cfg.FE_MASS_FRAC_SUN
    )


def compute_mgfe(metals_array):
    return compute_abundance_ratio(
        metals_array, cfg.ELEM_MG, cfg.ELEM_FE,
        cfg.MG_MASS_FRAC_SUN, cfg.FE_MASS_FRAC_SUN
    )


def compute_no(metals_array):
    return compute_abundance_ratio(
        metals_array, cfg.ELEM_N, cfg.ELEM_O,
        cfg.N_MASS_FRAC_SUN, cfg.O_MASS_FRAC_SUN
    )


def compute_co(metals_array):
    return compute_abundance_ratio(
        metals_array, cfg.ELEM_C, cfg.ELEM_O,
        cfg.C_MASS_FRAC_SUN, cfg.O_MASS_FRAC_SUN
    )


def compute_cn(metals_array):
    return compute_abundance_ratio(
        metals_array, cfg.ELEM_C, cfg.ELEM_N,
        cfg.C_MASS_FRAC_SUN, cfg.N_MASS_FRAC_SUN
    )


def radial_profile(values, weights, radii, r_half, n_bins, r_max_rhalf):
    r_norm = radii / r_half
    bin_edges = np.linspace(0, r_max_rhalf, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    profile = np.full(n_bins, np.nan)
    profile_err = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        mask = (r_norm >= bin_edges[i]) & (r_norm < bin_edges[i + 1])
        mask &= np.isfinite(values)
        n = mask.sum()
        counts[i] = n
        if n < 5:
            continue
        w = weights[mask]
        v = values[mask]
        w_sum = w.sum()
        if w_sum <= 0:
            continue
        wmean = np.average(v, weights=w)
        wvar = np.average((v - wmean) ** 2, weights=w)
        profile[i] = wmean
        profile_err[i] = np.sqrt(wvar / n)

    return bin_centers, profile, profile_err, counts


def fit_gradient(bins, profile, r_min, r_max):
    mask = (bins >= r_min) & (bins <= r_max) & np.isfinite(profile)
    if mask.sum() < 3:
        return np.nan, np.nan, np.nan
    x = bins[mask]
    y = profile[mask]
    A = np.vstack([x, np.ones(len(x))]).T
    result = np.linalg.lstsq(A, y, rcond=None)
    slope, intercept = result[0]
    residuals = y - (slope * x + intercept)
    n = len(x)
    if n > 2:
        s2 = np.sum(residuals ** 2) / (n - 2)
        slope_err = np.sqrt(s2 / np.sum((x - x.mean()) ** 2))
    else:
        slope_err = np.nan
    return slope, intercept, slope_err


def process_galaxy(filepath, subhalo_id, r_half_pkpc, z):
    a = 1.0 / (1.0 + z)
    age_universe = cosmo.age(z).to(u.Gyr).value

    try:
        f = h5py.File(filepath, "r")
    except Exception as e:
        print(f"  Cannot open {filepath}: {e}")
        return None

    result = {"subhalo_id": subhalo_id}

    if "PartType4" not in f:
        f.close()
        return None

    stars = f["PartType4"]
    coords = stars["Coordinates"][:] * a / cfg.h
    masses = stars["Masses"][:] * 1e10 / cfg.h
    metals = stars["GFM_Metals"][:]
    metallicity = stars["GFM_Metallicity"][:]
    a_form = stars["GFM_StellarFormationTime"][:]

    star_mask = a_form > 0
    coords = coords[star_mask]
    masses = masses[star_mask]
    metals = metals[star_mask]
    metallicity = metallicity[star_mask]
    a_form = a_form[star_mask]

    if len(masses) < 50:
        f.close()
        return None

    com = np.average(coords, axis=0, weights=masses)
    dx = coords - com
    radii = np.sqrt(np.sum(dx ** 2, axis=1))

    z_form = 1.0 / a_form - 1.0
    ages = age_universe - np.array([cosmo.age(zf).to(u.Gyr).value for zf in z_form])

    feh = compute_feh(metals)
    ofe = compute_ofe(metals)
    mgfe = compute_mgfe(metals)
    no_ratio = compute_no(metals)
    co_ratio = compute_co(metals)
    cn_ratio = compute_cn(metals)

    bins, prof_feh, err_feh, cnt = radial_profile(
        feh, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_ofe, err_ofe, _ = radial_profile(
        ofe, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_mgfe, err_mgfe, _ = radial_profile(
        mgfe, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_no, err_no, _ = radial_profile(
        no_ratio, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_co, err_co, _ = radial_profile(
        co_ratio, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_cn, err_cn, _ = radial_profile(
        cn_ratio, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_zstar, err_zstar, _ = radial_profile(
        np.log10(metallicity / cfg.Z_SUN + 1e-30), masses, radii,
        r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
    _, prof_age, err_age, _ = radial_profile(
        ages, masses, radii, r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)

    result["r_bins"] = bins
    result["star_counts"] = cnt
    result["FeH_star"] = prof_feh
    result["FeH_star_err"] = err_feh
    result["OFe_star"] = prof_ofe
    result["OFe_star_err"] = err_ofe
    result["MgFe_star"] = prof_mgfe
    result["MgFe_star_err"] = err_mgfe
    result["NO_star"] = prof_no
    result["NO_star_err"] = err_no
    result["CO_star"] = prof_co
    result["CO_star_err"] = err_co
    result["CN_star"] = prof_cn
    result["CN_star_err"] = err_cn
    result["Z_star"] = prof_zstar
    result["Z_star_err"] = err_zstar
    result["age_star"] = prof_age
    result["age_star_err"] = err_age

    # Global mass-weighted properties
    valid_feh = np.isfinite(feh)
    if valid_feh.sum() > 0:
        result["FeH_star_global"] = np.average(feh[valid_feh], weights=masses[valid_feh])
    valid_ofe = np.isfinite(ofe)
    if valid_ofe.sum() > 0:
        result["OFe_star_global"] = np.average(ofe[valid_ofe], weights=masses[valid_ofe])
    valid_mgfe = np.isfinite(mgfe)
    if valid_mgfe.sum() > 0:
        result["MgFe_star_global"] = np.average(mgfe[valid_mgfe], weights=masses[valid_mgfe])
    valid_no = np.isfinite(no_ratio)
    if valid_no.sum() > 0:
        result["NO_star_global"] = np.average(no_ratio[valid_no], weights=masses[valid_no])
    valid_co = np.isfinite(co_ratio)
    if valid_co.sum() > 0:
        result["CO_star_global"] = np.average(co_ratio[valid_co], weights=masses[valid_co])
    valid_cn = np.isfinite(cn_ratio)
    if valid_cn.sum() > 0:
        result["CN_star_global"] = np.average(cn_ratio[valid_cn], weights=masses[valid_cn])
    # Mass-weighted mean stellar age
    result["age_star_global"] = np.average(ages, weights=masses)

    if "PartType0" in f and f["PartType0"]["Coordinates"].shape[0] > 10:
        gas = f["PartType0"]
        g_coords = gas["Coordinates"][:] * a / cfg.h
        g_masses = gas["Masses"][:] * 1e10 / cfg.h
        g_metals = gas["GFM_Metals"][:]
        g_metallicity = gas["GFM_Metallicity"][:]
        g_sfr = gas["StarFormationRate"][:]

        g_dx = g_coords - com
        g_radii = np.sqrt(np.sum(g_dx ** 2, axis=1))
        g_feh = compute_feh(g_metals)

        _, prof_zgas, err_zgas, g_cnt = radial_profile(
            np.log10(g_metallicity / cfg.Z_SUN + 1e-30), g_masses, g_radii,
            r_half_pkpc, cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
        _, prof_feh_gas, err_feh_gas, _ = radial_profile(
            g_feh, g_masses, g_radii, r_half_pkpc,
            cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)
        _, prof_sfr, err_sfr, _ = radial_profile(
            g_sfr, g_masses, g_radii, r_half_pkpc,
            cfg.N_RADIAL_BINS, cfg.MAX_RADIUS_RHALF)

        result["gas_counts"] = g_cnt
        result["Z_gas"] = prof_zgas
        result["Z_gas_err"] = err_zgas
        result["FeH_gas"] = prof_feh_gas
        result["FeH_gas_err"] = err_feh_gas
        result["SFR_profile"] = prof_sfr
        result["SFR_profile_err"] = err_sfr
        result["has_gas"] = True
    else:
        result["has_gas"] = False

    f.close()
    return result
