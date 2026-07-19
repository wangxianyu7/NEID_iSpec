#!/usr/bin/env python
"""Build a coadded, rest-frame, normalized NEID spectrum for iSpec modeling.

Pipeline (all functions live in helper.py):
    deblaze (SCIBLAZE) -> rest frame (L2 header) -> cosmic removal + continuum
    normalization -> inverse-variance coadd.

Run in a terminal:

    conda activate normal
    cd /Users/wangxianyu/Program/Github/NEID_iSpec
    python Get_Coadded_Spectra.py

Reads data/neidL2_*.fits, writes output/coadd_norm.txt and diagnostic PNGs.
"""
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')            # headless: save figures, don't pop up windows
import matplotlib.pyplot as plt
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")
import helper
from helper import (deblaze_neid, to_rest_frame, clean_and_normalize,
                    coadd_spectra)

# --- edit this to your iSpec installation (or export ISPEC_DIR in the shell) ---
ISPEC_DIR = os.environ.get('ISPEC_DIR', '/Users/wangxianyu/Program/Github/iSpec_v20230804')
if not os.path.exists(ISPEC_DIR): # /content/iSpec_v20230804
    ISPEC_DIR = os.environ.get('ISPEC_DIR', '/content/iSpec_v20230804')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUT_DIR = os.path.join(BASE_DIR, 'output')
os.makedirs(OUT_DIR, exist_ok=True)

# Continuum normalization: 'template' picks a FIXED template a priori from the
# header Teff (QTEFF) -- like choosing a CCF mask by spectral type -- and divides
# out the line-blanketing so the continuum stops riding into crowded line forests.
# 'splines' is the model-independent fallback. Template logg/[M/H] are rough on
# purpose: they set the continuum shape, not the fit.
CONTINUUM = 'template'          # 'template' | 'splines'
TEMPLATE_LOGG = 4.0
TEMPLATE_MH = 0.0
RESOLUTION = 110000


def main():
    L2_FILES = sorted(glob.glob(os.path.join(DATA_DIR, 'neidL2_*.fits')))
    if not L2_FILES:
        raise SystemExit(f'No NEID L2 files in {DATA_DIR}')
    print(f'{len(L2_FILES)} NEID L2 files:')
    for f in L2_FILES:
        print('  ', os.path.basename(f), '->', fits.getheader(f)['OBJECT'])

    # a-priori continuum template from the header Teff (not from any fit)
    template = None
    if CONTINUUM == 'template':
        qteff = float(fits.getheader(L2_FILES[0])['QTEFF'])
        print(f'building continuum template at header QTEFF={qteff:.0f} K, '
              f'logg={TEMPLATE_LOGG}, [M/H]={TEMPLATE_MH} ...')
        template = helper.make_continuum_template(
            qteff, logg=TEMPLATE_LOGG, MH=TEMPLATE_MH, resolution=RESOLUTION)

    # deblaze -> rest frame -> clean + normalize, per exposure
    norm_specs = []
    orders0 = None
    for f in L2_FILES:
        name = os.path.basename(f)
        print(name)
        orders = deblaze_neid(f)
        if orders0 is None:
            orders0 = orders                      # keep first for the deblaze plot
        w, fl, er, info = to_rest_frame(orders, f)
        nspec, ncos = clean_and_normalize(w, fl, er, resolution=RESOLUTION,
                                          template=template)
        norm_specs.append((nspec, name))

    # coadd
    coadd = coadd_spectra(norm_specs)
    helper._ispec().write_spectrum(coadd, os.path.join(OUT_DIR, 'coadd_norm.txt'))
    print('saved output/coadd_norm.txt')

    _diagnostics(orders0, norm_specs, coadd)


def _diagnostics(orders, norm_specs, coadd, n_junction=60):
    """Save deblaze / normalize / coadd check figures to output/."""
    # 1) deblazed orders + one clean order junction
    fig, ax = plt.subplots(2, 1, figsize=(13, 6))
    for od in orders:
        ax[0].plot(od['wave_nm'], od['flux'], lw=0.4)
    ax[0].set_title('Deblazed + edge-trimmed orders (all)')
    ax[0].set_ylabel('flux / blaze')
    oa = [o for o in orders if o['order'] == n_junction][0]
    ob = [o for o in orders if o['order'] == n_junction + 1][0]
    wc = 0.5 * (oa['wave_nm'].max() + ob['wave_nm'].min())
    w0, w1 = wc - 1.5, wc + 1.5
    for od in orders:
        m = (od['wave_nm'] > w0) & (od['wave_nm'] < w1)
        if m.sum() > 10:
            ax[1].plot(od['wave_nm'][m], od['flux'][m], lw=0.7, label=f'ord {od["order"]}')
    ax[1].set_xlim(w0, w1)
    ax[1].set_xlabel('wavelength [nm, air]'); ax[1].set_ylabel('flux / blaze')
    ax[1].set_title(f'Order junction {n_junction}/{n_junction+1} (~{wc:.1f} nm)')
    ax[1].legend(fontsize=8, ncol=4)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'deblaze_check.png'), dpi=120, bbox_inches='tight')
    plt.close(fig)

    # 2) normalized exposures (full + Mg b zoom)
    fig, ax = plt.subplots(2, 1, figsize=(13, 6))
    for nspec, name in norm_specs:
        ax[0].plot(nspec['waveobs'], nspec['flux'], lw=0.4, label=name)
    ax[0].axhline(1, color='k', lw=0.6, ls='--'); ax[0].set_ylim(0, 1.2)
    ax[0].set_ylabel('normalized flux'); ax[0].set_title('Normalized (full range)')
    ax[0].legend(fontsize=8)
    for nspec, name in norm_specs:
        m = (nspec['waveobs'] > 516) & (nspec['waveobs'] < 519)
        ax[1].plot(nspec['waveobs'][m], nspec['flux'][m], lw=0.7, label=name)
    ax[1].axhline(1, color='k', lw=0.6, ls='--'); ax[1].set_ylim(0, 1.15)
    ax[1].set_xlabel('wavelength [nm, air, rest]'); ax[1].set_ylabel('normalized flux')
    ax[1].set_title('Zoom: Mg b (516-519 nm)'); ax[1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'normalize_check.png'), dpi=120, bbox_inches='tight')
    plt.close(fig)

    # 3) coadd vs individual (Mg b)
    fig, ax = plt.subplots(figsize=(13, 4))
    for nspec, name in norm_specs:
        m = (nspec['waveobs'] > 516) & (nspec['waveobs'] < 519)
        ax.plot(nspec['waveobs'][m], nspec['flux'][m], lw=0.5, alpha=0.5, label=name)
    m = (coadd['waveobs'] > 516) & (coadd['waveobs'] < 519)
    ax.plot(coadd['waveobs'][m], coadd['flux'][m], 'k-', lw=0.8, label='coadd')
    ax.axhline(1, color='gray', lw=0.5, ls='--'); ax.set_ylim(0, 1.15)
    ax.set_xlabel('wavelength [nm, air, rest]'); ax.set_ylabel('normalized flux')
    ax.set_title('Coadd vs individual exposures (Mg b)'); ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'coadd_check.png'), dpi=120, bbox_inches='tight')
    plt.close(fig)
    print('saved diagnostic figures to output/')


if __name__ == '__main__':
    main()
