#!/usr/bin/env python
"""Fit stellar parameters on the coadded NEID spectrum with iSpec.

Run in a terminal (iSpec streams its own log to stdout / ispec.log):

    conda activate normal
    cd /Users/wangxianyu/Program/Github/NEID_iSpec
    python run_modeling.py

Initial guess = solar. vmic is set from the empirical relation; Teff, logg,
[Fe/H], vsini and vmac are all free. Fitting regions are the H-alpha/H-beta/Mg I
wings + Fe I/Fe II lines (segments_feh_Halpha_Hbeta_MgI.txt).
"""
import os
import numpy as np
import helper
# ignore any warnings
import warnings
warnings.filterwarnings("ignore")

# --- edit this to your iSpec installation (or export ISPEC_DIR in the shell) ---
ISPEC_DIR = os.environ.get('ISPEC_DIR', '/Users/wangxianyu/Program/Github/iSpec_v20230804')
if not os.path.exists(ISPEC_DIR): # /content/iSpec_v20230804
    ISPEC_DIR = os.environ.get('ISPEC_DIR', '/Users/wangxianyu/Program/Github/iSpec_v20230804')

helper.set_ispec_dir(ISPEC_DIR)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUT_DIR = os.path.join(BASE_DIR, 'output')
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR, exist_ok=True)
COADD = os.path.join(OUT_DIR, 'coadd_norm.txt')

# --- initial guess; free teff/logg/MH/alpha/vmic/vsini (matches iSpec GUI) ---
# alpha_mode: 'free' | 'auto'   vmic_mode: 'free' | 'GES'   vmac_mode: 'free'|'GES'|'Doyle2014'
# R fixed at RESOLUTION, LDC fixed at 0.6, vrad fixed at 0.

# initial Teff from the NEID L2 header (QTEFF); logg/MH/vsini stay generic
QTEFF = helper.header_teff(DATA_DIR)
print(f'initial Teff from header QTEFF = {QTEFF:.0f} K')
INIT = dict(teff=QTEFF, logg=4.44, MH=0.0, vsini=5.0)
ALPHA_MODE = 'free'
VMIC_MODE = 'free'
VMAC_MODE = 'Doyle2014'
RESOLUTION = 110000

if __name__ == '__main__':
    ispec = helper._ispec()
    # if coadded file doesn't exist, run get_coadded_spectra.py first
    if not os.path.exists(COADD):
        # print(f'Coadded spectrum {COADD} not found; run get_coadded_spectra.py first.')
        os.system('python get_coadded_spectra.py')
    spec = ispec.read_spectrum(COADD)
    print(f'Loaded {COADD}: {len(spec)} pixels, '
          f'{spec["waveobs"].min():.1f}-{spec["waveobs"].max():.1f} nm')

    # SPECTRUM + MARCS + GESv6 loaded inside fit_stellar_params (resources=None)
    params, errors, obs, model, status = helper.fit_stellar_params(
        spec, alpha_mode=ALPHA_MODE, vmic_mode=VMIC_MODE, vmac_mode=VMAC_MODE,
        resolution=RESOLUTION, **INIT)

    # mpfit return code: >0 = converged, <=0 = problem
    conv_code = status.get('status', None)
    converged = conv_code is not None and conv_code > 0

    print('\n===== RESULT =====')
    for k in ['teff', 'logg', 'MH', 'vmic', 'vmac', 'vsini']:
        if k in params:
            print(f'{k:6s} = {params[k]:9.3f} +/- {errors.get(k, float("nan")):.3f}')
    print(f'rchisq = {status.get("rchisq", float("nan")):.3f}, '
          f'rms = {status.get("rms", float("nan")):.4f}, '
          f'niter = {status.get("niter", "?")}, dof = {status.get("dof", "?")}')
    print(f'convergence: code={conv_code} ({"CONVERGED" if converged else "CHECK"}), '
          f'msg="{status.get("error", "")}"')

    # --- save: parameter table (with convergence status) ---
    with open(os.path.join(OUT_DIR, 'stellar_params.txt'), 'w') as f:
        f.write('# NEID + iSpec stellar parameters (XO-3 coadd)\n')
        f.write(f'# initial guess: solar; free: teff,logg,MH,vsini; '
                f'vmic_mode={VMIC_MODE}; vmac_mode={VMAC_MODE}\n')
        f.write('param\tvalue\terror\n')
        for k in ['teff', 'logg', 'MH', 'vmic', 'vmac', 'vsini']:
            if k in params:
                f.write(f'{k}\t{params[k]:.4f}\t{errors.get(k, np.nan):.4f}\n')
        f.write('#\n# --- convergence / fit status ---\n')
        for key in ['status', 'error', 'niter', 'nsynthesis', 'dof',
                    'chisq', 'rchisq', 'wchisq', 'rwchisq', 'rms']:
            if key in status:
                f.write(f'# {key}\t{status[key]}\n')
        f.write(f'# converged\t{converged}\n')

    # --- save: coadded (observed) spectrum used, and best-fit model ---
    ispec.write_spectrum(obs, os.path.join(OUT_DIR, 'coadd_used.txt'))
    ispec.write_spectrum(model, os.path.join(OUT_DIR, 'bestfit_model.txt'))

    # --- interactive zoomable HTML of obs vs best-fit + residuals ---
    html = helper.save_fit_html(obs, model, os.path.join(OUT_DIR, 'fit_check.html'),
                                params=params, status=status)
    print('\nSaved: output/stellar_params.txt, output/coadd_used.txt, '
          'output/bestfit_model.txt, ' + os.path.basename(html))
