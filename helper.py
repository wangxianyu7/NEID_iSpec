"""Helper functions for the NEID + iSpec stellar-parameter pipeline.

Currently: NEID L2 deblazing (SCIBLAZE profile) with order-edge trimming and
the 77/78 junction offset correction. More stages (RV correction, continuum
normalization, coaddition, iSpec fitting) will be added as they are validated
in main.ipynb.
"""
import os
import logging
import numpy as np
from astropy.io import fits

# hush numexpr's "NumExpr defaulting to N threads" INFO line (emitted on import)
logging.getLogger('numexpr').setLevel(logging.WARNING)
import warnings
warnings.filterwarnings("ignore")
# NEID echelle orders covering ~480-680 nm (matches the SPECTRUM_MARCS grid)
ORD_LO, ORD_HI = 45, 82


def vac2air(wl_vac_A):
    """Vacuum -> air wavelength (Birch & Downs 1994 / Morton). Angstrom in/out."""
    s = 1e4 / wl_vac_A
    n = (1 + 0.0000834254 + 0.02406147 / (130 - s**2) + 0.00015998 / (38.9 - s**2))
    return wl_vac_A / n


def order_window(order):
    """Central pixel window per order: trims noisy blaze edges / order overlaps.

    Asymmetric, widening redward (empirical, NEID 9216-pix format; from the
    TOI-4468 pipeline). Returns (start_idx, end_idx).
    """
    start = int(4608 - 2300 - 500 - (order - 46) * 10)
    end = int(4608 + 2300 - 1000 + (order - 46) * 30)
    return start, end


def deblaze_neid(filename, ord_lo=ORD_LO, ord_hi=ORD_HI, stitch_offset=True):
    """Deblaze one NEID L2 with the SCIBLAZE profile, trim edges, and stitch.

    Special handling at order junctions (matches the prior TOI-4468 pipeline):
      1. keep only the central window of each order (drop blaze edges/overlaps);
      2. correct the known flux jump at the 77/78 junction by scaling orders >77
         to match orders <=77 in the 641.6-642.0 nm overlap.

    SCIWAVE is vacuum, barycentric-frame; RV/barycentric handling is downstream.

    Returns a list of per-order dicts: {order, wave_nm (air), flux, err},
    with flux = SCIFLUX / SCIBLAZE and err = sqrt(SCIVAR) / SCIBLAZE.
    """
    with fits.open(filename) as h:
        flux = h['SCIFLUX'].data
        var = h['SCIVAR'].data
        wave = h['SCIWAVE'].data          # vacuum, Angstrom
        blaze = h['SCIBLAZE'].data

    # --- 77/78 offset from the overlap (in deblazed flux) ---
    offset = 1.0
    if stitch_offset:
        def _ov(o):
            w = vac2air(wave[o]) / 10.0
            fdb = flux[o] / blaze[o]
            m = (w > 641.6) & (w < 642.0) & np.isfinite(fdb)
            if m.sum() < 5:
                return None
            wr = np.arange(641.6, 642.0, 1e-4)
            return np.interp(wr, w[m], fdb[m])
        f77, f78 = _ov(77), _ov(78)
        if f77 is not None and f78 is not None:
            offset = float(np.nanmedian(f77 / f78))

    orders = []
    for o in range(ord_lo, ord_hi + 1):
        w, fl, bl, vr = wave[o], flux[o], blaze[o], var[o]
        s, e = order_window(o)
        s, e = max(s, 0), min(e, len(w))
        w, fl, bl, vr = w[s:e], fl[s:e], bl[s:e], vr[s:e]
        good = (np.isfinite(w) & (w > 0) & np.isfinite(fl) &
                np.isfinite(bl) & (bl > 0))
        if good.sum() < 100:
            continue
        w, fl, bl, vr = w[good], fl[good], bl[good], vr[good]
        f_db = fl / bl
        e_db = np.sqrt(np.clip(vr, 0, None)) / bl
        if stitch_offset and o > 77:
            f_db = f_db * offset
            e_db = e_db * offset
        orders.append({'order': o, 'wave_nm': vac2air(w) / 10.0,
                       'flux': f_db, 'err': e_db})
    return orders


C_KMS = 299792.458   # speed of light [km/s]


def join_orders(orders):
    """Concatenate per-order deblaze output into sorted (wave_nm, flux, err) arrays."""
    w = np.concatenate([o['wave_nm'] for o in orders])
    f = np.concatenate([o['flux'] for o in orders])
    e = np.concatenate([o['err'] for o in orders])
    s = np.argsort(w)
    return w[s], f[s], e[s]


def to_rest_frame(orders, filename, verbose=True):
    """Shift the joined spectrum to the stellar rest frame using NEID L2 header values.

    NEID SCIWAVE is topocentric (observer-frame), vacuum. The header carries the
    barycentric redshift (SSBZ100) and the DRP systemic RV (QRV), so we shift
    analytically instead of running a cross-correlation:

        lambda_rest = lambda_obs * (1 + z_bary) / (1 + RV_sys/c)

    (multiply by 1+z_bary because the observer moves; divide by 1+z_sys because
    the star moves -- same convention as neidspecmatch.CombineNEIDSpectra.)

    Input: `orders` = output of deblaze_neid; `filename` = the NEID L2 FITS.
    Returns (wave_nm_rest, flux, err, info).
    """
    hdr = fits.getheader(filename)
    z_bary = float(hdr['SSBZ100'])            # barycentric redshift
    rv_sys = float(hdr['QRV'])                # DRP systemic RV [km/s]
    z_bulk = rv_sys / C_KMS
    factor = (1.0 + z_bary) / (1.0 + z_bulk)

    w, f, e = join_orders(orders)
    w_rest = w * factor
    info = {'z_bary': z_bary, 'berv_kms': z_bary * C_KMS,
            'rv_sys_kms': rv_sys, 'net_shift_kms': (factor - 1.0) * C_KMS}
    if verbose:
        print(f'  BERV = {info["berv_kms"]:+.3f} km/s, RV_sys = {rv_sys:+.3f} km/s '
              f'-> rest-frame shift {info["net_shift_kms"]:+.3f} km/s')
    return w_rest, f, e, info


def header_teff(path, key='QTEFF', default=5771.0):
    """Initial Teff from a NEID L2 header (QTEFF = queue/catalog Teff). `path` may
    be a FITS file or a directory of neidL2_*.fits (uses the first). Falls back to
    `default` (solar) if the file/keyword is missing or blank."""
    import glob
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, 'neidL2_*.fits')))
        if not files:
            return default
        path = files[0]
    try:
        v = fits.getheader(path).get(key, default)
        return float(v) if str(v).strip() != '' else default
    except (OSError, ValueError, TypeError):
        return default


# Path to the iSpec installation. Not hardcoded: set it from the calling script
# (helper.set_ispec_dir(...)) or via the ISPEC_DIR environment variable.
ISPEC_DIR = os.environ.get('ISPEC_DIR', '')

# iSpec input files needed for the paper's setup, RELATIVE to ISPEC_DIR (stable
# across installs); joined with ISPEC_DIR lazily so the setter takes effect.
_STRONG_LINES_REL = 'input/regions/relevant/relevant_line_masks.txt'
_MODEL_ATMOSPHERE_REL = 'input/atmospheres/MARCS.GES/'
_ATOMIC_LINELIST_REL = 'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv'
_ISOTOPE_REL = 'input/isotopes/SPECTRUM.lst'
_SOLAR_ABUNDANCES_REL = 'input/abundances/Grevesse.2007/stdatom.dat'


def set_ispec_dir(path):
    """Point the pipeline at an iSpec installation (call before any fit/synth)."""
    global ISPEC_DIR
    ISPEC_DIR = path


def resolve_ispec_dir(extra_candidates=()):
    """Locate the iSpec install and register it. Order: $ISPEC_DIR, any
    extra_candidates, then a few common locations (Colab /content, home, this
    machine). Registers it via set_ispec_dir AND exports $ISPEC_DIR so that
    child processes (e.g. os.system('python get_coadded_spectra.py')) inherit it.
    Returns the path; raises if none exist.
    """
    candidates = ([os.environ.get('ISPEC_DIR', '')] + list(extra_candidates) +
                  ['/content/iSpec_v20230804',
                   os.path.expanduser('~/iSpec_v20230804'),
                   '/Users/wangxianyu/Program/Github/iSpec_v20230804'])
    for c in candidates:
        if c and os.path.isdir(c):
            set_ispec_dir(c)
            os.environ['ISPEC_DIR'] = c
            return c
    raise RuntimeError(
        'iSpec install not found. Set the ISPEC_DIR environment variable, or pass '
        'helper.resolve_ispec_dir(extra_candidates=["/path/to/iSpec_v20230804"]).')


def _ispec_path(rel):
    """Full path to an iSpec input file, resolved from the current ISPEC_DIR."""
    return os.path.join(ISPEC_DIR, rel)


def _ispec():
    import sys
    if not ISPEC_DIR:
        raise RuntimeError(
            'ISPEC_DIR is not set. Call helper.set_ispec_dir("/path/to/iSpec") '
            'or set the ISPEC_DIR environment variable before fitting.')
    if ISPEC_DIR not in sys.path:
        sys.path.insert(0, ISPEC_DIR)
    import ispec
    return ispec


def make_continuum_template(teff, logg=4.0, MH=0.0, vsini=10.0, resolution=110000,
                            wave_base=480.0, wave_top=680.0, wave_step=0.001,
                            resources=None):
    """Synthesize a FIXED template for continuum normalization, chosen a priori by
    spectral type / header Teff (analogous to picking a CCF mask by spectral type).

    Because Teff/logg/[Fe/H] come from the header (QTEFF), not from the fit, this
    is not circular. The template carries the line-blanketing that a blind spline
    cannot recover in crowded regions (it stops the continuum from riding into
    deep line forests). Returns an iSpec spectrum on [wave_base, wave_top].
    """
    ispec = _ispec()
    if resources is None:
        common = np.arange(wave_base, wave_top, wave_step)
        dummy = ispec.create_spectrum_structure(common)
        resources = load_rt_resources(dummy)
    layers, linelist, isotopes, abund = resources
    alpha = ispec.determine_abundance_enchancements(MH)
    vmic = ispec.estimate_vmic(teff, logg, MH)
    vmac = ispec.estimate_vmac(teff, logg, MH, relation='Doyle2014')
    atm = ispec.interpolate_atmosphere_layers(
        layers, {'teff': teff, 'logg': logg, 'MH': MH, 'alpha': alpha}, code=RT_CODE)
    tmpl = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    tmpl['flux'] = ispec.generate_spectrum(
        tmpl['waveobs'], atm, teff, logg, MH, alpha, linelist, isotopes, abund,
        fixed_abundances=None, microturbulence_vel=vmic, macroturbulence=vmac,
        vsini=vsini, limb_darkening_coeff=0.6, R=resolution, verbose=0, code=RT_CODE)
    return tmpl


def clean_and_normalize(wave_nm, flux, err, resolution=110000, variation_limit=0.30,
                        template=None, verbose=True):
    """Cosmic-ray removal + continuum normalization.

    Cosmic filter needs a rough continuum first (its median+max filtering is
    spike-robust), so: rough continuum -> filter cosmics -> final continuum ->
    normalize to 1.

    Continuum model for the final fit:
      * template is None -> model-independent Splines (median+max, strong lines
        ignored). Fast, presumes no stellar params, but can ride into deep line
        forests in crowded regions.
      * template given   -> model='Template' against that FIXED template (chosen
        a priori from the header Teff, see make_continuum_template). Recovers the
        line-blanketing pseudo-continuum in crowded regions; not circular because
        the template is not the fit result.
    Returns (norm_spec, n_cosmics).
    """
    ispec = _ispec()
    spec = ispec.create_spectrum_structure(wave_nm, flux.astype(float), err.astype(float))
    spec = spec[np.isfinite(spec['flux'])]

    # 1) rough, cosmic-robust continuum
    rough = ispec.fit_continuum(
        spec, from_resolution=resolution, model='Splines', order='median+max',
        median_wave_range=0.05, max_wave_range=1.0,
        automatic_strong_line_detection=True, strong_line_probability=0.5,
        use_errors_for_fitting=True)

    # 2) drop positive cosmic-ray spikes
    cosmics = ispec.create_filter_cosmic_rays(
        spec, rough, resampling_wave_step=0.001, window_size=15,
        variation_limit=variation_limit)
    n_cosmics = int(np.sum(cosmics))
    spec = spec[~cosmics]

    # 3) final continuum on the cleaned spectrum
    strong = ispec.read_line_regions(_ispec_path(_STRONG_LINES_REL))
    if template is not None:
        cont = ispec.fit_continuum(
            spec, from_resolution=resolution, ignore=strong, nknots=None,
            median_wave_range=5, model='Template', template=template)
    else:
        cont = ispec.fit_continuum(
            spec, from_resolution=resolution, ignore=strong, model='Splines',
            order='median+max', median_wave_range=0.05, max_wave_range=1.0,
            automatic_strong_line_detection=True, strong_line_probability=0.5,
            use_errors_for_fitting=True)

    # 4) normalize
    norm = ispec.normalize_spectrum(spec, cont, consider_continuum_errors=False)
    if verbose:
        mode = 'template' if template is not None else 'splines'
        print(f'  removed {n_cosmics} cosmic-ray pixels; normalized ({mode})')
    return norm, n_cosmics


def coadd_spectra(norm_specs, wave_base=480.0, wave_top=680.0, wave_step=0.001,
                  verbose=True):
    """Resample (iSpec) each normalized spectrum onto a common grid and combine.

    Resampling uses ispec.resample_spectrum; the combine is an inverse-variance
    weighted mean (falls back to a plain mean where errors are zero/invalid).
    Since the inputs are continuum-normalized, the weighted mean stays at ~1.

    Input: norm_specs = list of (ispec_spectrum, name).
    Returns the coadded iSpec spectrum on `common_wave`.
    """
    ispec = _ispec()
    common = np.arange(wave_base, wave_top, wave_step)
    fluxes, weights = [], []
    for spec, _ in norm_specs:
        r = ispec.resample_spectrum(spec, common, method='linear', zero_edges=True)
        f = r['flux'].astype(float)
        e = r['err'].astype(float)
        w = np.where(np.isfinite(e) & (e > 0), 1.0 / e**2, 0.0)
        f = np.where(np.isfinite(f), f, 0.0)
        fluxes.append(f)
        weights.append(w)
    F = np.vstack(fluxes)
    W = np.vstack(weights)
    wsum = W.sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        coflux = np.where(wsum > 0, (F * W).sum(axis=0) / wsum,
                          np.nanmean(np.where(F != 0, F, np.nan), axis=0))
        coerr = np.where(wsum > 0, 1.0 / np.sqrt(wsum), np.nan)

    good = np.isfinite(coflux) & (coflux > 0) & (coflux < 1.3)
    co = ispec.create_spectrum_structure(common[good])
    co['flux'] = coflux[good]
    co['err'] = coerr[good]
    if verbose:
        print(f'  coadded {len(norm_specs)} spectra -> {good.sum()} px, '
              f'{common[good].min():.1f}-{common[good].max():.1f} nm')
    return co


# --- radiative-transfer setup (SPECTRUM + MARCS + GES v6, following the paper) ---
RT_CODE = 'spectrum'                                                   # Gray 1994
SEGMENTS_FILE = os.path.join(os.path.dirname(__file__), 'segments_feh_Halpha_Hbeta_MgI.txt')


def load_rt_resources(spec):
    """Load the MARCS atmosphere pack, GESv6 line list, isotopes and solar
    abundances needed for on-the-fly SPECTRUM synthesis over `spec`'s range."""
    ispec = _ispec()
    modeled_layers_pack = ispec.load_modeled_layers_pack(_ispec_path(_MODEL_ATMOSPHERE_REL))
    atomic_linelist = ispec.read_atomic_linelist(
        _ispec_path(_ATOMIC_LINELIST_REL),
        wave_base=spec['waveobs'].min(), wave_top=spec['waveobs'].max())
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01]
    isotopes = ispec.read_isotope_data(_ispec_path(_ISOTOPE_REL))
    solar_abundances = ispec.read_solar_abundances(_ispec_path(_SOLAR_ABUNDANCES_REL))
    return modeled_layers_pack, atomic_linelist, isotopes, solar_abundances


def _set_vmac_relation(relation):
    """Runtime-patch the estimate_vmac used by the mpfit 'tied' evaluation so the
    empirical-relation mode can be GES or Doyle2014 -- WITHOUT editing iSpec.

    mpfit.tie() eval's the tied string 'estimate_vmac(p[0],p[1],p[2])' against the
    estimate_vmac imported into ispec.modeling.mpfit; we swap that name for a
    wrapper with the chosen relation. relation=None restores the GES default.
    """
    _ispec()                       # ensure ISPEC_DIR is on sys.path
    import ispec.modeling.mpfit as _mp
    from ispec.common import estimate_vmac as _orig
    if relation in (None, 'GES'):
        _mp.estimate_vmac = _orig
    else:
        _mp.estimate_vmac = lambda t, g, m: _orig(t, g, m, relation=relation)


def fit_stellar_params(spec, segments_file=SEGMENTS_FILE, resources=None,
                       teff=5771, logg=4.44, MH=0.0, vsini=5.0,
                       alpha_mode='free', vmic_mode='free', vmac_mode='Doyle2014',
                       resolution=110000, use_errors=False,
                       max_iterations=20, verbose=True):
    """Fit stellar parameters on a normalized (coadded) spectrum via iSpec.

    Radiative transfer: SPECTRUM (Gray 1994) with the MARCS atmosphere model
    (Gustafsson 2008) and the GES v6 atomic line list (Heiter 2021), synthesized
    on the fly. Fitting regions (segments_file) are the wings of H-alpha, H-beta
    and the Mg I triplet (Teff, logg) plus Fe I / Fe II lines ([Fe/H], vsini).

    teff, logg, MH and vsini are always free; R is fixed at `resolution`
    (NEID is well known). alpha / vmic / vmac each take a mode:
      alpha_mode : 'free' | 'auto'   ('auto' ties [alpha/Fe] to [M/H])
      vmic_mode  : 'free' | 'GES'    (Doyle has no vmic relation)
      vmac_mode  : 'free' | 'GES' | 'Doyle2014'
    For 'GES'/'Doyle2014' the parameter is tied to teff/logg/MH and recomputed
    every iteration (matches the iSpec GUI 'Automatic from empirical relation');
    'Doyle2014' is enabled via a runtime patch, no iSpec source change.

    Defaults match the iSpec GUI settings used for the paper: free teff/logg/MH/
    alpha/vmic/vsini, vmac from the empirical relation (Doyle 2014), R fixed.

    use_errors=False (GUI 'Use errors for fitting' unchecked) minimizes the
    UNWEIGHTED sum((obs-model)^2); the printed CHI-SQUARE then matches the GUI.
    use_errors=True weights residuals by 1/sqrt(err) (per iSpec), which both
    changes the CHI-SQUARE scale (~1/<err> larger) and shrinks the formal
    parameter errors. Default False to reproduce the GUI workflow.

    Returns (params, errors, obs_spec, model_spec, status).
    """
    ispec = _ispec()
    if alpha_mode not in ('free', 'auto'):
        raise ValueError("alpha_mode must be 'free' or 'auto'")
    if vmic_mode not in ('free', 'GES'):
        raise ValueError("vmic_mode must be 'free' or 'GES' (no Doyle vmic relation)")
    if vmac_mode not in ('free', 'GES', 'Doyle2014'):
        raise ValueError("vmac_mode must be 'free', 'GES' or 'Doyle2014'")
    if resources is None:
        resources = load_rt_resources(spec)
    modeled_layers_pack, atomic_linelist, isotopes, solar_abundances = resources
    segments = ispec.read_segment_regions(segments_file)

    # build free_params + empirical-relation flags from the modes
    free_params = ['teff', 'logg', 'MH', 'vsini']
    enhance_abundances = (alpha_mode == 'auto')     # GUI 'Automatic alpha enhancement'
    vmic_from_rel = (vmic_mode == 'GES')
    vmac_from_rel = (vmac_mode in ('GES', 'Doyle2014'))
    if alpha_mode == 'free':
        free_params.append('alpha')
    if vmic_mode == 'free':
        free_params.append('vmic')
    if vmac_mode == 'free':
        free_params.append('vmac')
    _set_vmac_relation('Doyle2014' if vmac_mode == 'Doyle2014' else 'GES')

    cont = ispec.fit_continuum(spec, fixed_value=1.0, model='Fixed value')
    alpha = 0.0 if alpha_mode == 'free' else ispec.determine_abundance_enchancements(MH)
    vmic = ispec.estimate_vmic(teff, logg, MH)
    vmac = (ispec.estimate_vmac(teff, logg, MH, relation='Doyle2014')
            if vmac_mode == 'Doyle2014' else ispec.estimate_vmac(teff, logg, MH))

    obs, model, params, errors, abund, loggf, status, _ = ispec.model_spectrum(
        spec, cont, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances,
        free_abundances=None, linelist_free_loggf=None,
        initial_teff=teff, initial_logg=logg, initial_MH=MH,
        initial_alpha=alpha, initial_vmic=vmic, initial_vmac=vmac,
        initial_vsini=vsini, initial_limb_darkening_coeff=0.6,
        initial_R=resolution, initial_vrad=0.0,
        free_params=free_params, segments=segments, linemasks=None,
        enhance_abundances=enhance_abundances, scale=None,
        vmic_from_empirical_relation=vmic_from_rel,
        vmac_from_empirical_relation=vmac_from_rel,
        use_errors=use_errors, max_iterations=max_iterations,
        verbose=1 if verbose else 0,          # iSpec prints per-iteration params + chisq
        code=RT_CODE,
    )
    if verbose:
        print(f'  rchisq={status.get("rchisq", float("nan")):.3f}, '
              f'niter={status.get("niter", "?")}')
        for k in ['teff', 'logg', 'MH', 'vmic', 'vmac', 'vsini']:
            if k in params:
                print(f'  {k:6s} = {params[k]:9.3f} +/- {errors.get(k, float("nan")):.3f}')
    return params, errors, obs, model, status


def save_fit_html(obs, model, out_path, params=None, status=None,
                  segments_file=SEGMENTS_FILE, title='NEID + iSpec fit'):
    """Write a zoomable interactive HTML (Plotly) of the fit: observed vs best-fit
    model on top, residual (obs-model) below. The fit segments (mask) are shaded.
    Only pixels the model covers (flux>0, i.e. inside the segments) are drawn for
    the model/residual. WebGL traces handle the ~200k points.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    w = np.asarray(obs['waveobs']); fo = np.asarray(obs['flux'])
    fm = np.asarray(model['flux'])
    m = np.isfinite(fm) & (fm > 0)                 # model only exists in fit segments
    res = np.where(m, fo - fm, np.nan)
    rms = float(np.sqrt(np.nanmean(res[m] ** 2))) if m.any() else float('nan')

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.72, 0.28],
                        vertical_spacing=0.03)

    # shade the fit segments (mask). vrect shapes don't render over Scattergl's
    # WebGL canvas, so draw the bands as a filled Scattergl trace (same layer as
    # the data). Rectangles are built as one trace via 'toself' + None separators;
    # y extents overshoot each panel so the bands stay full-height under zoom.
    if segments_file and os.path.exists(segments_file):
        seg = np.genfromtxt(segments_file, names=True)
        bs = np.atleast_1d(seg['wave_base']); ts = np.atleast_1d(seg['wave_top'])

        def _bands(y0, y1, row, showlegend):
            xs, ys = [], []
            for b, t in zip(bs, ts):
                xs += [b, b, t, t, None]
                ys += [y0, y1, y1, y0, None]
            fig.add_trace(go.Scattergl(
                x=xs, y=ys, fill='toself', mode='lines',
                fillcolor='rgba(70,130,180,0.15)', line=dict(width=0),
                name='fit region (mask)', legendgroup='mask',
                showlegend=showlegend, hoverinfo='skip'), row=row, col=1)

        _bands(-1.0, 3.0, 1, True)      # top panel (flux ~0.1-1.2)
        _bands(-5.0, 5.0, 2, False)     # residual panel

    fig.add_trace(go.Scattergl(x=w, y=fo, name='observed', mode='lines',
                               line=dict(color='black', width=1)), row=1, col=1)
    fig.add_trace(go.Scattergl(x=w[m], y=fm[m], name='best-fit model', mode='lines',
                               line=dict(color='crimson', width=1)), row=1, col=1)
    fig.add_trace(go.Scattergl(x=w[m], y=res[m], name='obs - model', mode='lines',
                               line=dict(color='seagreen', width=1)), row=2, col=1)
    fig.add_hline(y=1.0, line=dict(color='gray', dash='dash', width=0.7), row=1, col=1)
    for y in (0.0, rms, -rms):
        fig.add_hline(y=y, line=dict(color='gray', dash='dot', width=0.6), row=2, col=1)

    sub = ''
    if params:
        sub = '  |  ' + '  '.join(
            f'{k}={params[k]:.2f}' for k in ['teff', 'logg', 'MH', 'vmic', 'vmac', 'vsini']
            if k in params)
    if status and 'rms' in status:
        sub += f'  |  rms={status["rms"]:.4f}'
    fig.update_layout(title=title + sub, template='plotly_white', hovermode='x unified',
                      legend=dict(orientation='h', y=1.02, x=0), height=720,
                      margin=dict(l=60, r=20, t=60, b=50))
    fig.update_yaxes(title_text='normalized flux', range=[0, 1.2], row=1, col=1)
    fig.update_yaxes(title_text='residual', range=[-0.5, 0.5], row=2, col=1)
    fig.update_xaxes(title_text='wavelength [nm, air, rest]', row=2, col=1)
    fig.write_html(out_path, include_plotlyjs=True, full_html=True)
    return out_path
