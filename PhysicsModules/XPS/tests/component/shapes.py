"""Fixed inputs for numerical component tests of peak/background math.

Everything here is the single source of truth shared by the tests and the
reference generator, so a reference file always matches the invocation that
produced it. Values are physically sane for the DeltaPux1000 U-4d-like
spectrum used across the test suite (BE grid 480 -> 410 descending).
"""

import numpy as np

# Descending BE grid, matching the repo's sample data convention.
X_GRID = np.arange(480.0, 410.0, -0.1)

# Common peak arguments: (BE, width/gamma, sigma, A, asym, asymD,
#                         singlet, coster_kronig, width_CK, asym_CK,
#                         branch, split)
PEAK_ARGS = dict(
    BE=421.2,
    width=0.5,
    sigma=1.0,
    A=1000.0,
    asym=1.0,
    asymD=0.2,
    singlet=True,
    coster_kronig=False,
    width_CK=0.6,
    asym_CK=1.1,
    branch=0.57,
    split=13.0,
)

# Per-peak 1-D [min, max, delta] ranges, as Individual hands them to `peak`.
# Every range must produce a non-empty np.arange (bug register #6).
PARAM_RANGE = {
    "Gaussian": [0.7, 1.3, 0.001],
    "Lorentzian": [0.3, 0.7, 0.001],
    "Lorentzian Coster-Kronig": [0.1, 1.5, 0.001],
    "Binding Energy": [-0.2, 0.2, 0.01],
    "Amplitude": [800.0, 1200.0, 0.05],
    "Asymmetry": [0.9, 1.1, 0.01],
    "Asymmetry Range Coster-Kronig": [0.0, 10.0, 0.01],
    "Asymmetry Doniach-Sunjic": [0.0, 1.0, 0.01],
    "spinOrbitSplitting": [-0.02, 0.02, 0.01],
    "br_range": [0.45, 0.55, 0.01],
    # background-only keys
    "k_range": [0.0, 0.1, 0.001],
    "Background": [0.0, 5000.0, 0.05],
    "Slope": [0.0, 200.0, 0.01],
    "baseline": [0.0, 100.0, 0.1],
    "BE": [421.2],
}

# name -> (peakType, arg overrides)
PEAK_VARIANTS = {
    "voigt_singlet": ("Voigt", {}),
    "voigt_doublet": ("Voigt", {"singlet": False}),
    "gaussian_singlet": ("Gaussian", {}),
    "lorentzian_singlet": ("Lorentzian", {}),
    "double_lorentzian_singlet": ("Double Lorentzian", {}),
    "double_lorentzian_doublet": ("Double Lorentzian", {"singlet": False}),
    "double_lorentzian_coster_kronig": (
        "Double Lorentzian",
        {"singlet": False, "coster_kronig": True},
    ),
    "doniach_sunjic_singlet": ("Doniach-Sunjic", {}),
}

# Deterministic state poked onto background instances after construction,
# replacing the values their __init__ draws from the RNG.
BKGN_STATE = {"k": 0.05, "baseline_value": 50.0, "background": 100.0, "slope": 1.2}

# name -> background type string (as accepted by background.__init__)
BKGN_VARIANTS = {
    "baseline": "Baseline",
    "linear": "Linear",
    "shirley": "Shirley",
    "svsc": "SVSC",
}

# A fixed synthetic "measured" y for background functions that need one:
# one broad Voigt-ish bump on a sloped base, deterministic closed form.
Y_MEASURED = (
    50.0
    + 0.5 * (X_GRID - 410.0)
    + 1000.0 * np.exp(-((X_GRID - 421.2) ** 2) / (2 * 1.2**2))
)

# What Individual.getFit passes as `tot_area`: the summed peak-fit *array*
# (see xps_individual.py:719/818), not a scalar area despite the name.
TOT_PEAK_FIT = 1000.0 * np.exp(-((X_GRID - 421.2) ** 2) / (2 * 1.2**2))


def make_peak(xps_fit, name, seed=1234):
    """Instantiate `peak` for a variant; RNG seeded so init draws are fixed
    (they are overwritten by explicit args when the shape function runs)."""
    import random

    peak_type, overrides = PEAK_VARIANTS[name]
    random.seed(seed)
    np.random.seed(seed)
    p = xps_fit.peak(PARAM_RANGE, peak_type)
    args = dict(PEAK_ARGS, **overrides)
    return p, args


def eval_peak(xps_fit, name):
    p, a = make_peak(xps_fit, name)
    y = p.getY(
        X_GRID,
        a["BE"],
        a["width"],
        a["sigma"],
        a["A"],
        a["asym"],
        a["asymD"],
        a["singlet"],
        a["coster_kronig"],
        a["width_CK"],
        a["asym_CK"],
        a["branch"],
        a["split"],
    )
    return np.asarray(y, dtype=np.float64)


def eval_background(xps_fit, name, seed=1234):
    import random

    random.seed(seed)
    np.random.seed(seed)
    b = xps_fit.background(PARAM_RANGE, BKGN_VARIANTS[name], "Voigt")
    for attr, val in BKGN_STATE.items():
        if hasattr(b, attr):
            setattr(b, attr, val)
    a = PEAK_ARGS
    y = b.getY(
        X_GRID,
        Y_MEASURED,
        a["BE"],
        a["width"],
        a["sigma"],
        a["A"],
        a["asym"],
        a["asymD"],
        a["singlet"],
        a["coster_kronig"],
        a["width_CK"],
        a["asym_CK"],
        a["branch"],
        a["split"],
        False,
        TOT_PEAK_FIT,
    )
    return np.asarray(y, dtype=np.float64)
