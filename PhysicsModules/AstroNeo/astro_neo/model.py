"""
The one X-ray CCD spectral model this module fits: an absorbed powerlaw
plus a two-temperature APEC plasma and an ACX2 charge-exchange component,
ported unchanged (same expression, same fixed/frozen parameter values,
same free-parameter set) from the original Astro_Neo `fitness.py`.

Unlike EXAFS (any set of FEFF paths) or NanoIndentation (any number of
Oliver-Pharr terms), this model expression is not a runtime choice -- it is
literally hardcoded, same as upstream. Generalizing to an arbitrary XSPEC
model string is future work, not something this port invents; see the
module README's "Known limitations".

XSPEC has no by-name parameter setter for Model.setPars(dict) -- only
1-based positional indices into the model expression below. FREE_PARAM_INDICES
and the keys of FIXED_PARAMS are therefore tied exactly to MODEL_EXPR; if the
expression ever changes, both must be regenerated (e.g. via `model.show()`).
"""

MODEL_EXPR = "TBabs(TBabs*powerlaw + lsmooth(vapec + vapec + zashift*vacx2))"

# Free parameters this module fits, in genome order.
PARAM_NAMES = [
    "TBabs_2_nH",
    "PhoIndex",
    "Pl_norm",
    "vapec_kT",
    "vapec_C",
    "vapec_N",
    "vapec_O",
    "vapec_Ne",
    "vapec_Mg",
    "vapec_Fe",
    "vapec_norm",
    "vapec_6_kT",
    "vapec_6_norm",
    "vacx2_collnpar",
    "vacx2_norm",
]

# (low, high) fit bounds per free parameter.
PARAM_RANGES = {
    "TBabs_2_nH": (0.00, 2e-5),
    "PhoIndex": (0.05, 1.1),
    "Pl_norm": (1e-4, 1e-3),
    "vapec_kT": (0.0808, 1),
    "vapec_C": (0.00, 1.00),
    "vapec_N": (0.00, 1.50),
    "vapec_O": (0.0, 1.00),
    "vapec_Ne": (0.0, 1.00),
    "vapec_Mg": (0.0, 2.00),
    "vapec_Fe": (0.0, 1.00),
    "vapec_norm": (0.0001, 0.0010),
    "vapec_6_kT": (0.0808, 0.75),
    "vapec_6_norm": (0, 0.0010),
    "vacx2_collnpar": (250, 500),
    "vacx2_norm": (2e-4, 5e-4),
}

# XSPEC's 1-based parameter index for each name above, in the same order --
# this is what Model.setPars({index: value}) needs.
FREE_PARAM_INDICES = [2, 3, 4, 7, 9, 10, 11, 12, 13, 19, 22, 23, 38, 41, 60]

# Every other (frozen/linked) parameter in MODEL_EXPR, by index, ported
# unchanged from the original fit.
FIXED_PARAMS = {
    1: 2.79000e-02,
    2: 1.87763e-05,
    3: 1.00455,
    4: 5.94711e-04,
    5: 2.00000e-02,
    6: 1.00000,
    7: 0.788251,
    8: 1.00000,
    9: 0.644919,
    10: 1.07904,
    11: 0.205951,
    12: 0.487054,
    13: 1.35956,
    14: 1.00000,
    15: 1.00000,
    16: 1.00000,
    17: 1.00000,
    18: 1.00000,
    19: 0.188908,
    20: 1.00000,
    21: 8.10000e-04,
    22: 3.41553e-04,
    23: 0.434492,
    24: 1.00000,
    25: 0.644919,
    26: 1.07904,
    27: 0.205951,
    28: 0.487054,
    29: 1.35956,
    30: 1.00000,
    31: 1.00000,
    32: 1.00000,
    33: 1.00000,
    34: 1.00000,
    35: 0.188908,
    36: 1.00000,
    37: 8.10000e-04,
    38: 7.27081e-04,
    39: 8.10000e-04,
    40: 0.434492,
    41: 272.624,
    42: 4.00000,
    43: 2.00000,
    44: 2.00000,
    45: 9.00000e-02,
    46: 1.00000,
    47: 1.00000,
    48: 0.644919,
    49: 1.07904,
    50: 0.205951,
    51: 0.487054,
    52: 1.35956,
    53: 1.00000,
    54: 1.00000,
    55: 1.00000,
    56: 1.00000,
    57: 1.00000,
    58: 0.188908,
    59: 1.00000,
    60: 2.43828e-04,
}


def build_xspec_model():
    """Build MODEL_EXPR with FIXED_PARAMS and the vapec_6/vacx2 parameter
    links set up. Requires PyXspec; caller imports xspec first."""
    import xspec

    model = xspec.Model(MODEL_EXPR, setPars=FIXED_PARAMS)

    model.vapec_6.kT.frozen = False
    model.vapec_6.C.link = model.vapec.C
    model.vapec_6.N.link = model.vapec.N
    model.vapec_6.O.link = model.vapec.O
    model.vapec_6.Ne.link = model.vapec.Ne
    model.vapec_6.Mg.link = model.vapec.Mg
    model.vapec_6.Fe.link = model.vapec.Fe
    model.vapec_6.Redshift.link = model.vapec.Redshift

    model.vacx2.temperature.link = model.vapec_6.kT
    model.vacx2.C.link = model.vapec.C
    model.vacx2.N.link = model.vapec.N
    model.vacx2.O.link = model.vapec.O
    model.vacx2.Ne.link = model.vapec.Ne
    model.vacx2.Mg.link = model.vapec.Mg
    model.vacx2.Fe.link = model.vapec.Fe

    return model


def genes_to_xspec_params(genes):
    """Map an ordered genome (PARAM_NAMES order) to the {xspec_index: value}
    dict Model.setPars() expects."""
    return dict(zip(FREE_PARAM_INDICES, (float(g) for g in genes)))
