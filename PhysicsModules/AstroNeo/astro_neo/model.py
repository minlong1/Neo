"""
XSPEC model definitions for AstroNeo, built around a `ModelSpec`: an XSPEC
model expression plus its free parameters (fit genome), fixed/frozen
overrides, an explicit-unfreeze list, and parameter links -- all addressed
by dotted `component.parameter` name (e.g. `"powerlaw.PhoIndex"`) rather
than XSPEC's raw positional parameter index.

Names, not indices, because indices are only meaningful for one specific
model expression (renumber if the expression changes at all) and can only
be discovered by building the model and reading `model.show()` -- fragile
and opaque for anyone defining a new spec. Every `xspec.Parameter` exposes
its own `.index`, so `build_xspec_model`/`resolve_free_param_indices`
resolve names to indices once, up front; the actual per-generation fitness
evaluation still uses the resolved indices through `Model.setPars({index:
value})`, XSPEC's fast bulk setter -- resolving by name on every call
instead measured ~14x slower in practice (bulk index-based setPars: ~0.4s/
eval; per-parameter name-based attribute sets: ~5.5s/eval), which would
undo the whole point of building the Model once per problem instance (see
`problem.py`).

`DEFAULT_MODEL_SPEC` is the one X-ray CCD spectral model originally ported
from Astro_Neo's `fitness.py` -- an absorbed powerlaw plus a two-temperature
APEC plasma and an ACX2 charge-exchange component -- with the exact same
expression, fixed/frozen parameter values, links, and free-parameter set,
just re-expressed by name instead of by the original's raw indices.
Fitting a *different* XSPEC model means constructing a different
`ModelSpec` (see the module README's "Generalizing to a different model")
and passing it to `AstroNeoProblem(model_spec=...)`; nothing else in this
module needs to change.
"""

from dataclasses import dataclass, field
from typing import Dict, Tuple


def _get_param(model, dotted_name: str):
    """Resolve 'component.parameter' to its xspec.Parameter object."""
    *components, leaf = dotted_name.split(".")
    obj = model
    for name in components:
        obj = getattr(obj, name)
    return getattr(obj, leaf)


def _set_param(model, dotted_name: str, value: float) -> None:
    *components, leaf = dotted_name.split(".")
    obj = model
    for name in components:
        obj = getattr(obj, name)
    setattr(obj, leaf, value)


@dataclass(frozen=True)
class ModelSpec:
    """Declares one fittable XSPEC model.

    `expr`: an XSPEC model expression string, e.g. "TBabs*powerlaw".
    `free_params`: dotted `component.parameter` names, in genome order.
    `param_ranges`: (low, high) fit bounds per name in `free_params`.
    `fixed_params`: dotted name -> value, applied once at model build time
        (before `unfrozen`/`links`, so a link or explicit unfreeze below
        can override a value/frozen-state set here -- same precedence as
        the original hand-written model construction this replaces).
    `unfrozen`: dotted names to explicitly set `.frozen = False` on after
        `fixed_params` is applied (XSPEC freezes some repeated-component
        parameters by default; this is how the original model un-freezes
        its second APEC temperature).
    `links`: (target, source) dotted-name pairs; `target.link = source`,
        applied last (so a linked parameter's value/frozen state always
        wins over anything set for it in `fixed_params`).
    """

    expr: str
    free_params: Tuple[str, ...]
    param_ranges: Dict[str, Tuple[float, float]]
    fixed_params: Dict[str, float] = field(default_factory=dict)
    unfrozen: Tuple[str, ...] = ()
    links: Tuple[Tuple[str, str], ...] = ()


DEFAULT_MODEL_SPEC = ModelSpec(
    expr="TBabs(TBabs*powerlaw + lsmooth(vapec + vapec + zashift*vacx2))",
    free_params=(
        "TBabs_2.nH",
        "powerlaw.PhoIndex",
        "powerlaw.norm",
        "vapec.kT",
        "vapec.C",
        "vapec.N",
        "vapec.O",
        "vapec.Ne",
        "vapec.Mg",
        "vapec.Fe",
        "vapec.norm",
        "vapec_6.kT",
        "vapec_6.norm",
        "vacx2.collnpar",
        "vacx2.norm",
    ),
    param_ranges={
        "TBabs_2.nH": (0.00, 2e-5),
        "powerlaw.PhoIndex": (0.05, 1.1),
        "powerlaw.norm": (1e-4, 1e-3),
        "vapec.kT": (0.0808, 1),
        "vapec.C": (0.00, 1.00),
        "vapec.N": (0.00, 1.50),
        "vapec.O": (0.0, 1.00),
        "vapec.Ne": (0.0, 1.00),
        "vapec.Mg": (0.0, 2.00),
        "vapec.Fe": (0.0, 1.00),
        "vapec.norm": (0.0001, 0.0010),
        "vapec_6.kT": (0.0808, 0.75),
        "vapec_6.norm": (0, 0.0010),
        "vacx2.collnpar": (250, 500),
        "vacx2.norm": (2e-4, 5e-4),
    },
    fixed_params={
        "TBabs.nH": 2.79000e-02,
        "TBabs_2.nH": 1.87763e-05,
        "powerlaw.PhoIndex": 1.00455,
        "powerlaw.norm": 5.94711e-04,
        "lsmooth.Sig_6keV": 2.00000e-02,
        "lsmooth.Index": 1.00000,
        "vapec.kT": 0.788251,
        "vapec.He": 1.00000,
        "vapec.C": 0.644919,
        "vapec.N": 1.07904,
        "vapec.O": 0.205951,
        "vapec.Ne": 0.487054,
        "vapec.Mg": 1.35956,
        "vapec.Al": 1.00000,
        "vapec.Si": 1.00000,
        "vapec.S": 1.00000,
        "vapec.Ar": 1.00000,
        "vapec.Ca": 1.00000,
        "vapec.Fe": 0.188908,
        "vapec.Ni": 1.00000,
        "vapec.Redshift": 8.10000e-04,
        "vapec.norm": 3.41553e-04,
        "vapec_6.kT": 0.434492,
        "vapec_6.He": 1.00000,
        "vapec_6.C": 0.644919,
        "vapec_6.N": 1.07904,
        "vapec_6.O": 0.205951,
        "vapec_6.Ne": 0.487054,
        "vapec_6.Mg": 1.35956,
        "vapec_6.Al": 1.00000,
        "vapec_6.Si": 1.00000,
        "vapec_6.S": 1.00000,
        "vapec_6.Ar": 1.00000,
        "vapec_6.Ca": 1.00000,
        "vapec_6.Fe": 0.188908,
        "vapec_6.Ni": 1.00000,
        "vapec_6.Redshift": 8.10000e-04,
        "vapec_6.norm": 7.27081e-04,
        "zashift.Redshift": 8.10000e-04,
        "vacx2.temperature": 0.434492,
        "vacx2.collnpar": 272.624,
        "vacx2.collntype": 4.00000,
        "vacx2.acxmodel": 2.00000,
        "vacx2.recombtype": 2.00000,
        "vacx2.Hefrac": 9.00000e-02,
        "vacx2.H": 1.00000,
        "vacx2.He": 1.00000,
        "vacx2.C": 0.644919,
        "vacx2.N": 1.07904,
        "vacx2.O": 0.205951,
        "vacx2.Ne": 0.487054,
        "vacx2.Mg": 1.35956,
        "vacx2.Al": 1.00000,
        "vacx2.Si": 1.00000,
        "vacx2.S": 1.00000,
        "vacx2.Ar": 1.00000,
        "vacx2.Ca": 1.00000,
        "vacx2.Fe": 0.188908,
        "vacx2.Ni": 1.00000,
        "vacx2.norm": 2.43828e-04,
    },
    unfrozen=("vapec_6.kT",),
    links=(
        ("vapec_6.C", "vapec.C"),
        ("vapec_6.N", "vapec.N"),
        ("vapec_6.O", "vapec.O"),
        ("vapec_6.Ne", "vapec.Ne"),
        ("vapec_6.Mg", "vapec.Mg"),
        ("vapec_6.Fe", "vapec.Fe"),
        ("vapec_6.Redshift", "vapec.Redshift"),
        ("vacx2.temperature", "vapec_6.kT"),
        ("vacx2.C", "vapec.C"),
        ("vacx2.N", "vapec.N"),
        ("vacx2.O", "vapec.O"),
        ("vacx2.Ne", "vapec.Ne"),
        ("vacx2.Mg", "vapec.Mg"),
        ("vacx2.Fe", "vapec.Fe"),
    ),
)

# Backwards-compatible aliases: model.PARAM_NAMES/PARAM_RANGES used to be
# the only spec (dotted names now, previously flat e.g. "TBabs_2_nH").
PARAM_NAMES = DEFAULT_MODEL_SPEC.free_params
PARAM_RANGES = DEFAULT_MODEL_SPEC.param_ranges


def build_xspec_model(spec: ModelSpec = DEFAULT_MODEL_SPEC):
    """Build `spec.expr`, apply its fixed params, unfreeze list, and
    parameter links, in that order. Requires PyXspec; caller imports
    xspec (and registers any local models the expression needs, e.g.
    ACX2) first."""
    import xspec

    model = xspec.Model(spec.expr)
    for name, value in spec.fixed_params.items():
        _set_param(model, name, value)
    for name in spec.unfrozen:
        _get_param(model, name).frozen = False
    for target, source in spec.links:
        _get_param(model, target).link = _get_param(model, source)
    return model


def resolve_free_param_indices(model, spec: ModelSpec = DEFAULT_MODEL_SPEC):
    """XSPEC's 1-based parameter index for each of spec.free_params, in
    order -- resolved once against a built model, then reused every
    fitness() call for the fast bulk Model.setPars({index: value})."""
    return [_get_param(model, name).index for name in spec.free_params]


def genes_to_xspec_params(genes, free_param_indices):
    """Map an ordered genome (spec.free_params order) to the
    {xspec_index: value} dict Model.setPars() expects."""
    return dict(zip(free_param_indices, (float(g) for g in genes)))
