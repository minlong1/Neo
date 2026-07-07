"""
Author: Alaina Humiston
email: athompson9@hawk.iit.edu

Periodic Table of information into the BE and spin-orbit splitting/branching
ratios of all the elements, which the user will be able to select.
Goal is to create less user error by having default values for the GA to use
as an initial guess.

User simply picks the element and the algorithm examines the data and sets the
bounds accordingly.

Values for BE and spin-orbit splitting taken from phi XPS data handbook:
Jill Chastain and Roger C King Jr. "Handbook of X-ray photoelectron
spectroscopy". In: Perkin-Elmer Corporation 40 (1992), ISBN: 0-9627026-2-5
Values for BE_lit, width, and rec_width taken from
http://dx.doi.org/10.1006/adnd.2000.0848
Values for BE_alt and alt_width taken from
http://dx.doi.org/10.1007/978-3-540-28619-6 This list also includes the data
from the first, making it a greater culmination of data.
All values take the average of the BE found and the peak width.

Will want to add in data from XPS Oasis.
Need to also add in information about typical background subtraction/peak
shape to choose for each element i.e. Double Lorentzian for asymmetric ones
and such.

Note: This body of information is only as accurate as its sources. Recent
studies into x-ray data may be more applicable than the information presented
here.

BE are given from the pure elemental state.

The reference values themselves live in xps_neo/periodic_table_data.py
(generated; one entry per (element, photoelectron line), sparse against the
LineParams defaults below). This module provides the lookup:

- get_line_params(element, line) -> LineParams: the clean API.
- ElementData: backward-compatible adapter with the original getParams
  signature (10-slot lists in, 12-tuple of 10-slot lists out), used by the
  GUI's periodic table tab.
"""

# Issue: This is in BE not KE
# Should we store information about KLL edges and satellites?

from dataclasses import dataclass
from typing import Optional

from PhysicsModules.XPS.xps_neo.periodic_table_data import ELEMENT_LINES


@dataclass(frozen=True)
class LineParams:
    """Literature parameters for one (element, photoelectron line).

    Defaults here define what a sparse ELEMENT_LINES entry omits; an unknown
    element yields exactly these defaults (be_lit 0.0, Voigt singlet).
    """

    be_lit: float = 0.0
    be_alt: float = 0.0
    width: float = 0.05
    rec_width: float = 0.0
    alt_width: float = 0.0
    width_range: float = 0.05
    is_singlet: bool = True
    so_split: float = 0.0
    branching_ratio: float = 0.0
    peak_type: str = "Voigt"
    is_coster_kronig: bool = False
    # Coster-Kronig width from the literature; recorded in the data table but
    # not returned by the ElementData.getParams compat API (the original
    # implementation assigned it to a dead local).
    ck_width: Optional[float] = None


def get_line_params(element: str, line: str) -> LineParams:
    """Look up literature values for an element's photoelectron line.

    Falls back to the element's default line (its (element, None) entry,
    typically the strongest line) when the requested line is not listed, and
    to all-defaults when the element is unknown.
    """
    entry = ELEMENT_LINES.get((element, line))
    if entry is None:
        entry = ELEMENT_LINES.get((element, None), {})
    return LineParams(**entry)


class ElementData:  # Will need to change these values for KE input
    """Backward-compatible adapter for the GUI's periodic table tab."""

    def __init__(self, element, photoelectronLine):
        self.element = element
        self.photoelectronLine = photoelectronLine
        # Should Hydrogen be an option? Some XPS data can see He so I guess
        # that should stay

    def getParams(self, element, photoelectronLine):
        """Original 10-slot list API: returns (BE_lit, is_singlet, so_split,
        branching_ratio, BE_alt, alt_width, width_range, width, rec_width,
        Default, peakTypes, is_coster_kronig)."""
        params = [get_line_params(element[i], photoelectronLine[i]) for i in range(10)]
        return (
            [p.be_lit for p in params],
            [p.is_singlet for p in params],
            [p.so_split for p in params],
            [p.branching_ratio for p in params],
            [p.be_alt for p in params],
            [p.alt_width for p in params],
            [p.width_range for p in params],
            [p.width for p in params],
            [p.rec_width for p in params],
            False,  # Default: never set by the original implementation
            [p.peak_type for p in params],
            [p.is_coster_kronig for p in params],
        )
