"""
Per-band parameter search grids.

Raman analogue of the EXAFS module's ``pathrange.py``. Each band has three
searchable per-band parameters (amplitude, width, centre shift); the global
wavenumber shift (the E0 analogue) is one gene shared by all bands, with its
own grid (``DEFAULT_SHIFT_RANGE``). problem.py turns these grids into the
flat-genome ``ParameterSpace``.

Custom range file format (one band per row, 10 comma-separated columns)::

    band, amp_low, amp_high, amp_dt, width_low, width_high, width_dt,
          cshift_low, cshift_high, cshift_dt

Ported from the standalone ``raman_neo`` package (bandrange.py).
"""

import numpy as np

# Default (low, high, step) grids per band parameter
DEFAULT_AMP_RANGE = (0.02, 2.0, 0.02)  # a leading 0 is prepended (band "off")
DEFAULT_WIDTH_RANGE = (2.0, 40.0, 1.0)
DEFAULT_CSHIFT_RANGE = (-15.0, 15.1, 1.0)
# Global shift (E0 analogue): searched as a gene and swept at mid/end of run
DEFAULT_SHIFT_RANGE = (-10.0, 10.05, 0.1)


def check_range(val_low, val_high, delta_val):
    """(low, high, step) -> np.arange grid, rejecting a step wider than the span."""
    if delta_val <= 0 or np.abs(val_high - val_low) / delta_val < 1:
        raise ValueError(f"Delta spacing too high for range ({val_low}, {val_high}, {delta_val})")
    return np.arange(val_low, val_high, delta_val)


def shift_grid(shift_range=()):
    """Grid for the shared global-shift gene (default ±10 cm^-1, 0.1 step)."""
    if len(shift_range) == 0:
        shift_range = DEFAULT_SHIFT_RANGE
    return check_range(*shift_range)


class Bandrange_limits:
    """Search grids for one band's (amp, width, center_shift) genes, either
    the defaults or one row of a custom range file."""

    def __init__(self, _band_, custom_range=None):
        self._band_ = _band_

        if custom_range is None:
            self.Amp_min, self.Amp_max, self.Amp_dt = DEFAULT_AMP_RANGE
            self.Width_min, self.Width_max, self.Width_dt = DEFAULT_WIDTH_RANGE
            self.Cshift_min, self.Cshift_max, self.Cshift_dt = DEFAULT_CSHIFT_RANGE

            self.rangeAmp = np.insert(np.arange(*DEFAULT_AMP_RANGE), 0, 0.0)
            self.Amp_min = 0.0
        else:
            self.Amp_min, self.Amp_max, self.Amp_dt = custom_range[1:4]
            self.Width_min, self.Width_max, self.Width_dt = custom_range[4:7]
            self.Cshift_min, self.Cshift_max, self.Cshift_dt = custom_range[7:10]

            self.rangeAmp = check_range(self.Amp_min, self.Amp_max, self.Amp_dt)

        self.rangeWidth = check_range(self.Width_min, self.Width_max, self.Width_dt)
        self.rangeCshift = check_range(self.Cshift_min, self.Cshift_max, self.Cshift_dt)

    def get_lim(self):
        return (
            self._band_,
            self.get_lim_Amp(),
            self.get_lim_Width(),
            self.get_lim_Cshift(),
        )

    def get_lim_Amp(self):
        return self.Amp_min, self.Amp_max, self.Amp_dt

    def get_lim_Width(self):
        return self.Width_min, self.Width_max, self.Width_dt

    def get_lim_Cshift(self):
        return self.Cshift_min, self.Cshift_max, self.Cshift_dt

    def get_band(self):
        return self._band_

    def getrange_Amp(self):
        return self.rangeAmp

    def getrange_Width(self):
        return self.rangeWidth

    def getrange_Cshift(self):
        return self.rangeCshift

    def getrange(self):
        return self.getrange_Amp(), self.getrange_Width(), self.getrange_Cshift()


def read_bandrange_file(file, num_band):
    """Read a custom per-band range file (10 columns per band)."""
    raw_data = np.genfromtxt(file, delimiter=",", comments="#")
    if raw_data.ndim == 1:
        raw_data = raw_data.reshape(1, -1)
    if raw_data.shape[0] != num_band:
        raise ValueError(
            f"Mismatch in number of bands: range file has {raw_data.shape[0]}, fit has {num_band}"
        )
    if raw_data.shape[1] != 10:
        raise ValueError(f"Range file needs 10 columns per band, got {raw_data.shape[1]}")
    for row in raw_data:
        check_range(row[1], row[2], row[3])
        check_range(row[4], row[5], row[6])
        check_range(row[7], row[8], row[9])
    return raw_data


def build_bandrange_pars(num_band, bandrange_file=None):
    """One Bandrange_limits per fitted band: defaults, or rows of a range file."""
    if bandrange_file is None:
        return [Bandrange_limits(i) for i in range(num_band)]
    rows = read_bandrange_file(bandrange_file, num_band)
    return [Bandrange_limits(i, rows[i, :]) for i in range(num_band)]
