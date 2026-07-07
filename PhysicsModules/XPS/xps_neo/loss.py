"""The canonical XPS_Neo loss function (Phase 3 unification).

Extracted op-for-op from XPS_GA.fitness so the GA (xps_neo/xps.py) and the
GUI post-analysis (gui/xps_analysis2.py) score candidates with the SAME
math. The CLI definition is canonical (design.md Phase 3; decision log in
docs/fork-divergence.md): background window N = 50, and penalties apply
only while apply_penalties is True (the GA passes genNum < gen_alt; the
post-analysis passes False to report the loss the GA optimizes in its
final phase).

Phase 4.1 (order-changing tolerance class, rtol-verified, user-approved):
the per-point Python loop is vectorized. Summation order changes, so
values may differ from the loop version at the ~1e-15 relative level.
Dead computations the loop carried (peak_diff, overfit_penalty,
y_low_BE_penalty — computed but never added to the loss) are dropped;
this also removes the loop's potential UnboundLocalError on spectra that
never rise in scan direction (bug register #1, fixed here).
"""

import numpy as np


def compute_loss(x_array, y_array, yTotal, y_left_avg, y_right_avg, apply_penalties):
    """Return (loss, residual) for a candidate curve yTotal against data.

    x_array, y_array: measured spectrum (x may be BE-descending or
        KE-ascending; only the y edge averages depend on direction and
        the caller supplies those).
    yTotal: candidate fit evaluated on x_array.
    y_left_avg, y_right_avg: data edge averages, computed by the caller
        the same way as before (N=30 windows in both callers).
    apply_penalties: include the envelope/background penalty terms
        (the GA's early-generation mode, genNum < gen_alt).
    """
    yTotal = np.asarray(yTotal, dtype=np.float64)
    y_arr = np.asarray(y_array, dtype=np.float64)

    difference = yTotal - y_arr
    residual = difference.copy()
    sigma = np.sqrt(np.abs(y_arr))

    if not apply_penalties:
        # After gen_alt generations the loss focuses only on the fit
        # quality; the population should already be culled by then.
        loss = np.sum(difference**2 * sigma)
        return float(loss), residual

    N = 50  # background window; changed from 30 (see fork-divergence.md)
    if y_left_avg > y_right_avg:  # BE condition
        window_avg = np.sum(yTotal[:N]) / N
        edge_avg = y_left_avg
    else:  # KE
        window_avg = np.sum(yTotal[-N:]) / N
        edge_avg = y_right_avg
    # Penalty for the background window sitting above OR below the data
    # edge (the loop's two mutually-exclusive branches used the same
    # magnitude in both directions).
    bkgn_diff = abs(window_avg - edge_avg)

    # Envelope penalties: overshoot beyond 0.5*sigma, and any undershoot
    # (difference < 0 always satisfies the loop's `<= 0.5*sigma` check).
    penalty = np.where(
        (difference > 0) & (difference >= 0.5 * sigma), 1000.0 * sigma, 0.0
    )
    penalty_under = np.where(difference < 0, 500.0 * sigma, 0.0)

    loss = np.sum(
        difference**2 * sigma
        + penalty * sigma
        + penalty_under * sigma
        + 3000.0 * sigma * bkgn_diff
    )
    return float(loss), residual
