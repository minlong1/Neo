"""Tolerance-aware comparison of canonical golden texts (design.md Phase 4,
order-changing tolerance class).

Splits both texts into tokens; numeric tokens must agree within
rtol/atol, everything else must match exactly. Reports the worst
relative deviation so a change can be classified against the policy
(rtol <= 1e-9) before rebaselining.
"""

import math
import re

NUM_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def _tokens(text):
    parts = []
    pos = 0
    for m in NUM_RE.finditer(text):
        if m.start() > pos:
            parts.append(("text", text[pos : m.start()]))
        parts.append(("num", m.group()))
        pos = m.end()
    if pos < len(text):
        parts.append(("text", text[pos:]))
    return parts


def compare(expected, actual, rtol=1e-9, atol=1e-12):
    """Return (ok, worst_rel, message). Numeric tokens compared within
    tolerance; positional/textual structure must match exactly."""
    exp_toks = _tokens(expected)
    act_toks = _tokens(actual)
    if len(exp_toks) != len(act_toks):
        return (
            False,
            math.inf,
            (
                f"token count differs: {len(exp_toks)} vs {len(act_toks)} "
                f"(structural change, not a rounding difference)"
            ),
        )
    worst = 0.0
    for (ek, ev), (ak, av) in zip(exp_toks, act_toks):
        if ek != ak:
            return False, math.inf, f"token kind mismatch: {ev!r} vs {av!r}"
        if ek == "text":
            if ev != av:
                return False, math.inf, f"text differs: {ev!r} vs {av!r}"
            continue
        e, a = float(ev), float(av)
        if e == a:
            continue
        err = abs(a - e)
        rel = err / max(abs(e), abs(a))
        if err > atol + rtol * max(abs(e), abs(a)):
            return (
                False,
                rel,
                (
                    f"numeric deviation beyond tolerance: {ev} vs {av} "
                    f"(rel {rel:.3e})"
                ),
            )
        worst = max(worst, rel)
    return True, worst, f"within tolerance (worst rel {worst:.3e})"
