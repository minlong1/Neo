# Known-broken cases

These INIs exercise features that are **broken in the current CLI** (verified
2026-07-04 on the pinned baseline environment). They are excluded from the
golden matrix until the underlying bugs are fixed — at which point they move
back to `cases/` and get captured with a CHANGELOG entry.

| Case | Failure |
|---|---|
| `voigt_tougaard2.ini`, `voigt_tougaard3.ini` | numba `TypingError`: `background.peak_shirley` is decorated `@nb.jit` but takes `self` (a `background` instance), which numba cannot type. Any Tougaard background crashes in `get_Background` (`xps_neo/xps_fit.py:1448` → `:2011`). |
| `voigt_peak_adding.ini` | `XPS_GA.addPeak` appends to every per-peak list except `asym_limited` (and the corresponding asym bookkeeping), so the rebuild loop crashes with `IndexError` at `xps_neo/xps.py:768` the first time a peak is added. Peak addition is unusable from the CLI today. |

See the bug register in `design.md` §7. Do not fix these during phases 0–1;
characterize-as-is comes first.
