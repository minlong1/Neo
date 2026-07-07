"""Shared test setup for XPS_Neo.

Since Phase 2 the package imports cleanly with no CLI arguments;
`load_xps_module` remains as the single import point used by the component
tests and the reference generator (it was an argv-injection shim before
Phase 2 and keeps its name so callers didn't change).
"""

import pathlib

TESTS_DIR = pathlib.Path(__file__).parent


def load_xps_module(name):
    """Import and return a xps_neo submodule."""
    import importlib

    return importlib.import_module(name)


def load_case_config(case):
    """Parse a golden case INI into a config dict via the real pipeline."""
    from PhysicsModules.XPS.xps_neo.parser2 import read_input_file
    from PhysicsModules.XPS.xps_neo.ini_parser import load_config

    ini = TESTS_DIR / "golden" / "cases" / f"{case}.ini"
    return load_config(read_input_file(str(ini)))
