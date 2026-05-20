"""Test package bootstrap.

Ensure the repository root is on ``sys.path`` so imports like
``PhysicsModules.EXAFS...`` resolve correctly when tests are discovered from
this subdirectory.
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
	sys.path.insert(0, str(ROOT))

