"""Packaging bounds are coherent (audit finding 34): the lower bounds of setup.py can be pinned,
the oldest supported Python is the one the docs and the CI matrix state, and the unused
SigProfilerExtractor (it drags in torch) is not required. The resolution of the minimum
environment itself needs the network and runs in CI (.github/scripts/minimum_environment.py).
Run with: pytest tests/test_minimum_environment.py -v
"""

import importlib.util
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / ".github" / "scripts" / "minimum_environment.py"

spec = importlib.util.spec_from_file_location("minimum_environment", str(SCRIPT))
minimum_environment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(minimum_environment)

DECLARED = minimum_environment.declared_requirements()
FLOOR = minimum_environment.python_floor(DECLARED["python_requires"])


class TestScriptHelpers:
    def test_python_floor(self):
        assert minimum_environment.python_floor(">=3.11") == "3.11"
        assert minimum_environment.python_floor(">= 3.9.2") == "3.9"

    def test_python_floor_rejects_other_specifiers(self):
        with pytest.raises(ValueError):
            minimum_environment.python_floor(">=3.9,<4")

    def test_lower_bounds_become_pins(self):
        assert minimum_environment.minimum_pins(["numpy>=2.0.0", "seaborn >= 0.13.2"]) == [
            "numpy==2.0.0",
            "seaborn==0.13.2",
        ]

    def test_a_requirement_without_a_lower_bound_is_an_error(self):
        with pytest.raises(ValueError):
            minimum_environment.minimum_pins(["numpy"])


class TestDeclaredBounds:
    def test_every_requirement_has_a_lower_bound(self):
        assert len(minimum_environment.minimum_pins(DECLARED["install_requires"])) == len(
            DECLARED["install_requires"]
        )

    def test_numpy_two_is_required(self):
        assert "numpy==2.0.0" in minimum_environment.minimum_pins(DECLARED["install_requires"])

    def test_unused_signature_extractor_is_not_required(self):
        assert not any("extractor" in r.lower() for r in DECLARED["install_requires"])

    def test_the_package_does_not_import_the_extractor(self):
        for path in (ROOT / "SigProfilerClusters").rglob("*.py"):
            code = [l for l in path.read_text().splitlines() if not l.lstrip().startswith("#")]
            assert not any("SigProfilerExtractor" in l for l in code), str(path)


class TestSupportedPythonIsConsistent:
    def test_floor_is_at_least_3_11(self):
        assert tuple(int(x) for x in FLOOR.split(".")) >= (3, 11)

    def test_ci_matrix_starts_at_the_floor(self):
        ci = (ROOT / ".github" / "workflows" / "ci.yml").read_text()
        versions = re.search(r"python-version: \[([^\]]*)\]", ci).group(1)
        listed = [v.strip().strip("'\"") for v in versions.split(",")]
        assert listed[0] == FLOOR

    def test_ci_resolves_the_minimum_environment(self):
        ci = (ROOT / ".github" / "workflows" / "ci.yml").read_text()
        assert ".github/scripts/minimum_environment.py" in ci

    def test_docs_state_the_same_floor(self):
        docs = (ROOT / "docs" / "1_installation.md").read_text()
        assert re.search(r"version >= " + re.escape(FLOOR), docs)
