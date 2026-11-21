"""Every module shipped in the SigProfilerClusters package must be importable.
A module that cannot be imported is dead code that still gets packaged.
No genome reference required.
Run with: pytest tests/test_package_imports.py -v
"""

import importlib
import pkgutil

import pytest

import SigProfilerClusters

MODULES = sorted(
    info.name
    for info in pkgutil.walk_packages(
        SigProfilerClusters.__path__, prefix="SigProfilerClusters."
    )
)


def test_package_has_modules():
    assert "SigProfilerClusters.SigProfilerClusters" in MODULES
    assert "SigProfilerClusters.hotspot" in MODULES


@pytest.mark.parametrize("module", MODULES)
def test_module_is_importable(module):
    importlib.import_module(module)
