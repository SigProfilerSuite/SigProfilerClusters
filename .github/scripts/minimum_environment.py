"""Check that the lower bounds declared in setup.py form an installable environment.

Every `name>=x` of install_requires becomes `name==x` and pip resolves them for the oldest
supported Python (python_requires) with `pip install --dry-run`, using wheels only, so that
nothing is installed and no compiler is needed. A conflict between two lower bounds, or a
bound with no wheel for the oldest Python, fails the check (audit finding 34).

Run from anywhere: python .github/scripts/minimum_environment.py
"""

import ast
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def declared_requirements(setup_path=ROOT / "setup.py"):
    """python_requires and install_requires, read from setup.py without executing it."""
    tree = ast.parse(Path(setup_path).read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "setup":
            wanted = ("python_requires", "install_requires")
            return {kw.arg: ast.literal_eval(kw.value) for kw in node.keywords if kw.arg in wanted}
    raise ValueError("no setup() call found in " + str(setup_path))


def python_floor(python_requires):
    """'>=3.11' -> '3.11'."""
    match = re.fullmatch(r"\s*>=\s*(\d+\.\d+)(?:\.\d+)?\s*", python_requires)
    if match is None:
        raise ValueError("python_requires must be a plain '>=X.Y' bound, got %r" % python_requires)
    return match.group(1)


def minimum_pins(requirements):
    """Turn every 'name>=x' into 'name==x'; a requirement without a lower bound is an error."""
    pins = []
    for requirement in requirements:
        match = re.fullmatch(r"\s*([A-Za-z0-9_.-]+)\s*>=\s*([A-Za-z0-9_.]+)\s*", requirement)
        if match is None:
            raise ValueError("requirement without a plain '>=' lower bound: %r" % requirement)
        pins.append("%s==%s" % (match.group(1), match.group(2)))
    return pins


def main():
    declared = declared_requirements()
    floor = python_floor(declared["python_requires"])
    pins = minimum_pins(declared["install_requires"])
    print("oldest supported Python: " + floor)
    print("minimum environment: " + " ".join(pins))
    with tempfile.TemporaryDirectory() as target:
        command = [
            sys.executable, "-m", "pip", "install", "--dry-run", "--ignore-installed",
            "--python-version", floor, "--only-binary=:all:", "--target", target,
        ] + pins
        return subprocess.run(command).returncode


if __name__ == "__main__":
    sys.exit(main())
