"""The default values documented for analysis() must match the code (audit finding #37).
The documented table is parsed from docs/4_using_the_tool_input.md and compared with the
signature of SigProfilerClusters.analysis and with the command line parser.
Run with: pytest tests/test_docs_defaults.py -v
"""

import inspect
import re
from pathlib import Path

import pytest

from SigProfilerClusters.SigProfilerClusters import analysis
from SigProfilerClusters.controllers.cli_controller import parse_arguments_clusters

DOCS = Path(__file__).resolve().parent.parent / "docs" / "4_using_the_tool_input.md"
ROW = re.compile(r"^\|\s*`(?P<name>\w+)`\s*\|.*?Default:\s*`(?P<default>[^`]+)`", re.M)


def _documented_defaults():
    return {m.group("name"): m.group("default") for m in ROW.finditer(DOCS.read_text())}


def _literal(text):
    return {"True": True, "False": False}.get(text, text.strip('"'))


DOCUMENTED = _documented_defaults()
SIGNATURE = {
    name: p.default
    for name, p in inspect.signature(analysis).parameters.items()
    if p.default is not inspect.Parameter.empty
}
CLI = vars(parse_arguments_clusters(["p", "GRCh37", "96", "96", "/d/"], "d"))


def test_the_table_was_parsed():
    assert len(DOCUMENTED) >= 10


@pytest.mark.parametrize("name", sorted(DOCUMENTED))
def test_documented_default_matches_the_python_default(name):
    assert name in SIGNATURE, f"{name} is documented but analysis() has no such parameter"
    documented = _literal(DOCUMENTED[name])
    if isinstance(SIGNATURE[name], (int, float)) and not isinstance(SIGNATURE[name], bool):
        documented = type(SIGNATURE[name])(documented)
    assert SIGNATURE[name] == documented


@pytest.mark.parametrize("name", sorted(set(DOCUMENTED) & set(CLI)))
def test_documented_default_matches_the_command_line_default(name):
    documented = _literal(DOCUMENTED[name])
    if isinstance(CLI[name], (int, float)) and not isinstance(CLI[name], bool):
        documented = type(CLI[name])(documented)
    assert CLI[name] == documented
