"""CI and test infrastructure consistency (audit finding #36): the extras the workflow installs
are defined, the workflow does not run an import-only script, and the integration test does not
depend on the name of the checkout directory.
Run with: pytest tests/test_ci_setup.py -v
"""

import ast
import inspect
import re
from pathlib import Path

import test_variant_caller

ROOT = Path(__file__).resolve().parent.parent
CI = (ROOT / ".github" / "workflows" / "ci.yml").read_text()


def _setup_keywords():
    tree = ast.parse((ROOT / "setup.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "setup":
            return {kw.arg: kw.value for kw in node.keywords}
    raise AssertionError("no setup() call found in setup.py")


class TestExtras:
    def test_tests_extra_is_defined(self):
        extras = ast.literal_eval(_setup_keywords()["extras_require"])
        assert "pytest" in extras["tests"]

    def test_every_extra_the_workflow_installs_exists(self):
        extras = ast.literal_eval(_setup_keywords()["extras_require"])
        requested = re.findall(r"pip install [^\n]*\.\[([\w,]+)\]", CI)
        assert requested, "the workflow installs no extras: update this test"
        for group in requested:
            for extra in group.split(","):
                assert extra in extras, f"ci.yml installs the undefined extra '{extra}'"


class TestWorkflow:
    def test_workflow_does_not_run_an_import_only_script(self):
        for script in re.findall(r"python3? (\S+\.py)", CI):
            path = ROOT / script
            assert path.exists(), f"ci.yml runs {script}, which does not exist"
            body = [l for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
            assert not all(l.startswith(("import ", "from ")) for l in body), f"{script} only imports"

    def test_pytest_runs_the_tests_folder(self):
        assert re.search(r"^\s*pytest [^\n]*\btests\b", CI, re.M)


class TestIntegrationTestPaths:
    def test_data_directory_does_not_depend_on_the_checkout_name(self):
        source = inspect.getsource(test_variant_caller.run_one_genome_test)
        assert '"SigProfilerClusters"' not in source

    def test_data_directory_is_the_tests_folder(self):
        source = inspect.getsource(test_variant_caller.run_one_genome_test)
        assert "os.path.dirname(os.path.abspath(__file__))" in source
        assert (Path(test_variant_caller.__file__).parent / "input_vcf_sim").is_dir()
