"""Version and release metadata consistency (audit finding #35).
The version has one source (SigProfilerClusters/version.py); setup.py reads it and must not
rewrite it, and the CHANGELOG must list releases in chronological order.
Run with: pytest tests/test_version_metadata.py -v
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

import SigProfilerClusters

ROOT = Path(__file__).resolve().parent.parent
HEADING = re.compile(r"^## \[(\d+\.\d+\.\d+)\] - (\d{4}-\d{2}-\d{2})", re.M)


def _version_tuple(text):
    return tuple(int(x) for x in text.split("."))


def _setup_version(directory):
    result = subprocess.run(
        [sys.executable, "setup.py", "--version"],
        cwd=str(directory),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=True,
    )
    return result.stdout.strip().splitlines()[-1]


def _copy_packaging_files(tmp_path):
    (tmp_path / "SigProfilerClusters").mkdir()
    shutil.copy(str(ROOT / "setup.py"), str(tmp_path / "setup.py"))
    shutil.copy(str(ROOT / "README.md"), str(tmp_path / "README.md"))
    shutil.copy(
        str(ROOT / "SigProfilerClusters" / "version.py"),
        str(tmp_path / "SigProfilerClusters" / "version.py"),
    )
    return tmp_path


class TestVersionSource:
    def test_setup_and_package_report_the_same_version(self, tmp_path):
        assert _setup_version(_copy_packaging_files(tmp_path)) == SigProfilerClusters.__version__

    def test_setup_does_not_rewrite_version_py(self, tmp_path):
        copy = _copy_packaging_files(tmp_path)
        target = copy / "SigProfilerClusters" / "version.py"
        target.write_text("short_version = '7.8.9'\nversion = '7.8.9'\n")
        before = target.read_text()
        assert _setup_version(copy) == "7.8.9"
        assert target.read_text() == before

    def test_version_py_is_not_marked_as_generated(self):
        text = (ROOT / "SigProfilerClusters" / "version.py").read_text()
        assert "GENERATED" not in text


class TestChangelogChronology:
    def _entries(self):
        return HEADING.findall((ROOT / "CHANGELOG.md").read_text())

    def test_the_changelog_has_releases(self):
        assert len(self._entries()) >= 3

    def test_releases_are_listed_newest_first(self):
        dates = [date for _, date in self._entries()]
        assert dates == sorted(dates, reverse=True)

    def test_versions_decrease_down_the_file(self):
        versions = [_version_tuple(v) for v, _ in self._entries()]
        assert versions == sorted(versions, reverse=True)

    def test_latest_entry_is_not_newer_than_the_package(self):
        latest = self._entries()[0][0]
        assert _version_tuple(latest) <= _version_tuple(SigProfilerClusters.__version__)
