from setuptools import setup, find_packages
import os
import re
import shutil


# remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")


def readme():
    this_directory = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(this_directory, "README.md"), encoding="latin-1") as f:
        long_description = f.read()
        return long_description
    # with open('README.rst') as f:
    # 	return(f.read())


def read_version():
    """The version lives only in SigProfilerClusters/version.py (the file that
    sp-release/bump_version.py updates); setup.py never rewrites it."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SigProfilerClusters", "version.py")
    with open(path) as f:
        match = re.search(r"^short_version\s*=\s*['\"]([^'\"]+)['\"]", f.read(), re.M)
    if match is None:
        raise RuntimeError("short_version not found in " + path)
    return match.group(1)


VERSION = read_version()

setup(
    name="SigProfilerClusters",
    version=VERSION,
    description="SigProfilerClusters tool",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="",
    author="Erik Bergstrom",
    author_email="ebergstr@eng.ucsd.edu",
    license="UCSD",
    packages=find_packages(),  # ['SigProfilerMatrixGenerator'],
    python_requires=">=3.9",
    install_requires=[
        "matplotlib>=3.5",
        "sigProfilerPlotting>=1.4.0",
        "SigProfilerMatrixGenerator>=1.3.2",
        "SigProfilerSimulator>=1.2.1",
        "statsmodels>=0.9.0",
        "scipy>=1.1.0",
        "pandas>=2.0.0",
        "numpy>=2.0.0",
        "SigProfilerExtractor>=1.2.0",
        "seaborn>=0.13.0",
    ],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "SigProfilerClusters=SigProfilerClusters.sigprofilerclusters_cli:main_function",
        ],
    },
    extras_require={"tests": ["pytest"]},
    zip_safe=False,
)
