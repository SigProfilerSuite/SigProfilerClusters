"""Errors raised by SigProfilerClusters.

They are ordinary exceptions (not SystemExit), so that code calling the package can catch
them; the command line reports them and exits with a code different from 0.
"""


class SigProfilerClustersError(Exception):
    """Base class of the errors that SigProfilerClusters raises on purpose."""


class InvalidParameterError(SigProfilerClustersError, ValueError):
    """A parameter has a value that the analysis does not support."""


class MissingSimulationsError(SigProfilerClustersError):
    """The simulations needed by the analysis are missing or too few."""
