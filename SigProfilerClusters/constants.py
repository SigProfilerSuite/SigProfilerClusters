"""Scientific thresholds of SigProfilerClusters, in one place (audit finding 22).

These values were literals scattered through the code. Naming them does not change any result.
They are method constants, not user options: modules read them as `constants.NAME` at the point
of use, so a test can change one and observe the effect.
"""

# Sample-specific IMD cutoff (hotspot.first_run and hotspot.refineIMD)
IMD_CUTOFF_CAP = 10000
"""The sample IMD cutoff never exceeds this distance (bp)."""

CLUSTERED_FRACTION = 0.90
"""Minimum fraction of observed mutations among observed + simulated ones up to the cutoff."""

CLUSTERED_FRACTION_CORRECTED = 0.95
"""The same fraction for the cutoffs of dense regions (regional density correction)."""

Q_VALUE_THRESHOLD = 0.01
"""Significance threshold (FDR q-value) of the enrichment over the simulations."""

# Grouping and subclassification of clustered mutations (classifyFunctions)
EVENT_GROUPING_DISTANCE = 10001
"""Consecutive clustered mutations closer than this distance (bp) belong to the same event."""

VAF_MAX_DIFFERENCE = 0.1
"""Largest VAF difference between adjacent mutations of an event."""

CCF_MAX_DIFFERENCE = 0.25
"""Largest CCF difference between adjacent mutations of an event (includedCCFs=True)."""

PROCESSIVE_RUN_FRACTION = 0.75
"""A run of one base class longer than this fraction of the event makes it processive."""

PROCESSIVE_PURITY = 0.9
"""Fraction of one base class (C/G or T/A) that makes an event processive."""

# Regional density correction (hotspot.densityCorrection and densityCorrectionOriginal)
DENSITY_FOLD_CHANGE = 1.25
"""Minimum fold change of the observed over the simulated density to report a dense region."""

DENSITY_MIN_MUTATIONS = 50
"""Minimum number of observed mutations in a dense window."""

DENSITY_PEAK_THRESHOLD = 9.0
"""Peak height used by densityCorrection to find dense windows."""

# Event probability (SigProfilerClusters.eventProbability)
EVENT_PROBABILITY_FLOOR = 1e-05
"""Smallest score written to the EVENT_PROBABILITY column."""

# Values that analysis() accepts (audit finding 29)
SUPPORTED_ANALYSES = ("all", "hotspot", "subClassify")
# Names of the mutation contexts that analysis() accepts as input (audit findings 118, 119).
# The standard names are SBS96 and ID83; the others are accepted as aliases. DBS is not an
# input: it is one of the outputs of SigProfilerClusters.
SBS_CONTEXT = "SBS96"
INDEL_CONTEXT = "ID83"
CONTEXT_ALIASES = {
    SBS_CONTEXT: ("SBS96", "96", "SBS"),
    INDEL_CONTEXT: ("ID83", "ID", "83", "INDEL"),
}


def canonical_context(value):
    """The standard name (SBS96 or ID83) of a context given with any of its aliases, or None
    when the value is not a mutation context that analysis() supports."""
    if not isinstance(value, str):
        return None
    for standard, aliases in CONTEXT_ALIASES.items():
        if value in aliases:
            return standard
    return None


def is_indel_context(value):
    """True for ID83 and its aliases (ID, 83, INDEL)"""
    return canonical_context(value) == INDEL_CONTEXT


def is_sbs_context(value):
    """True for SBS96 and its aliases (96, SBS)"""
    return canonical_context(value) == SBS_CONTEXT


VARIANT_CALLERS = ("standard", "caveman", "mutect2", "custom")
INPUT_FILE_EXTENSIONS = ("txt", "vcf", "maf", "tsv")
