"""Central configuration for the ZrO2/Na3LaCl6 interface workflow."""
from pathlib import Path

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "inputs"
OUTPUT_DIR = ROOT / "outputs"
DIRECT_DIR = OUTPUT_DIR / "direct"
DME_DIR = OUTPUT_DIR / "dme"
VALIDATION_DIR = OUTPUT_DIR / "validation"

for p in (OUTPUT_DIR, DIRECT_DIR, DME_DIR, VALIDATION_DIR):
    p.mkdir(parents=True, exist_ok=True)

ZRO2_CIF = INPUT_DIR / "ZrO2.cif"
NA3LACL6_CIF = INPUT_DIR / "Na3LaCl6.cif"

# -----------------------------------------------------------------------------
# Reference lattice values used for reporting/checks.
# The actual CIF values are used for geometry generation.
# -----------------------------------------------------------------------------
A_ZRO2_EXPECTED = 5.09
A_NA3LACL6_EXPECTED = 7.521

# -----------------------------------------------------------------------------
# Interface geometry
# -----------------------------------------------------------------------------
ZRO2_MILLER = (1, 1, 1)
NA3LACL6_MILLER = (0, 0, 1)

ZRO2_LAYERS = 3
NA3LACL6_LAYERS = 2

INTERFACE_GAP = 2.8  # Angstrom
VACUUM = 15.0        # Angstrom, one periodic vacuum gap

# False -> current ASE (111) surface side is used (for this CIF: O-terminated)
# True  -> flip the slab so the opposite side faces the interface.
ZRO2_FLIP_INTERFACE_SIDE = False

# -----------------------------------------------------------------------------
# In-plane strain convention
# -----------------------------------------------------------------------------
# "ZRO2_REFERENCE": keep ZrO2 supercell fixed and strain Na3LaCl6.
# "NA3LACL6_REFERENCE": keep Na3LaCl6 fixed and strain ZrO2.
# "SYMMETRIC": distribute the difference by taking the Cartesian midpoint cell.
STRAIN_MODEL = "ZRO2_REFERENCE"

MAX_APPLIED_STRAIN = 0.05  # 5%, hard safety limit for generated DFT models

# -----------------------------------------------------------------------------
# Registry sampling
# Fractional translations of the upper Na3LaCl6 slab relative to the common
# in-plane cell. All registries have the same atom count.
# -----------------------------------------------------------------------------
REGISTRY_OFFSETS = (
    (0.0, 0.0),
    (1.0 / 3.0, 1.0 / 3.0),
    (0.5, 0.0),
)

# -----------------------------------------------------------------------------
# DME / 2D coincidence search
# -----------------------------------------------------------------------------
# The search enumerates 2D integer transformation matrices with determinant
# <= this value, then filters using the FINAL slab atom count.
DME_MAX_DET = 6
DME_COEFF_LIMIT = 6
DME_MAX_ATOMS = 130
DME_MAX_LENGTH_MISMATCH = 0.05
DME_MAX_AREA_MISMATCH = 0.10
DME_MAX_ANGLE_DIFF_DEG = 5.0
DME_MAX_VECTOR_LENGTH_A = 30.0
DME_MAX_SUPERCELL_AREA_A2 = 160.0
DME_REQUIRE_IMPROVEMENT_OVER_DIRECT = True
DME_MIN_IMPROVEMENT = 0.001  # absolute fractional improvement = 0.1 percentage point
DME_ALLOW_NONIMPROVING_BUILD = False
DME_TOP_N = 30

# DME literature-style directional m:n diagnostic.
# This is only a diagnostic; it is NOT automatically built into the DFT model.
DME_1D_MAX_MULTIPLIER = 30

# -----------------------------------------------------------------------------
# Geometry validation
# -----------------------------------------------------------------------------
MIN_INTERATOMIC_DISTANCE = 1.50  # Angstrom
MAX_MODEL_ATOMS = 130
