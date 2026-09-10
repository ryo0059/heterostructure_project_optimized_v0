# ZrO2(111) / Na3LaCl6(0001) Heterostructure Workflow

This project builds a small periodic heterointerface for computational study:

**ZrO2 (111) // Na3LaCl6 (0001)**

The workflow is intentionally designed around a **maximum final model size of 130 atoms**, which is the current practical server limit for the intended calculations.

## Scientific strategy

The workflow separates two ideas that were previously mixed together:

1. **Direct coherent 1x1 interface** — the primary DFT model.
2. **DME / 2D coincidence search** — a screening step that looks for small integer superlattices with lower residual mismatch.

The DME literature describes domain matching as matching integral multiples of lattice planes/domains across an interface. For small mismatch, conventional coherent lattice matching is treated as a limiting/special case of the broader epitaxy framework. The 2D search here follows the same general logic as Zur-McGill-style interface matching: enumerate integer superlattices and compare in-plane area, vector lengths, and angle. It is intentionally constrained by the final atom count.

## Important result for the present CIFs

For the supplied structures:

- ZrO2 bulk lattice parameter is approximately 5.09 Å.
- Na3LaCl6 has a = b = 7.521 Å and gamma = 120°.
- The ZrO2(111) in-plane primitive length is:

  `a_111 = sqrt(2) * a_ZrO2`

  which is about 7.199 Å.

- Direct 1x1 mismatch is about **-4.29%**.

With 3 ZrO2 surface layers and 2 Na3LaCl6 surface layers, the generated direct model contains **64 atoms** before any DFT relaxation.

The current atom limit is 130 atoms.

The DME search shows that a substantially lower residual mismatch requires a much larger domain ratio. For example, the first very small-residual directional ratios are around 23:22 and 24:23, but a symmetric 2D repetition of those domains would require tens of thousands of atoms for the current slab thickness. Therefore, those large DME cells are not suitable for the present DFT server constraint.

The correct scientific outcome is therefore **not to force a large DME cell into the 130-atom model**. The 64-atom direct interface is the primary model, while DME is retained as a documented screening/diagnostic result.

## Files

```text
heterostructure_project/
│
├── inputs/
│   ├── ZrO2.cif
│   └── Na3LaCl6.cif
│
├── outputs/
│   ├── direct/
│   ├── dme/
│   └── validation/
│
├── config.py
├── hetero_utils.py
├── 01_build_direct.py
├── 02_search_dme.py
├── 03_build_dme_hetero.py
├── 04_validate_structures.py
├── requirements.txt
└── README.md
```

## Installation

Use a fresh virtual environment when possible.

### Windows

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

### Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Configuration

All principal settings are in `config.py`.

Important values:

```python
ZRO2_LAYERS = 3
NA3LACL6_LAYERS = 2

INTERFACE_GAP = 2.8
VACUUM = 15.0

DME_MAX_ATOMS = 130
DME_MAX_DET = 6
DME_MAX_LENGTH_MISMATCH = 0.05
DME_MAX_AREA_MISMATCH = 0.10
DME_MAX_ANGLE_DIFF_DEG = 5.0
```

### Strain convention

The default is:

```python
STRAIN_MODEL = "ZRO2_REFERENCE"
```

This means:

- ZrO2 is treated as the reference in-plane lattice.
- Na3LaCl6 is strained to the common in-plane cell.

This is a **modeling assumption**, not a universal standard. If a different physical assumption is desired, the configuration supports alternative strain conventions.

## Step 1 — Build direct interface

Run:

```bash
python 01_build_direct.py
```

The script:

- reads both CIF files;
- creates ZrO2(111) and Na3LaCl6(0001) slabs;
- converts the ZrO2(111) 60° primitive basis to an equivalent 120° basis using `B -> B - A`;
- applies the selected in-plane strain convention;
- constructs the interface with a 2.8 Å initial gap;
- adds a 15 Å periodic vacuum gap;
- checks the minimum periodic interatomic distance;
- generates several representative lateral registries.

The three registry offsets are:

```text
(0, 0)
(1/3, 1/3)
(1/2, 0)
```

All have the same atom count. They are useful because interface energy can depend on the lateral registry after relaxation.

### Direct outputs

```text
outputs/direct/
├── Hetero_Direct_ZrO2_111_Na3LaCl6_0001.cif
├── Hetero_Direct_ZrO2_111_Na3LaCl6_0001.xyz
├── Hetero_Direct_registry_01.cif
├── Hetero_Direct_registry_02.cif
├── Hetero_Direct_registry_03.cif
└── direct_summary.json
```

Open the CIF files in VESTA before doing any DFT calculation.

## Step 2 — Search for a small DME/coincidence cell

Run:

```bash
python 02_search_dme.py
```

The script does **not** automatically build a heterostructure.

It searches small 2D integer transformations while enforcing the **final slab atom count**. This fixes a major problem in the earlier workflow, where the search used one-layer atom counts and could therefore select a candidate that exceeded the true 3-layer/2-layer DFT size.

The search also rejects unnecessarily long/skew supercells using maximum vector length and area limits.

### DME outputs

```text
outputs/dme/
├── matching_report.txt
├── dme_report.json
└── best_small_candidate.json
```

If an improved DME cell is found within the limits, `best_match.json` is also created.

If `best_match.json` does not exist, that is an expected scientific result when no improved DME cell fits the atom limit.

## Step 3 — Build DME interface

Run:

```bash
python 03_build_dme_hetero.py
```

By default, the script only builds an **improved** DME candidate.

If no improved candidate exists within 130 atoms, the script intentionally stops rather than silently building a non-improving or excessively strained cell.

To intentionally build the best small non-improving coincidence candidate for comparison, change:

```python
DME_ALLOW_NONIMPROVING_BUILD = False
```

to:

```python
DME_ALLOW_NONIMPROVING_BUILD = True
```

That option should be treated as a comparison/diagnostic model, not automatically as the preferred physical model.

## Step 4 — Validate generated structures

Run:

```bash
python 04_validate_structures.py
```

The validator checks:

- atom count <= 130;
- minimum periodic interatomic distance;
- periodic boundary conditions;
- fractional x/y coordinates inside the cell;
- z-axis cell geometry;
- periodic vacuum gap.

Outputs:

```text
outputs/validation/
├── validation_report.txt
└── validation_report.json
```

The validator exits with a non-zero status when an existing generated structure fails a required check.

## VESTA inspection

For the direct interface, first open:

```text
outputs/direct/Hetero_Direct_ZrO2_111_Na3LaCl6_0001.cif
```

Check:

- ZrO2 is below Na3LaCl6;
- the interface is parallel to x-y;
- the initial interface separation is approximately 2.8 Å;
- the periodic vacuum region is approximately 15 Å;
- there are no visibly overlapping atoms;
- the periodic cell is reasonable for DFT.

Also inspect `registry_02` and `registry_03` when studying interface registry effects.

## Why the 1x1 direct model is retained

For the present lattice parameters, direct coherent matching already gives a mismatch of about 4.29%.

That is small enough to make a direct coherent model scientifically reasonable as the primary reference model. The DME literature explicitly places conventional lattice matching within the same broader epitaxy framework for small misfit systems.

A near-zero-residual DME ratio exists only at much larger integer domains. Those domains rapidly increase the periodic interface area and therefore the atom count. Under the current **130-atom limit**, forcing such a DME model would defeat the computational purpose.

## Recommended research workflow

```text
CIF
 ↓
01_build_direct.py
 ↓
VESTA inspection
 ↓
04_validate_structures.py
 ↓
Direct 64-atom reference model
 ↓
02_search_dme.py
 ↓
Check whether an improved DME model fits <=130 atoms
 ↓
If yes → 03_build_dme_hetero.py
If no  → keep Direct 1x1 as primary model
 ↓
Final VESTA inspection
 ↓
DFT relaxation
 ↓
SCF / electronic structure / DOS / PDOS / charge analysis
```

## Important methodological note

The initial interface gap and lattice strain are **starting-model parameters**. They are not final equilibrium values. The final interfacial geometry should be obtained from structural relaxation in the chosen DFT method.

Likewise, the choice of which material receives the coherent in-plane strain is a modeling assumption and should be reported explicitly in a thesis or paper.

## Literature/method references

- Narayan and Larson, *Domain epitaxy: A unified paradigm for thin film growth*, Journal of Applied Physics 93, 278–285 (2003), DOI: 10.1063/1.1528301.
- Zur and McGill, lattice matching methodology for heterostructural interfaces, Journal of Applied Physics 55 (1984), 378, DOI: 10.1063/1.333084.
- ASE surface construction follows the documented Miller-index slab workflow.
- Pymatgen's modern interface tools implement Zur-McGill-style 2D lattice matching and can be used as an independent cross-check of the matching search.
