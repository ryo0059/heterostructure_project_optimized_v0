# Local Test Result

The optimized scripts were syntax-checked and the geometry workflow was exercised using the supplied CIFs with ASE 3.29.0 source and the working NumPy environment.

## Direct model

- Orientation: ZrO2(111) // Na3LaCl6(0001)
- ZrO2 layers: 3
- Na3LaCl6 layers: 2
- Total atoms: 64
- Common in-plane cell: 7.198553 x 7.198553 Å, gamma = 120°
- Initial interface gap: 2.8 Å
- Periodic vacuum gap: 15 Å
- Maximum coherent strain under `ZRO2_REFERENCE`: 4.287%
- Minimum periodic interatomic distance in the generated direct model: 1.855 Å

## DME search

The search is performed with the final 3-layer/2-layer atom counts and a hard limit of 130 atoms.

For the present lattice parameters, no nontrivial DME/coincidence cell within the configured limits improves the direct 1x1 residual mismatch.

The directional DME diagnostic identifies very large near-zero-residual integer ratios (for example 23:22 and 24:23), but a symmetric 2D repetition of those domains would require tens of thousands of atoms and is therefore not suitable for the present DFT resource constraint.

The intended primary DFT structure is therefore the 64-atom direct interface.
