"""Build the selected small DME/coincidence interface.

Run 02_search_dme.py first. By default this script refuses to build a DME cell
unless the selected candidate improves the direct 1x1 residual mismatch and
stays within the 130-atom and 5%-strain limits.
"""
from __future__ import annotations

import json
import sys

import numpy as np
from ase.io import write

from config import (
    DME_DIR,
    DME_MAX_ATOMS,
    DME_ALLOW_NONIMPROVING_BUILD,
    DME_MAX_LENGTH_MISMATCH,
    INTERFACE_GAP,
    MAX_APPLIED_STRAIN,
    MIN_INTERATOMIC_DISTANCE,
    NA3LACL6_CIF,
    NA3LACL6_LAYERS,
    STRAIN_MODEL,
    VACUUM,
    ZRO2_CIF,
    ZRO2_LAYERS,
)
from hetero_utils import (
    apply_fractional_xy_translation,
    apply_strain_model,
    build_2d_supercell,
    build_surface_slab,
    combine_two_slabs,
    convert_hex_surface_basis_to_120,
    count_species,
    load_atoms,
    minimum_periodic_distance,
    surface_cell_metrics,
    write_json,
)


BEST_MATCH = DME_DIR / "best_match.json"
BEST_SMALL = DME_DIR / "best_small_candidate.json"


def build_dme():
    selected_path = BEST_MATCH
    if not selected_path.exists():
        if DME_ALLOW_NONIMPROVING_BUILD and BEST_SMALL.exists():
            selected_path = BEST_SMALL
        else:
            raise RuntimeError(
                "No improved DME candidate exists within the 130-atom constraint. "
                "Run 02_search_dme.py and use Direct 1x1 as the primary model. "
                "Set DME_ALLOW_NONIMPROVING_BUILD=True only if you intentionally want "
                "to build the best non-improving small coincidence cell for comparison."
            )

    best = json.loads(selected_path.read_text(encoding="utf-8"))

    if best["total_atoms"] > DME_MAX_ATOMS:
        raise RuntimeError(
            f"Selected candidate has {best['total_atoms']} atoms > {DME_MAX_ATOMS}."
        )

    zro2_bulk = load_atoms(ZRO2_CIF)
    hex_bulk = load_atoms(NA3LACL6_CIF)

    zro2 = build_surface_slab(zro2_bulk, (1, 1, 1), ZRO2_LAYERS)
    hex = build_surface_slab(hex_bulk, (0, 0, 1), NA3LACL6_LAYERS)
    zro2 = convert_hex_surface_basis_to_120(zro2)

    zro2 = build_2d_supercell(zro2, best["z_matrix"])
    hex = build_2d_supercell(hex, best["h_matrix"])

    # Apply the selected strain convention. For the default ZRO2_REFERENCE,
    # Na3LaCl6 is rotated into the ZrO2 in-plane orientation and then strained.
    zro2, hex, strain = apply_strain_model(
        reference=zro2,
        film=hex,
        model=STRAIN_MODEL,
    )

    if strain["max_abs_strain"] > MAX_APPLIED_STRAIN + 1e-12:
        raise RuntimeError(
            f"Applied strain {strain['max_abs_strain']*100:.3f}% exceeds "
            f"the configured limit of {MAX_APPLIED_STRAIN*100:.2f}%."
        )

    hetero, physical_height = combine_two_slabs(
        lower=zro2,
        upper=hex,
        gap=INTERFACE_GAP,
        vacuum=VACUUM,
    )

    min_dist = minimum_periodic_distance(hetero)
    if min_dist < MIN_INTERATOMIC_DISTANCE:
        raise RuntimeError(
            f"Minimum periodic distance {min_dist:.3f} A is below "
            f"{MIN_INTERATOMIC_DISTANCE:.3f} A."
        )

    metrics = surface_cell_metrics(hetero)

    metadata = {
        "model": "DME/coincidence interface",
        "orientation": "ZrO2(111)//Na3LaCl6(0001)",
        "z_matrix": list(best["z_matrix"]),
        "h_matrix": list(best["h_matrix"]),
        "z_det": best["z_det"],
        "h_det": best["h_det"],
        "input_search_length_mismatch": best["max_abs_misfit"],
        "input_search_area_mismatch": best["area_mismatch"],
        "input_search_angle_difference_deg": best["angle_difference_deg"],
        "atom_count": int(len(hetero)),
        "species": count_species(hetero),
        "cell_a_A": metrics["a"],
        "cell_b_A": metrics["b"],
        "cell_gamma_deg": metrics["gamma"],
        "cell_area_A2": metrics["area"],
        "physical_stack_height_A": physical_height,
        "vacuum_gap_A": VACUUM,
        "interface_gap_A": INTERFACE_GAP,
        "minimum_periodic_distance_A": min_dist,
        "strain_model": STRAIN_MODEL,
        **strain,
        "pbc": [True, True, True],
    }

    output = DME_DIR / "Hetero_DME.cif"
    output_xyz = DME_DIR / "Hetero_DME.xyz"
    output_json = DME_DIR / "Hetero_DME.json"

    write(output, hetero)
    write(output_xyz, hetero)
    write_json(output_json, metadata)

    return output, metadata


def main():
    print("=" * 78)
    print("BUILD SELECTED DME / COINCIDENCE INTERFACE")
    print("=" * 78)

    output, meta = build_dme()

    print(f"Atoms            : {meta['atom_count']}")
    print(f"Cell a,b         : {meta['cell_a_A']:.6f}, {meta['cell_b_A']:.6f} A")
    print(f"Cell gamma       : {meta['cell_gamma_deg']:.3f} deg")
    print(f"Applied max strain: {meta['max_abs_strain']*100:.3f}%")
    print(f"Minimum distance : {meta['minimum_periodic_distance_A']:.3f} A")
    print(f"Output           : {output}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
