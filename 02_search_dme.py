"""Search small 2D coincidence/domain-matching candidates.

The search follows the logic used in Zur-McGill-style interface matching:
small integer 2D superlattices are generated and compared using area, vector
lengths, and included angle. The final slab atom count is applied during the
search so the 130-atom DFT limit is not accidentally exceeded.

A separate 1D DME diagnostic evaluates the classic m:n domain relation
m*a_film ~= n*a_substrate. That diagnostic can identify scientifically good
DME ratios even when the corresponding 2D periodic DFT cell is too large.
"""
from __future__ import annotations

import json
import math
from collections import defaultdict

import numpy as np

from ase.build import surface
from ase.io import read

from config import (
    A_NA3LACL6_EXPECTED,
    A_ZRO2_EXPECTED,
    DME_1D_MAX_MULTIPLIER,
    DME_ALLOW_NONIMPROVING_BUILD,
    DME_COEFF_LIMIT,
    DME_DIR,
    DME_MAX_ANGLE_DIFF_DEG,
    DME_MAX_AREA_MISMATCH,
    DME_MAX_ATOMS,
    DME_MAX_DET,
    DME_MAX_LENGTH_MISMATCH,
    DME_MAX_VECTOR_LENGTH_A,
    DME_MAX_SUPERCELL_AREA_A2,
    DME_MIN_IMPROVEMENT,
    DME_REQUIRE_IMPROVEMENT_OVER_DIRECT,
    DME_TOP_N,
    NA3LACL6_CIF,
    NA3LACL6_LAYERS,
    ZRO2_CIF,
    ZRO2_LAYERS,
)
from hetero_utils import (
    build_surface_slab,
    convert_hex_surface_basis_to_120,
    load_atoms,
)


def matrix_det(m: tuple[int, int, int, int]) -> int:
    a, b, c, d = m
    return a * d - b * c


def transform_vectors(cell: np.ndarray, m: tuple[int, int, int, int]):
    a11, a12, a21, a22 = m
    v1 = a11 * cell[0] + a12 * cell[1]
    v2 = a21 * cell[0] + a22 * cell[1]
    return v1, v2


def metrics(v1: np.ndarray, v2: np.ndarray) -> tuple[float, float, float, float]:
    l1 = float(np.linalg.norm(v1))
    l2 = float(np.linalg.norm(v2))
    if l1 < 1e-12 or l2 < 1e-12:
        raise ValueError("Degenerate surface vector")
    cos_gamma = np.dot(v1, v2) / (l1 * l2)
    gamma = float(np.degrees(np.arccos(np.clip(cos_gamma, -1.0, 1.0))))
    area = float(np.linalg.norm(np.cross(v1, v2)))
    return l1, l2, gamma, area


def generate_matrices(max_det: int, coeff_limit: int):
    """Generate compact positive-determinant 2D integer matrices."""
    matrices = []
    for a in range(1, coeff_limit + 1):
        for b in range(-coeff_limit, coeff_limit + 1):
            for c in range(-coeff_limit, coeff_limit + 1):
                for d in range(1, coeff_limit + 1):
                    m = (a, b, c, d)
                    det = matrix_det(m)
                    if 0 < det <= max_det:
                        matrices.append((m, det))
    return matrices


def one_d_dme(a_film: float, a_sub: float, max_multiplier: int, z_atoms: int, h_atoms: int):
    results = []
    for m in range(1, max_multiplier + 1):
        for n in range(1, max_multiplier + 1):
            residual = (m * a_film) / (n * a_sub) - 1.0
            results.append(
                {
                    "m": m,
                    "n": n,
                    "film_length_A": m * a_film,
                    "substrate_length_A": n * a_sub,
                    "residual_misfit": residual,
                    "estimated_2d_atoms_if_repeated_both_axes": (
                        z_atoms * m * m + h_atoms * n * n
                    ),
                }
            )
    results.sort(key=lambda x: abs(x["residual_misfit"]))
    return results[:15]


def main() -> None:
    print("=" * 78)
    print("SMALL 2D DME / COINCIDENCE SEARCH")
    print("ZrO2(111) // Na3LaCl6(0001)")
    print("=" * 78)

    zro2 = load_atoms(ZRO2_CIF)
    hex_ = load_atoms(NA3LACL6_CIF)

    z_slab = build_surface_slab(zro2, (1, 1, 1), ZRO2_LAYERS)
    h_slab = build_surface_slab(hex_, (0, 0, 1), NA3LACL6_LAYERS)
    z_slab = convert_hex_surface_basis_to_120(z_slab)

    zcell = z_slab.cell.array
    hcell = h_slab.cell.array
    z_atoms = len(z_slab)
    h_atoms = len(h_slab)

    z111 = A_ZRO2_EXPECTED * math.sqrt(2.0)
    direct_mismatch = (z111 - A_NA3LACL6_EXPECTED) / A_NA3LACL6_EXPECTED

    print(f"Base slab atoms: ZrO2={z_atoms}, Na3LaCl6={h_atoms}")
    print(f"Direct theoretical mismatch: {direct_mismatch * 100:+.4f}%")
    print(f"Final atom limit: {DME_MAX_ATOMS}")

    matrices = generate_matrices(DME_MAX_DET, DME_COEFF_LIMIT)
    grouped = defaultdict(list)

    # Precompute valid surface superlattices using final slab layer counts.
    for matrix, det in matrices:
        atoms = det * z_atoms
        if atoms <= DME_MAX_ATOMS:
            v1, v2 = transform_vectors(zcell, matrix)
            l1, l2, gamma, area = metrics(v1, v2)
            if max(l1, l2) <= DME_MAX_VECTOR_LENGTH_A and area <= DME_MAX_SUPERCELL_AREA_A2:
                grouped["z"].append(
                    {
                        "matrix": matrix,
                        "det": det,
                        "l1": l1,
                        "l2": l2,
                        "gamma": gamma,
                        "area": area,
                        "atoms": atoms,
                    }
                )

        atoms = det * h_atoms
        if atoms <= DME_MAX_ATOMS:
            v1, v2 = transform_vectors(hcell, matrix)
            l1, l2, gamma, area = metrics(v1, v2)
            if max(l1, l2) <= DME_MAX_VECTOR_LENGTH_A and area <= DME_MAX_SUPERCELL_AREA_A2:
                grouped["h"].append(
                    {
                        "matrix": matrix,
                        "det": det,
                        "l1": l1,
                        "l2": l2,
                        "gamma": gamma,
                        "area": area,
                        "atoms": atoms,
                    }
                )

    candidates = []

    for z in grouped["z"]:
        for h in grouped["h"]:
            # det=1 on both sides is only a primitive-basis change and is the
            # direct baseline, not a genuine enlarged DME/coincidence cell.
            if z["det"] == 1 and h["det"] == 1:
                continue

            total_atoms = z["atoms"] + h["atoms"]
            if total_atoms > DME_MAX_ATOMS:
                continue

            angle_diff = abs(z["gamma"] - h["gamma"])
            if angle_diff > DME_MAX_ANGLE_DIFF_DEG:
                continue

            misfit_1 = z["l1"] / h["l1"] - 1.0
            misfit_2 = z["l2"] / h["l2"] - 1.0
            max_abs_misfit = max(abs(misfit_1), abs(misfit_2))
            mean_abs_misfit = 0.5 * (abs(misfit_1) + abs(misfit_2))
            if max_abs_misfit > DME_MAX_LENGTH_MISMATCH:
                continue

            area_mismatch = z["area"] / h["area"] - 1.0
            if abs(area_mismatch) > DME_MAX_AREA_MISMATCH:
                continue

            # Ranking: residual mismatch dominates. Area/angle are secondary;
            # atom count is a tie-breaker, not the primary scientific criterion.
            score = (
                max_abs_misfit
                + 0.30 * abs(area_mismatch)
                + 0.10 * (angle_diff / 180.0)
                + 0.02 * (total_atoms / DME_MAX_ATOMS)
            )

            candidates.append(
                {
                    "z_matrix": list(z["matrix"]),
                    "h_matrix": list(h["matrix"]),
                    "z_det": z["det"],
                    "h_det": h["det"],
                    "z_l1": z["l1"],
                    "z_l2": z["l2"],
                    "h_l1": h["l1"],
                    "h_l2": h["l2"],
                    "misfit_1": misfit_1,
                    "misfit_2": misfit_2,
                    "max_abs_misfit": max_abs_misfit,
                    "mean_abs_misfit": mean_abs_misfit,
                    "area_mismatch": area_mismatch,
                    "angle_difference_deg": angle_diff,
                    "total_atoms": total_atoms,
                    "score": score,
                }
            )

    candidates.sort(key=lambda x: (x["score"], x["total_atoms"]))

    improved = [
        c
        for c in candidates
        if (
            not DME_REQUIRE_IMPROVEMENT_OVER_DIRECT
            or c["max_abs_misfit"] < abs(direct_mismatch) - DME_MIN_IMPROVEMENT
        )
    ]

    best_improved = improved[0] if improved else None
    best_small = candidates[0] if candidates else None

    diagnostic = one_d_dme(
        z111,
        A_NA3LACL6_EXPECTED,
        DME_1D_MAX_MULTIPLIER,
        z_atoms,
        h_atoms,
    )

    status = (
        "VALID_IMPROVED_DME_FOUND"
        if best_improved is not None
        else "NO_IMPROVED_DME_WITHIN_ATOM_LIMIT"
    )

    DME_DIR.mkdir(parents=True, exist_ok=True)

    report = {
        "status": status,
        "orientation": "ZrO2(111)//Na3LaCl6(0001)",
        "direct_theoretical_mismatch": direct_mismatch,
        "direct_theoretical_mismatch_percent": direct_mismatch * 100.0,
        "final_slab_base_atoms": {
            "ZrO2": z_atoms,
            "Na3LaCl6": h_atoms,
        },
        "limits": {
            "max_atoms": DME_MAX_ATOMS,
            "max_det": DME_MAX_DET,
            "max_coeff": DME_COEFF_LIMIT,
            "max_length_mismatch": DME_MAX_LENGTH_MISMATCH,
            "max_area_mismatch": DME_MAX_AREA_MISMATCH,
            "max_angle_difference_deg": DME_MAX_ANGLE_DIFF_DEG,
            "max_vector_length_A": DME_MAX_VECTOR_LENGTH_A,
            "max_supercell_area_A2": DME_MAX_SUPERCELL_AREA_A2,
            "require_improvement_over_direct": DME_REQUIRE_IMPROVEMENT_OVER_DIRECT,
            "min_improvement": DME_MIN_IMPROVEMENT,
            "allow_nonimproving_build": DME_ALLOW_NONIMPROVING_BUILD,
        },
        "best_improved_candidate": best_improved,
        "best_nontrivial_small_candidate": best_small,
        "top_candidates": candidates[:DME_TOP_N],
        "directional_dme_diagnostic": diagnostic,
    }

    (DME_DIR / "dme_report.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    with (DME_DIR / "matching_report.txt").open("w", encoding="utf-8") as f:
        f.write("DME / 2D COINCIDENCE MATCHING REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Status: {status}\n")
        f.write(f"Direct mismatch: {direct_mismatch * 100:+.6f}%\n")
        f.write(f"Final atom limit: {DME_MAX_ATOMS}\n\n")

        f.write("TOP NONTRIVIAL 2D CANDIDATES\n")
        f.write("-" * 80 + "\n")
        if not candidates:
            f.write("No candidates passed the geometric and atom-count filters.\n")
        else:
            for i, c in enumerate(candidates[:DME_TOP_N], 1):
                f.write(
                    f"#{i:02d} Z={tuple(c['z_matrix'])} H={tuple(c['h_matrix'])} "
                    f"atoms={c['total_atoms']} "
                    f"misfit1={c['misfit_1']*100:+.4f}% "
                    f"misfit2={c['misfit_2']*100:+.4f}% "
                    f"area={c['area_mismatch']*100:+.4f}% "
                    f"angle={c['angle_difference_deg']:.3f} deg "
                    f"score={c['score']:.7f}\n"
                )

        f.write("\nDIRECTIONAL DME DIAGNOSTIC\n")
        f.write("-" * 80 + "\n")
        f.write("m*a_ZrO2(111) versus n*a_Na3LaCl6(0001)\n")
        f.write("The estimated 2D atom count assumes the m:n repeat is used in both in-plane directions.\n\n")
        for item in diagnostic:
            f.write(
                f"m={item['m']:2d}, n={item['n']:2d}, "
                f"residual={item['residual_misfit']*100:+.5f}%, "
                f"lengths={item['film_length_A']:.3f}/{item['substrate_length_A']:.3f} A, "
                f"estimated_2D_atoms={item['estimated_2d_atoms_if_repeated_both_axes']}\n"
            )

        if best_improved is None:
            f.write("\nNO DME CANDIDATE IMPROVES ON DIRECT 1x1 WITHIN THE 130-ATOM LIMIT.\n")
            if best_small is not None:
                f.write(
                    "Best nontrivial small candidate: "
                    f"Z={tuple(best_small['z_matrix'])}, "
                    f"H={tuple(best_small['h_matrix'])}, "
                    f"atoms={best_small['total_atoms']}, "
                    f"max_abs_misfit={best_small['max_abs_misfit']*100:.4f}%\n"
                )
            f.write("Keep Direct 1x1 as the primary DFT model.\n")
        else:
            f.write("\nBEST IMPROVED DME CANDIDATE\n")
            f.write(json.dumps(best_improved, indent=2))
            f.write("\n")

    if best_improved is not None:
        (DME_DIR / "best_match.json").write_text(
            json.dumps(best_improved, indent=2),
            encoding="utf-8",
        )
        print("\nBest improved DME candidate:")
        print(json.dumps(best_improved, indent=2))
    else:
        stale = DME_DIR / "best_match.json"
        if stale.exists():
            stale.unlink()
        print("\nNo improved DME candidate found under the 130-atom limit.")

    if best_small is not None:
        (DME_DIR / "best_small_candidate.json").write_text(
            json.dumps(best_small, indent=2),
            encoding="utf-8",
        )

    print(f"Candidates passing filters: {len(candidates)}")
    print(f"Report: {DME_DIR / 'matching_report.txt'}")
    print(f"JSON  : {DME_DIR / 'dme_report.json'}")


if __name__ == "__main__":
    main()
