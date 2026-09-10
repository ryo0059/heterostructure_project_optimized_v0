"""
04_validate_structures.py

Validasi struktur heterostruktur:
ZrO2(111) // Na3LaCl6(0001)

Yang dicek:
1. Jumlah atom
2. Cell lengths dan angles
3. PBC
4. Minimum global interatomic distance
5. Geometric interlayer gap
6. Nearest cross-interface atom distance
7. Estimasi total vacuum

Catatan:
Rentang 3–5 Å digunakan sebagai screening geometrik awal.
Ini BUKAN bukti bahwa interface merupakan pure van der Waals.
Jarak final harus diperiksa kembali setelah structural relaxation.
"""

from pathlib import Path

import numpy as np
from ase.io import read

from config import (
    OUTPUT_DIR,
    VALIDATION_DIR,
    MIN_INTERATOMIC_DISTANCE,
    DME_MAX_ATOMS,
)


# ============================================================
# ELEMENT GROUP
# ============================================================

ZRO2_ELEMENTS = {"Zr", "O"}
NA3LACL6_ELEMENTS = {"Na", "La", "Cl"}


# Screening range untuk initial interface gap
INTERFACE_GAP_MIN = 3.0
INTERFACE_GAP_MAX = 5.0


# ============================================================
# BASIC CELL INFORMATION
# ============================================================

def cell_summary(atoms):

    lengths = atoms.cell.lengths()
    angles = atoms.cell.angles()

    return lengths, angles


# ============================================================
# GLOBAL MINIMUM DISTANCE
# ============================================================

def min_global_distance(atoms):

    n = len(atoms)

    if n < 2:
        return np.inf

    dmin = np.inf

    for i in range(n):

        indices = np.arange(n)

        distances = atoms.get_distances(
            i,
            indices,
            mic=True
        )

        # Jangan hitung atom terhadap dirinya sendiri
        distances[i] = np.inf

        local_min = np.min(distances)

        if local_min < dmin:
            dmin = local_min

    return float(dmin)


# ============================================================
# IDENTIFY MATERIALS
# ============================================================

def material_indices(atoms):

    symbols = np.array(
        atoms.get_chemical_symbols()
    )

    zro2_indices = np.array(
        [
            i
            for i, symbol in enumerate(symbols)
            if symbol in ZRO2_ELEMENTS
        ],
        dtype=int
    )

    na3lacl6_indices = np.array(
        [
            i
            for i, symbol in enumerate(symbols)
            if symbol in NA3LACL6_ELEMENTS
        ],
        dtype=int
    )

    return zro2_indices, na3lacl6_indices


# ============================================================
# INTERFACE ANALYSIS
# ============================================================

def interface_analysis(atoms):

    zro2_indices, na3lacl6_indices = material_indices(atoms)

    if (
        len(zro2_indices) == 0
        or len(na3lacl6_indices) == 0
    ):
        return None

    positions = atoms.get_positions()

    z_zro2 = positions[
        zro2_indices,
        2
    ]

    z_na3lacl6 = positions[
        na3lacl6_indices,
        2
    ]

    # --------------------------------------------------------
    # GEOMETRIC INTERLAYER GAP
    #
    # Diasumsikan:
    #
    # ZrO2 berada di bawah
    # Na3LaCl6 berada di atas
    #
    # gap = bottom Na3LaCl6 - top ZrO2
    # --------------------------------------------------------

    top_zro2 = np.max(z_zro2)

    bottom_na3lacl6 = np.min(z_na3lacl6)

    geometric_gap = (
        bottom_na3lacl6
        - top_zro2
    )

    # --------------------------------------------------------
    # ATOM TERDEKAT DI INTERFACE
    #
    # Kita hanya mengambil:
    #
    # Zr/O bagian atas
    # Na/La/Cl bagian bawah
    #
    # sehingga tidak tertukar dengan atom di bagian dalam slab.
    # --------------------------------------------------------

    zro2_surface_cutoff = np.median(z_zro2)

    na3lacl6_surface_cutoff = np.median(
        z_na3lacl6
    )

    zro2_surface = zro2_indices[
        z_zro2 >= zro2_surface_cutoff
    ]

    na3lacl6_surface = na3lacl6_indices[
        z_na3lacl6 <= na3lacl6_surface_cutoff
    ]

    nearest_distance = None
    nearest_pair = None

    for i in zro2_surface:

        distances = atoms.get_distances(
            i,
            na3lacl6_surface,
            mic=True
        )

        if len(distances) == 0:
            continue

        local_index = int(
            np.argmin(distances)
        )

        distance = float(
            distances[local_index]
        )

        if (
            nearest_distance is None
            or distance < nearest_distance
        ):

            nearest_distance = distance

            nearest_pair = (
                int(i),
                int(
                    na3lacl6_surface[
                        local_index
                    ]
                )
            )

    return {
        "geometric_gap": float(
            geometric_gap
        ),
        "cross_interface_distance": nearest_distance,
        "pair": nearest_pair,
    }


# ============================================================
# VACUUM ANALYSIS
# ============================================================

def vacuum_analysis(atoms):

    positions = atoms.get_positions()

    c_length = float(
        atoms.cell.lengths()[2]
    )

    z_min = np.min(
        positions[:, 2]
    )

    z_max = np.max(
        positions[:, 2]
    )

    z_extent = (
        z_max
        - z_min
    )

    total_vacuum = (
        c_length
        - z_extent
    )

    return (
        float(z_extent),
        float(total_vacuum)
    )


# ============================================================
# GAP CLASSIFICATION
# ============================================================

def classify_gap(gap):

    if gap is None:
        return "UNKNOWN"

    if gap < INTERFACE_GAP_MIN:
        return "BELOW_3_A"

    if gap <= INTERFACE_GAP_MAX:
        return "WITHIN_3_TO_5_A"

    return "ABOVE_5_A"


# ============================================================
# VALIDATE ONE FILE
# ============================================================

def validate_file(path):

    atoms = read(path)

    lengths, angles = cell_summary(
        atoms
    )

    global_dmin = min_global_distance(
        atoms
    )

    interface = interface_analysis(
        atoms
    )

    z_extent, total_vacuum = vacuum_analysis(
        atoms
    )

    # --------------------------------------------------------
    # HEADER
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("HETEROSTRUCTURE VALIDATION")
    print("=" * 72)

    print(f"File                   : {path}")

    # --------------------------------------------------------
    # BASIC STRUCTURE
    # --------------------------------------------------------

    print()
    print("STRUCTURE")
    print("-" * 72)

    print(
        f"Total atoms            : "
        f"{len(atoms)}"
    )

    print(
        "Cell lengths            : "
        f"a={lengths[0]:.4f} Å, "
        f"b={lengths[1]:.4f} Å, "
        f"c={lengths[2]:.4f} Å"
    )

    print(
        "Cell angles             : "
        f"α={angles[0]:.3f}°, "
        f"β={angles[1]:.3f}°, "
        f"γ={angles[2]:.3f}°"
    )

    print(
        f"PBC                    : "
        f"{atoms.pbc}"
    )

    # --------------------------------------------------------
    # ATOM COUNT
    # --------------------------------------------------------

    print()
    print("ATOM COUNT")
    print("-" * 72)

    if len(atoms) <= DME_MAX_ATOMS:

        print(
            f"[OK] {len(atoms)} atoms "
            f"<= limit {DME_MAX_ATOMS}"
        )

    else:

        print(
            f"[WARNING] {len(atoms)} atoms "
            f"> limit {DME_MAX_ATOMS}"
        )

    # --------------------------------------------------------
    # GLOBAL DISTANCE
    # --------------------------------------------------------

    print()
    print("GLOBAL INTERATOMIC DISTANCE")
    print("-" * 72)

    print(
        f"Minimum global distance : "
        f"{global_dmin:.4f} Å"
    )

    if global_dmin < MIN_INTERATOMIC_DISTANCE:

        print(
            f"[FAIL] Distance is below "
            f"{MIN_INTERATOMIC_DISTANCE:.2f} Å"
        )

    else:

        print(
            f"[OK] Distance >= "
            f"{MIN_INTERATOMIC_DISTANCE:.2f} Å"
        )

    # --------------------------------------------------------
    # INTERFACE
    # --------------------------------------------------------

    print()
    print("INTERFACE ANALYSIS")
    print("-" * 72)

    if interface is None:

        print(
            "[FAIL] Could not identify "
            "Zr/O and Na/La/Cl groups."
        )

    else:

        gap = interface[
            "geometric_gap"
        ]

        cross_distance = interface[
            "cross_interface_distance"
        ]

        pair = interface[
            "pair"
        ]

        print(
            f"Top ZrO2 → bottom Na3LaCl6 "
            f"gap                    : "
            f"{gap:.4f} Å"
        )

        print(
            f"Gap classification      : "
            f"{classify_gap(gap)}"
        )

        if cross_distance is not None:

            symbols = (
                atoms.get_chemical_symbols()
            )

            i, j = pair

            print(
                f"Nearest cross-interface : "
                f"{cross_distance:.4f} Å"
            )

            print(
                f"Nearest atom pair       : "
                f"{symbols[i]}-{symbols[j]}"
            )

        else:

            print(
                "Nearest cross-interface : "
                "N/A"
            )

        # ----------------------------------------------------
        # GAP WARNING
        # ----------------------------------------------------

        if gap < INTERFACE_GAP_MIN:

            print(
                "[WARNING] Interface gap < 3 Å."
            )

            print(
                "Inspect the interface in VESTA "
                "for possible atomic overlap."
            )

        elif gap <= INTERFACE_GAP_MAX:

            print(
                "[OK] Initial geometric gap "
                "is within 3–5 Å."
            )

        else:

            print(
                "[INFO] Interface gap > 5 Å."
            )

            print(
                "This is not automatically wrong; "
                "inspect the interface and later "
                "compare after relaxation."
            )

    # --------------------------------------------------------
    # VACUUM
    # --------------------------------------------------------

    print()
    print("VACUUM")
    print("-" * 72)

    print(
        f"Slab z-extent           : "
        f"{z_extent:.4f} Å"
    )

    print(
        f"Estimated total vacuum  : "
        f"{total_vacuum:.4f} Å"
    )

    # --------------------------------------------------------
    # SCIENTIFIC NOTE
    # --------------------------------------------------------

    print()
    print("SCIENTIFIC NOTE")
    print("-" * 72)

    print(
        "The 3–5 Å range is used only as an "
        "initial geometric screening criterion."
    )

    print(
        "It does NOT by itself prove a pure "
        "van der Waals interface."
    )

    print(
        "The final interface distance should "
        "be checked again after structural "
        "relaxation."
    )

    print("=" * 72)

    return {
        "file": str(path),
        "atoms": len(atoms),
        "minimum_global_distance": global_dmin,
        "interface_gap": (
            None
            if interface is None
            else interface["geometric_gap"]
        ),
        "nearest_cross_interface_distance": (
            None
            if interface is None
            else interface[
                "cross_interface_distance"
            ]
        ),
        "vacuum_total": total_vacuum,
        "cell_lengths": lengths.tolist(),
        "cell_angles": angles.tolist(),
    }


# ============================================================
# MAIN
# ============================================================

def main():

    output_root = Path(
        OUTPUT_DIR
    )

    validation_root = Path(
        VALIDATION_DIR
    )

    validation_root.mkdir(
        parents=True,
        exist_ok=True
    )

    cif_files = sorted(
        output_root.rglob("*.cif")
    )

    if not cif_files:

        print(
            "Tidak ditemukan file CIF "
            "di dalam outputs/."
        )

        return

    results = []

    for path in cif_files:

        try:

            result = validate_file(
                path
            )

            results.append(
                result
            )

        except Exception as exc:

            print()
            print("=" * 72)
            print(
                f"ERROR VALIDATING: {path}"
            )
            print("=" * 72)

            print(exc)

    # ========================================================
    # SAVE REPORT
    # ========================================================

    report_path = (
        validation_root
        / "validation_report.txt"
    )

    with report_path.open(
        "w",
        encoding="utf-8"
    ) as f:

        f.write(
            "HETEROSTRUCTURE VALIDATION REPORT\n"
        )

        f.write(
            "=" * 72
            + "\n\n"
        )

        for result in results:

            f.write(
                f"FILE: {result['file']}\n"
            )

            f.write(
                f"Atoms: "
                f"{result['atoms']}\n"
            )

            f.write(
                "Minimum global distance: "
                f"{result['minimum_global_distance']:.4f} Å\n"
            )

            if result[
                "interface_gap"
            ] is not None:

                f.write(
                    "Geometric interface gap: "
                    f"{result['interface_gap']:.4f} Å\n"
                )

            else:

                f.write(
                    "Geometric interface gap: N/A\n"
                )

            if result[
                "nearest_cross_interface_distance"
            ] is not None:

                f.write(
                    "Nearest cross-interface distance: "
                    f"{result['nearest_cross_interface_distance']:.4f} Å\n"
                )

            else:

                f.write(
                    "Nearest cross-interface distance: N/A\n"
                )

            f.write(
                "Estimated total vacuum: "
                f"{result['vacuum_total']:.4f} Å\n"
            )

            f.write(
                "Cell lengths: "
                + ", ".join(
                    f"{x:.4f}"
                    for x in result[
                        "cell_lengths"
                    ]
                )
                + " Å\n"
            )

            f.write(
                "Cell angles: "
                + ", ".join(
                    f"{x:.3f}"
                    for x in result[
                        "cell_angles"
                    ]
                )
                + " deg\n"
            )

            f.write(
                "\n"
            )

    print()
    print("=" * 72)
    print(
        "Validation report saved:"
    )
    print(
        report_path
    )
    print("=" * 72)


if __name__ == "__main__":
    main()