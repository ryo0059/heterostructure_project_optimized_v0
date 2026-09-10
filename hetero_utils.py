"""Shared geometry utilities for the heterostructure workflow."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from ase import Atoms
from ase.build import surface, make_supercell
from ase.io import read, write


def load_atoms(path: Path | str) -> Atoms:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input structure not found: {path}")
    atoms = read(path)
    if len(atoms) == 0:
        raise ValueError(f"Structure contains zero atoms: {path}")
    return atoms


def build_surface_slab(
    bulk: Atoms,
    miller: Sequence[int],
    layers: int,
    temporary_vacuum: float = 4.0,
) -> Atoms:
    """Build an ASE surface slab and give it a nonzero temporary c vector."""
    slab = surface(bulk, tuple(miller), layers=layers, vacuum=None)

    if len(slab) == 0:
        raise ValueError(f"Surface {tuple(miller)} generated zero atoms")

    zmin = float(np.min(slab.positions[:, 2]))
    zmax = float(np.max(slab.positions[:, 2]))
    thickness = zmax - zmin
    c_len = thickness + 2.0 * temporary_vacuum

    cell = slab.cell.array.copy()
    cell[2] = np.array([0.0, 0.0, c_len])
    slab.set_cell(cell, scale_atoms=False)
    slab.positions[:, 2] -= zmin
    slab.positions[:, 2] += temporary_vacuum
    slab.pbc = (True, True, False)
    slab.wrap(pbc=[True, True, False])
    return slab


def flip_slab_z(slab: Atoms) -> Atoms:
    result = slab.copy()
    z = result.positions[:, 2]
    zmin, zmax = float(z.min()), float(z.max())
    result.positions[:, 2] = zmin + zmax - z
    return result


def surface_cell_metrics(atoms: Atoms) -> dict:
    cell = atoms.cell.array
    a_vec = cell[0].copy()
    b_vec = cell[1].copy()
    a = float(np.linalg.norm(a_vec))
    b = float(np.linalg.norm(b_vec))
    if a < 1e-12 or b < 1e-12:
        raise ValueError("Invalid in-plane cell")
    gamma = float(
        np.degrees(
            np.arccos(
                np.clip(np.dot(a_vec, b_vec) / (a * b), -1.0, 1.0)
            )
        )
    )
    area = float(np.linalg.norm(np.cross(a_vec, b_vec)))
    return {"a": a, "b": b, "gamma": gamma, "area": area}


def convert_hex_surface_basis_to_120(slab: Atoms) -> Atoms:
    """Convert a 60-degree (111) primitive basis to an equivalent 120-degree basis.

    For a cubic (111) surface produced by ASE, the two primitive vectors have
    equal length and gamma=60 deg. Replacing B by B-A gives an equally valid
    primitive surface basis with the same vector length and gamma=120 deg.
    """
    result = slab.copy()
    cell = result.cell.array.copy()
    A = cell[0].copy()
    B = cell[1].copy()
    C = cell[2].copy()

    new_A = A
    new_B = B - A

    scaled = result.get_scaled_positions(wrap=False)
    # old r = x*A + y*B = (x+y)*A + y*(B-A)
    scaled[:, 0] = scaled[:, 0] + scaled[:, 1]

    new_cell = np.array([new_A, new_B, C], dtype=float)
    result.set_cell(new_cell, scale_atoms=False)
    result.set_scaled_positions(scaled)
    result.wrap(pbc=[True, True, False])
    return result


def set_inplane_cell(atoms: Atoms, target_cell: np.ndarray, scale_atoms: bool = True) -> Atoms:
    result = atoms.copy()
    new_cell = result.cell.array.copy()
    new_cell[0] = np.asarray(target_cell[0], dtype=float)
    new_cell[1] = np.asarray(target_cell[1], dtype=float)
    result.set_cell(new_cell, scale_atoms=scale_atoms)
    return result


def rotate_inplane_to_match(source: Atoms, target_cell: np.ndarray) -> Atoms:
    """Rotate source around z so its first in-plane vector points like target[0]."""
    result = source.copy()
    s = result.cell.array[0, :2]
    t = np.asarray(target_cell[0])[:2]
    if np.linalg.norm(s) < 1e-12 or np.linalg.norm(t) < 1e-12:
        return result
    ang_s = np.arctan2(s[1], s[0])
    ang_t = np.arctan2(t[1], t[0])
    angle_deg = float(np.degrees(ang_t - ang_s))
    result.rotate(angle_deg, "z", center=(0.0, 0.0, 0.0), rotate_cell=True)
    return result


def build_2d_supercell(atoms: Atoms, matrix4: Sequence[int]) -> Atoms:
    m = tuple(int(x) for x in matrix4)
    a11, a12, a21, a22 = m
    det = a11 * a22 - a12 * a21
    if det <= 0:
        raise ValueError(f"Transformation determinant must be positive, got {det}: {m}")
    P = np.array(
        [[a11, a12, 0], [a21, a22, 0], [0, 0, 1]],
        dtype=int,
    )
    return make_supercell(atoms, P)


def apply_strain_model(reference: Atoms, film: Atoms, model: str) -> tuple[Atoms, Atoms, dict]:
    """Return reference, film with a common in-plane cell plus strain metadata."""
    ref = reference.copy()
    fil = film.copy()
    ref_cell_before = ref.cell.array.copy()
    fil_cell_before = fil.cell.array.copy()

    if model == "ZRO2_REFERENCE":
        common = ref_cell_before.copy()
        fil = rotate_inplane_to_match(fil, common)
        fil = set_inplane_cell(fil, common, scale_atoms=True)
        reference_strain = np.array([0.0, 0.0])
        film_strain = np.array([
            np.linalg.norm(common[0]) / np.linalg.norm(fil_cell_before[0]) - 1.0,
            np.linalg.norm(common[1]) / np.linalg.norm(fil_cell_before[1]) - 1.0,
        ])
    elif model == "NA3LACL6_REFERENCE":
        common = fil_cell_before.copy()
        ref = rotate_inplane_to_match(ref, common)
        ref = set_inplane_cell(ref, common, scale_atoms=True)
        reference_strain = np.array([
            np.linalg.norm(common[0]) / np.linalg.norm(ref_cell_before[0]) - 1.0,
            np.linalg.norm(common[1]) / np.linalg.norm(ref_cell_before[1]) - 1.0,
        ])
        film_strain = np.array([0.0, 0.0])
    elif model == "SYMMETRIC":
        fil = rotate_inplane_to_match(fil, ref_cell_before)
        common = ref_cell_before.copy()
        common[0] = 0.5 * (ref_cell_before[0] + fil.cell.array[0])
        common[1] = 0.5 * (ref_cell_before[1] + fil.cell.array[1])
        ref = set_inplane_cell(ref, common, scale_atoms=True)
        fil = set_inplane_cell(fil, common, scale_atoms=True)
        reference_strain = np.array([
            np.linalg.norm(common[0]) / np.linalg.norm(ref_cell_before[0]) - 1.0,
            np.linalg.norm(common[1]) / np.linalg.norm(ref_cell_before[1]) - 1.0,
        ])
        film_strain = np.array([
            np.linalg.norm(common[0]) / np.linalg.norm(fil_cell_before[0]) - 1.0,
            np.linalg.norm(common[1]) / np.linalg.norm(fil_cell_before[1]) - 1.0,
        ])
    else:
        raise ValueError(f"Unknown STRAIN_MODEL: {model}")

    metadata = {
        "reference_strain_x": float(reference_strain[0]),
        "reference_strain_y": float(reference_strain[1]),
        "film_strain_x": float(film_strain[0]),
        "film_strain_y": float(film_strain[1]),
        "max_abs_strain": float(max(np.max(np.abs(reference_strain)), np.max(np.abs(film_strain)))),
    }
    return ref, fil, metadata


def shift_bottom_to(atoms: Atoms, z_bottom: float) -> Atoms:
    result = atoms.copy()
    result.positions[:, 2] += z_bottom - float(np.min(result.positions[:, 2]))
    return result


def apply_fractional_xy_translation(atoms: Atoms, fx: float, fy: float) -> Atoms:
    result = atoms.copy()
    shift = fx * result.cell.array[0] + fy * result.cell.array[1]
    result.positions += shift
    return result


def combine_two_slabs(lower: Atoms, upper: Atoms, gap: float, vacuum: float) -> tuple[Atoms, float]:
    lower = shift_bottom_to(lower, 0.0)
    upper = shift_bottom_to(upper, float(np.max(lower.positions[:, 2]) + gap))

    combined = lower + upper
    max_z = float(np.max(combined.positions[:, 2]))
    min_z = float(np.min(combined.positions[:, 2]))
    combined.positions[:, 2] -= min_z
    total_height = max_z - min_z

    cell = combined.cell.array.copy()
    cell[0] = lower.cell.array[0]
    cell[1] = lower.cell.array[1]
    cell[2] = np.array([0.0, 0.0, total_height + vacuum])

    combined.set_cell(cell, scale_atoms=False)
    combined.pbc = (True, True, True)
    combined.wrap()
    return combined, total_height


def minimum_periodic_distance(atoms: Atoms) -> float:
    if len(atoms) < 2:
        return float("inf")
    dist = atoms.get_all_distances(mic=True)
    np.fill_diagonal(dist, np.inf)
    return float(np.min(dist))


def count_species(atoms: Atoms) -> dict:
    result: dict[str, int] = {}
    for symbol in atoms.get_chemical_symbols():
        result[symbol] = result.get(symbol, 0) + 1
    return dict(sorted(result.items()))


def check_expected_lattice(atoms: Atoms, expected_a: float, label: str, tolerance: float = 0.02) -> str:
    a = float(np.linalg.norm(atoms.cell.array[0]))
    rel = abs(a - expected_a) / expected_a
    if rel > tolerance:
        return f"WARNING: {label} a={a:.6f} A differs from expected {expected_a:.6f} A by {rel*100:.3f}%"
    return f"PASS: {label} a={a:.6f} A"


def write_json(path: Path | str, data: dict) -> None:
    path = Path(path)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def read_matrix_json(path: Path | str) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))
