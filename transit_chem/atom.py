from __future__ import annotations
import attr

import numpy as np

from transit_chem.elements import Element


@attr.s
class Atom:
    """Atom class containing coordinates, basis and mass."""

    element: Element = attr.ib()
    x: float = attr.ib()
    y: float = attr.ib()
    z: float = attr.ib()


def translate(atom: Atom, x=0.0, y=0.0, z=0.0) -> Atom:
    return attr.evolve(atom, x=atom.x + x, y=atom.y + y, z=atom.z + z)


def scale(atom: Atom, factor: float) -> Atom:
    return attr.evolve(atom, x=factor * atom.x, y=factor * atom.y, z=factor * atom.z)


def rotate(atom: Atom, R: np.ndarray) -> Atom:
    if not R.shape == (3, 3):
        raise ValueError(f"Rotation matrix R must be 3x3, got {R.shape}")
    x, y, z = R @ np.array([atom.x, atom.y, atom.z])
    return attr.evolve(atom, x=x, y=y, z=z)


def flip_x(atom: Atom) -> Atom:
    return attr.evolve(atom, x=-atom.x)
