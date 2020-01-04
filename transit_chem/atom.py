from __future__ import annotations

from math import cos, isclose, sin, sqrt

import attr
import numpy as np

from transit_chem.config import conf
from transit_chem.elements import Element, element_from_string
from transit_chem.utils import angle


@attr.s(frozen=True, cmp=False)
class Atom:
    """Atom class containing coordinates, basis and mass."""

    element: Element = attr.ib(converter=element_from_string)
    x: float = attr.ib(converter=float)
    y: float = attr.ib(converter=float)
    z: float = attr.ib(converter=float)

    @property
    def mass(self):
        return self.element.z

    @property
    def coords(self):
        return np.array([self.x, self.y, self.z])

    def with_coords(self, x: float, y: float, z: float) -> "Atom":
        return attr.evolve(self, x=x, y=y, z=z)

    def translated(self, x=0.0, y=0.0, z=0.0) -> Atom:
        return attr.evolve(self, x=self.x + x, y=self.y + y, z=self.z + z)

    def scaled(self, factor: float) -> Atom:
        return attr.evolve(self, x=factor * self.x, y=factor * self.y, z=factor * self.z)

    def rotated(self, r: np.ndarray) -> Atom:
        if not r.shape == (3, 3):
            raise ValueError(f"Rotation matrix R must be 3x3, got {r.shape}")
        x, y, z = r @ np.array([self.x, self.y, self.z])
        return attr.evolve(self, x=x, y=y, z=z)

    def flipped_x(self) -> Atom:
        return attr.evolve(self, x=-self.x)

    def distance(self, other: Atom) -> float:
        return sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2)

    @property
    def normalized_coords(self) -> np.array:
        return self.coords / np.linalg.norm(self.coords)

    def rotation_matrix_to(self, x: float, y: float, z: float):
        """https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d"""
        a = self.normalized_coords
        b = np.array([x, y, z])
        b = b / np.linalg.norm(b)
        v = np.cross(a, b)
        s = np.linalg.norm(v)
        c = np.dot(a, b)
        I = np.eye(*a.shape)  # noqa: E741

        v_x = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return I + v_x + v_x @ v_x * (1 - c) / s ** 2

    def angle_to_xy_plane(self) -> float:
        return angle(self.coords, [self.x, self.y, 0])

    def angle_to_xz_plane(self) -> float:
        return angle(self.coords, [self.x, 0, self.z])

    def angle_to_yz_plane(self) -> float:
        return angle(self.coords, [0, self.y, self.z])

    def rotated_about_x(self, angle: float) -> "Atom":
        r = np.array([[1, 0, 0], [0, cos(angle), -sin(angle)], [0, sin(angle), cos(angle)]])
        return self.with_coords(*r @ self.coords)

    def rotated_about_z(self, angle: float) -> "Atom":
        r = np.array([[cos(angle), -sin(angle), 0], [sin(angle), cos(angle), 0], [0, 0, 1]])
        return self.with_coords(*r @ self.coords)

    def rotated_about_y(self, angle: float) -> "Atom":
        r = np.array([[cos(angle), 0, -sin(angle)], [0, 1, 0], [sin(angle), 0, cos(angle)]])
        return self.with_coords(*r @ self.coords)

    def __eq__(self, other: "Atom") -> bool:
        return (
            self.element is other.element
            and isclose(self.x, other.x, abs_tol=conf.small_number)
            and isclose(self.y, other.y, abs_tol=conf.small_number)
            and isclose(self.z, other.z, abs_tol=conf.small_number)
        )
