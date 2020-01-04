from typing import Tuple, Union

import attr
import numpy as np

from transit_chem.config import conf
from transit_chem.utils import Parabola
from transit_chem.validation import Range


@attr.s(frozen=True)
class TripleWell:
    center1: Tuple[float, float] = attr.ib()
    barrier12: Tuple[float, float] = attr.ib()
    center2: Tuple[float, float] = attr.ib()
    barrier23: Tuple[float, float] = attr.ib()
    center3: Tuple[float, float] = attr.ib()
    well1: Parabola = attr.ib()
    well2: Parabola = attr.ib()
    well3: Parabola = attr.ib()
    """Construct a well potential from the points.

    x                                                         x
    x                                                         x
    x                                                         x
    x                                                        x
    x                                                        x
     x          barrier12                                   x
     x         x         x            barrier23            x
      x       x          x           x         x          x
       center1            x         x           x       x
                           x       x             center3
                            center2
    """

    def __attrs_post_init__(self):

        # Wells are in expected left to right ordering.
        if (
            not self.center1[0]
            < self.barrier12[0]
            < self.center2[0]
            < self.barrier23[0]
            < self.center3[0]
        ):
            raise ValueError("Points are not in ascending x-value order.")

        # Wells are below the barriers
        assert self.center1[1] < self.barrier12[1]
        assert self.center2[1] < self.barrier12[1]
        assert self.center2[1] < self.barrier23[1]
        assert self.center3[1] < self.barrier23[1]

    def _call_numpy(self, x: np.ndarray) -> np.ndarray:
        y = np.zeros_like(x)

        mask1 = np.where(x < self.barrier12[0])
        y[mask1] = self.well1(x[mask1])

        mask2 = np.where((self.barrier12[0] <= x) & (x <= self.barrier23[0]))
        y[mask2] = self.well2(x[mask2])

        mask3 = np.where(x > self.barrier23[0])
        y[mask3] = self.well3(x[mask3])

        return y

    def _call_scalar(self, x: float) -> float:
        if x < self.barrier12[0]:
            return self.well1(x)
        elif self.barrier12[0] <= x <= self.barrier23[0]:
            return self.well2(x)
        elif x > self.barrier23[0]:
            return self.well3(x)
        else:
            raise ValueError(
                f"Value {x} is not valid for {self}. "
                "Probably the well/barrier parameters are invalid"
            )

    def __call__(self, x: Union[np.ndarray, float]) -> Union[np.ndarray, float]:
        if isinstance(x, np.ndarray):
            return self._call_numpy(x)
        else:
            return self._call_scalar(x)

    @staticmethod
    def from_params(
        well1_depth: float,
        well1_halfwidth: float,
        bridge_length: float,
        bridge_depth: float,
        well3_halfwidth: float,
        well3_depth: float,
    ):
        center1 = (0, 0)
        barrier1 = (center1[0] + well1_halfwidth, center1[1] + well1_depth)
        center2 = (barrier1[0] + 0.5 * bridge_length, barrier1[1] - bridge_depth)
        barrier2 = (barrier1[0] + bridge_length, barrier1[1])
        center3 = (barrier2[0] + well3_halfwidth, barrier2[1] - well3_depth)

        def fit_well(center, barrier):
            center_x, center_y = center
            barrier_x, barrier_y = barrier
            # Third point is reflecting barrier about center
            x = -barrier_x + 2 * center_x
            y = barrier_y
            return Parabola.from_points(center, barrier, (x, y))

        well1 = fit_well(center1, barrier1)
        well2 = Parabola.from_points(barrier1, center2, barrier2)
        well3 = fit_well(center3, barrier2)

        return TripleWell(center1, barrier1, center2, barrier2, center3, well1, well2, well3)


@attr.s(frozen=True)
class Harmonic:
    center: float = attr.ib(validator=[Range(-conf.large_number, conf.large_number)])
    mass: float = attr.ib(default=1.0, validator=Range(conf.small_number, conf.large_number))
    omega: float = attr.ib(default=1.0, validator=Range(conf.small_number, conf.large_number))

    def __call__(self, x: float) -> float:
        return 0.5 * self.mass * (self.omega ** 2) * ((x - self.center) ** 2)
