from math import isinf, isnan
from typing import Any, Callable

from quantized.attr_wrapped import attrs, attrib
from attr import Attribute

Validator = Callable[[Any, Attribute], Any]


__all__ = ["Range", "not_inf", "not_nan"]


class ValidationError(ValueError):
    pass


@attrs()
class Range:
    """Attrs Validator to ensure values between a min and a max value.

    """

    min = attrib(desc="The minimum allowed value")
    max = attrib(desc="The maximum allowed value")

    @min.validator
    def check(self, attribute, value):
        if not self.min <= self.max:
            raise ValidationError(f"Min must be below max got min={self.min}, max={self.max}")

    def __call__(self, instance: Any, attribute: Attribute, value: Any) -> Any:
        if not self.min <= value <= self.max:
            raise ValidationError(
                f"Validation error for {instance}. "
                f"{attribute.name} out of range, must be "
                f"between {self.min} and {self.max}, got {value}"
            )


def not_nan(instance: Any, attribute: Attribute, value: float):
    if isnan(value):
        raise ValidationError(
            f"Validation error for {instance}. {attribute.name} is NaN but that's not allowed"
        )


def not_inf(instance: Any, attribute: Attribute, value: float):
    if isinf(value):
        raise ValidationError(
            f"Validation error for {instance}. {attribute.name} is Infinite but that's not allowed"
        )


class _Positive:
    def __call__(self, instance: Any, attribute: Attribute, value: float):
        if not value > 0:
            raise ValidationError(
                f"Validation error for {instance}. {attribute.name} is not positive"
            )

    def __repr__(self):
        return "Positive"


positive = _Positive()
