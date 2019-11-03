from math import isnan, isinf

import attr


__all__ = ["Range", "not_inf", "not_nan"]


class ValidationError(ValueError):
    pass


@attr.s
class Range:
    """Attrs Validator to ensure values between a min and a max value.

    Parameters
    -----------
    min
        The minimum allowed value
    max
        The maximum allowed value

    """

    min = attr.ib()
    max = attr.ib()

    @min.validator
    def check(self, attribute, value):
        if not self.min <= self.max:
            raise ValidationError(
                f"Min must be below max got min={self.min}, max={self.max}"
            )

    def __call__(self, instance, attribute, value):
        if not self.min <= value <= self.max:
            raise ValidationError(
                f"Validation error for {instance}. "
                f"{attribute.name} out of range, must be "
                f"between {self.min} and {self.max}, got {value}"
            )


def not_nan(instance, attribute, value):
    if isnan(value):
        raise ValidationError(
            f"Validation error for {instance}. {attribute.name} is NaN but that's not allowed"
        )


def not_inf(instance, attribute, value):
    if isinf(value):
        raise ValidationError(
            f"Validation error for {instance}. {attribute.name} is Infinite but that's not allowed"
        )
