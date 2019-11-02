import attr


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
            raise ValueError(
                f"Min must be below max got min={self.min}, max={self.max}"
            )

    def __call__(self, instance, attribute, value):
        if not self.min <= value <= self.max:
            raise ValueError(
                f"Validation error for {instance}. "
                f"{attribute.name} out of range, must be "
                f"between {self.min} and {self.max}, got {value}"
            )
