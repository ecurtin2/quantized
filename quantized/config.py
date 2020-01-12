import json
from pathlib import Path
from shutil import rmtree
from typing import Union

from quantized.attr_wrapped import attrs, attrib
import attr
import cattr


default_conf_path: Path = Path.home() / ".quantized/conf.yaml"


__all__ = ["to_bool", "Config", "conf"]


def to_bool(x: Union[str, bool]) -> bool:
    """Convert x to a bool

    For string arguments, "true" and "yes" are considered
    Truthy, and are compared in a case-insensitive manner.

    """
    if isinstance(x, str):
        return x.lower() in ("true", "yes")
    elif isinstance(x, bool):
        return x
    else:
        raise TypeError(f"Variable {x} of type {type(x)} can't be converted to bool")


@attrs(frozen=True)
class Config:
    """Class to hold configuration for quantized"""

    harmonic_oscillator_max_n: int = attrib(
        default=50, converter=int, desc="Limit of n for the harmonic oscillator"
    )
    small_number: float = attrib(
        default=1e-8, converter=float, desc="A value to be used to ensure quantities are nonzero"
    )
    large_number: float = attrib(
        default=1000.0,
        converter=float,
        desc="A cutoff value that's used to guard against errors, where a number is not normally expected to be large.",
    )
    float_tol: float = attrib(
        default=1e-6,
        converter=float,
        desc="Two floats are considered equal if they are within this value.",
    )
    enable_progressbar: bool = attrib(
        default=False,
        converter=to_bool,
        desc="If True, certain methods will print a progress bar to the screen",
    )
    cache_dir: Path = attrib(
        default=".quantized/cache",
        converter=Path,
        desc="The directory where cached objects are stored",
    )
    joblib_verbosity: int = attrib(
        default=0, converter=int, desc="Verbosity level for joblib's cache"
    )


cattr.register_structure_hook(Path, lambda s, t: Path(s))
cattr.register_unstructure_hook(Path, lambda p: str(p))
default_conf = Config()


def load(p: Path = default_conf_path) -> Config:
    """Load a configuration from a json file"""

    conf_file_d = json.loads(p.read_text())
    return attr.evolve(default_conf, **conf_file_d)


try:
    conf: Config = load()
except FileNotFoundError:
    conf: Config = Config()


def set_items(path: Path = default_conf_path, **kwargs):
    """Set the configuration """
    new_conf = attr.evolve(conf, **kwargs)
    set_conf(new_conf, path=path)


def show(c: Config = conf):
    """Print the current config to stdout"""
    name_width = 30
    value_width = 40
    type_width = 10

    print(
        "Name".ljust(name_width)
        + "Value".ljust(value_width)
        + "Type".ljust(type_width)
        + "Default".ljust(value_width)
    )
    total_width = name_width + value_width * 2 + type_width
    print("-" * total_width)
    for field in attr.fields(Config):
        print(
            field.name.ljust(name_width)
            + str(getattr(c, field.name)).ljust(value_width)
            + field.type.__name__.ljust(type_width)
            + str(field.default).ljust(value_width)
        )


def restore_defaults(path: Path = default_conf_path):
    """Return the config to the defaults"""
    set_conf(c=default_conf, path=path)
    show(default_conf)


def set_conf(c: Config, path: Path = default_conf_path):
    """Set the config to the given one"""
    d = cattr.unstructure(c)
    path.write_text(json.dumps(d))
    show(c)
    print(f"Updating config in {path}")


def clear_cache(c: Config = conf):
    """Remove all objects from the cache"""
    rmtree(c.cache_dir)
    print(f"Cleared {c.cache_dir}")
