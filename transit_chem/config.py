import json
from pathlib import Path
from shutil import rmtree
from typing import Union

import attr
import cattr

default_conf_path: Path = Path.home() / ".transit_chem/conf.yaml"


def to_bool(x: Union[str, bool]) -> bool:
    if isinstance(x, str):
        return x.lower() in ("true", "yes")
    elif isinstance(x, bool):
        return x
    else:
        raise TypeError(f"Variable {x} of type {type(x)} can't be converted to bool")


@attr.s()
class Config:
    harmonic_oscillator_max_n: int = attr.ib(default=50, converter=int)
    small_number: float = attr.ib(default=1e-8, converter=float)
    large_number: float = attr.ib(default=1000.0, converter=float)
    float_tol: float = attr.ib(default=1e-6, converter=float)
    enable_progressbar: bool = attr.ib(default=False, converter=to_bool)
    cache_dir: Path = attr.ib(default=Path.home() / ".transit_chem/cache", converter=Path)
    joblib_verbosity: int = attr.ib(default=0, converter=int)


cattr.register_structure_hook(Path, lambda s, t: Path(s))
cattr.register_unstructure_hook(Path, lambda p: str(p))
default_conf = Config()


def load(p: Path = default_conf_path) -> Config:

    conf_file_d = json.loads(p.read_text())
    return attr.evolve(default_conf, **conf_file_d)


try:
    conf: Config = load()
except FileNotFoundError:
    conf: Config = Config()


def set_items(path: Path = default_conf_path, **kwargs):
    new_conf = attr.evolve(conf, **kwargs)
    set_conf(new_conf, path=path)


def show(c: Config = conf):
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
    set_conf(c=default_conf, path=path)
    show(default_conf)


def set_conf(c: Config, path: Path = default_conf_path):
    d = cattr.unstructure(c)
    path.write_text(json.dumps(d))
    show(c)
    print(f"Updating config in {path}")


def clear_cache(c: Config = conf):
    rmtree(c.cache_dir)
    print(f"Cleared {c.cache_dir}")
