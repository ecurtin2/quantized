from textwrap import dedent
from functools import wraps

import attr

DESC_KEY = "QUANTIZED_FIELD_DESC"


def attrs(*args, **kwargs):
    @wraps(attr.s)
    def wrapper(cls):
        new_cls = attr.s(*args, **kwargs)(cls)
        new_cls.__doc__ = parse_attr(cls)
        return new_cls

    return wrapper


def document_me(f):
    f.include_in_docs = True
    return f


def repr_validator(validator) -> str:
    if isinstance(validator, attr.validators._OptionalValidator):
        return repr(validator.validator)
    else:
        return repr(validator)


def parse_attr(cls):
    if not attr.has(cls):
        raise TypeError(f"{cls} is not an attrs class")

    fields = []
    for f in attr.fields(cls):
        if hasattr(f.type, "__args__"):
            t = " | ".join(str(i.__name__) for i in f.type.__args__)
        elif hasattr(f.type, "__name__"):
            t = f.type.__name__
        else:
            t = f.type

        fields.append(f" **{f.name} ({t})** {get_desc(f)}")
        if f.default is not attr.NOTHING:
            fields.append(f"> Default: {f.default}")
        if f.validator is not None:
            try:
                v = ", ".join(
                    map(repr_validator, f.validator._validators)
                )  # and validator - multiple validators
            except AttributeError:
                v = repr_validator(f.validator)
            fields.append(f"> Constraints:  {v}")

    if cls.__doc__ is None:
        docstring = [""]
    else:
        docstring = cls.__doc__.splitlines()
    short_desc, *rest = docstring
    fields = f"\n\n#### Fields" + "\n\n" + "\n\n".join(fields) + "\n\n"
    return short_desc + dedent("\n".join(rest)) + fields


def attrib(*args, desc: str = None, **kwargs) -> attr.Attribute:
    meta = kwargs.pop("metadata", {})
    new_meta = {DESC_KEY: desc, **meta}
    return attr.ib(*args, metadata=new_meta, **kwargs)


def get_desc(a: attr.Attribute) -> str:
    desc = a.metadata.get(DESC_KEY, "")
    if desc is None:
        return ""
    else:
        return desc
