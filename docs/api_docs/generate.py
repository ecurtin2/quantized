from typing import Dict, List, Any

from quantized import elements


def markdown_table(d: Dict[str, List[Any]], width: int = 10, none: str = "-") -> str:
    cstr = lambda x: str(x).center(width) if x is not None else none.center(width)

    header = "|" + "|".join(map(cstr, d)) + "|"
    br = "|" + "|".join("-" * width for _ in d) + "|"
    rows = ["|" + "|".join(map(cstr, row)) + "|" for row in zip(*d.values())]
    return "\n".join([header, br] + rows)


attributes = ["name", "z", "voie", "alpha"]
d = {a: [getattr(e, a) for e in elements.all_elements] for a in attributes}
print(markdown_table(d))
