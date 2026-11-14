"""miPyRNA package."""

import collections
import collections.abc as abc

__version__ = "1.0.0"

for _name in ("Callable", "Iterable", "Mapping", "MutableMapping", "MutableSet", "Sequence"):
    if not hasattr(collections, _name):
        setattr(collections, _name, getattr(abc, _name))
