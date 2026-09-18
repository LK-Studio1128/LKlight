"""NumPy backward-compat shim (auto-loaded via sitecustomize).

LightDock 0.9.4 depends on ProDy, which still calls removed NumPy aliases and
module paths (``numpy.alltrue``, ``from numpy.lib.arraysetops import isin``,
etc.). Prepend this directory to PYTHONPATH so Python imports this module at
startup.  ``from numpy.lib.<mod> import ...`` statements require the shim
modules to be registered in ``sys.modules``, not merely attached as package
attributes.
"""
import sys as _sys

try:
    import numpy as _np
    import numpy.lib as _nplib

    # removed function aliases
    for _old, _new in (("alltrue", "all"), ("sometrue", "any"),
                       ("cumproduct", "cumprod"), ("product", "prod")):
        if not hasattr(_np, _old) and hasattr(_np, _new):
            setattr(_np, _old, getattr(_np, _new))

    # numpy 2.x moved lib submodules to *_impl; restore old names in both
    # the package namespace and sys.modules so `from numpy.lib.X import Y` works
    for _mod in ("arraysetops", "function_base", "nanfunctions", "shape_base",
                 "twodim_base", "type_check", "index_tricks", "ufunclike",
                 "arrayterator", "npyio", "polynomial", "scimath"):
        _full = f"numpy.lib.{_mod}"
        if not hasattr(_nplib, _mod) or _full not in _sys.modules:
            for _cand in (f"numpy.lib._{_mod}_impl", _full):
                try:
                    _impl = __import__(_cand, fromlist=[_mod])
                    setattr(_nplib, _mod, _impl)
                    _sys.modules[_full] = _impl
                    break
                except Exception:
                    continue
except Exception:  # pragma: no cover - best effort only
    pass
