#!/usr/bin/env python3

"""miPyRNA Command-Line Entry Point

Allows the package to be run as ``python -m mipyrna``. Delegates to
:func:`mipyrna.cli.main` and exits with its return code.

:Author: Naveen Duhan
"""

from mipyrna.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
