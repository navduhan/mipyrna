"""miPyRNA Slurm Array Worker

Runs one task of a Slurm job array submitted by :func:`mipyrna.utility.run_array`.
It is an internal entry point and is not listed in ``mipyrna --help``.

Configuration:
    ``python -m mipyrna._worker <tasks.pkl>``; the task is chosen by ``SLURM_ARRAY_TASK_ID``.

Dependencies:
    - pickle: Reads the task list written by the submitting process

Classes / Functions / Exceptions:
    - main: Run one task.

:Author: Naveen Duhan
"""

import os
import pickle
import sys


def main(argv=None):
    """Run the task of this array element.

    Args:
        argv (list, optional): Arguments; the first is the task list file. Defaults to ``sys.argv[1:]``.

    Returns:
        int: Process exit code.
    """
    argv = sys.argv[1:] if argv is None else argv
    with open(argv[0], "rb") as fh:
        tasks = pickle.load(fh)
    obj, method, args = tasks[int(os.environ["SLURM_ARRAY_TASK_ID"])]
    getattr(obj, method)(*args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
