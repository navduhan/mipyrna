"""miPyRNA Logging Module

This module provides the logger factory shared by all miPyRNA modules. Loggers
write to the console, and to a log file once :func:`set_log_file` has been called.
Earlier messages are held in memory and written then.

Features:
    - Named per-module loggers that never propagate to the root logger
    - Handlers are attached only once per logger name
    - One shared log file, created only when a run directory is set
    - Console logging at DEBUG level with a ``HH:MM:SS`` timestamp format

Configuration:
    The logger name is passed to :func:`MiPyRNALogger`. The log file location is set
    with :func:`set_log_file`; until then messages are only held in memory, so commands such
    as ``--help`` and ``--check`` leave no files behind.

Dependencies:
    - logging: Python standard library

Classes / Functions / Exceptions:
    - MiPyRNALogger: Return a configured logger for a module.
    - set_log_file: Start writing all miPyRNA logs to a file.

:Author: Naveen Duhan
"""

import logging
import os

FORMAT = '[%(asctime)s]  %(module)s :: %(levelname)s : %(message)s'
DATEFMT = '%H:%M:%S'


class _SharedFileHandler(logging.FileHandler):
    """File handler that stays silent and creates no file until a path is set."""

    def __init__(self):
        super().__init__(os.devnull, mode='a', delay=True)
        self.enabled = False
        self.pending = []

    def emit(self, record):
        """Write the record, or hold it until a log file has been set.

        Args:
            record (logging.LogRecord): The record to write.
        """
        if self.enabled:
            super().emit(record)
        else:
            self.pending.append(record)


_file_handler = _SharedFileHandler()
_file_handler.setFormatter(logging.Formatter(FORMAT, datefmt=DATEFMT))


def set_log_file(path):
    """Write all miPyRNA logs to ``path`` from now on.

    The parent directory is created if needed and the file is appended to. Messages
    logged earlier in the run are written first.

    Args:
        path (str): Log file path.
    """
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    _file_handler.close()
    _file_handler.baseFilename = path
    _file_handler.enabled = True
    for record in _file_handler.pending:
        _file_handler.emit(record)
    _file_handler.pending.clear()


def MiPyRNALogger(mode, log):
    """Create or fetch the logger for a miPyRNA module.

    The logger writes to the console and to the shared log file set by
    :func:`set_log_file`. If the named logger already has handlers it is returned
    unchanged.

    Args:
        mode (str): Unused, kept for backward compatibility.
        log (str): Logger name, normally the module name.

    Returns:
        logging.Logger: The configured logger.
    """
    if logging.getLoggerClass() is not logging.Logger:
        logging.setLoggerClass(logging.Logger)
    existing = logging.Logger.manager.loggerDict.get(log)
    if isinstance(existing, logging.Logger) and type(existing) is not logging.Logger:
        del logging.Logger.manager.loggerDict[log]
    logger = logging.getLogger(log)
    logger.propagate = False
    if logger.handlers:
        return logger

    logger.setLevel(logging.DEBUG)

    logger.addHandler(_file_handler)

    consoleHandler = logging.StreamHandler()

    consoleHandler.setFormatter(logging.Formatter(FORMAT, datefmt=DATEFMT))

    consoleHandler.setLevel(logging.DEBUG)

    logger.addHandler(consoleHandler)

    return logger
