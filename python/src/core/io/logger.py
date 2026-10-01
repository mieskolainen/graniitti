# Shared logging helpers for icetune-style command-line tools
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import logging
import socket
import sys
from typing import TextIO

DEFAULT_LOG_FORMAT = "[%(hostname)s:%(asctime)s] %(levelname)s %(name)s: %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


class _ContextFilter(logging.Filter):
    """Inject shared logging context into every record."""

    # Cache immutable host context in each configured handler
    def __init__(self):
        super().__init__()
        self.hostname = hostname()

    # Add host context to one log record
    def filter(self, record: logging.LogRecord) -> bool:
        record.hostname = self.hostname
        return True


def configure(*, level: int | str = logging.INFO, stream: TextIO | None = None, force: bool = False) -> logging.Logger:
    """Configure root logging with the shared hostname and timestamp format"""
    if isinstance(level, str):
        level = logging._nameToLevel.get(level.upper(), logging.INFO)
    handler = logging.StreamHandler(stream or sys.stdout)
    handler.setFormatter(logging.Formatter(DEFAULT_LOG_FORMAT, datefmt=DEFAULT_DATE_FORMAT))
    handler.addFilter(_ContextFilter())
    logging.basicConfig(level=level, handlers=[handler], force=force)
    logger = logging.getLogger()
    logger.setLevel(level)
    return logger


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a standard logger without changing process-wide configuration"""
    return logging.getLogger(name)


def level_from_verbosity(verbose: int | None) -> int:
    if verbose is None:
        return logging.INFO
    return logging.WARNING if verbose <= 0 else logging.INFO if verbose == 1 else logging.DEBUG


def hostname() -> str:
    """Return the short local hostname used in log records."""
    return socket.gethostname().split(".")[0]
