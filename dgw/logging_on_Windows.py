"""Set up logging like we do in DIALS, to test behaviour on Windows"""

import procrunner
import logging.config
import os
import sys
import time

try:
    from colorlog import ColoredFormatter
except ImportError:
    ColoredFormatter = None

# https://stackoverflow.com/questions/25194864/python-logging-time-since-start-of-program/25196134#25196134
class DialsLogfileFormatter:
    """A formatter for log files that prepends messages with the elapsed time
    or messages at warning level or above with 'WARNING:'"""

    def __init__(self, timed):
        self.timed = timed
        self.start_time = time.time()
        self.prefix = ""

    def format(self, record):
        if self.timed:
            elapsed_seconds = record.created - self.start_time
            prefix = "{:6.1f}: ".format(elapsed_seconds)
        else:
            prefix = ""
        indent = len(prefix)
        msg = record.getMessage()

        if record.levelno >= logging.WARNING:
            prefix = "{prefix:>{indent}s}".format(indent=indent, prefix="WARN: ")

        msg = msg.replace("\n", "\n" + " " * indent)
        if prefix == self.prefix:
            return " " * indent + msg
        else:
            self.prefix = prefix
            return prefix + msg


def config(verbosity=0, logfile=None):
    """
    Configure the logging.

    :param verbosity: Verbosity level of log output. Possible values:
                        * 0: Info log output to stdout/logfile
                        * 1: Info & debug log output to stdout/logfile
    :type verbosity: int
    :param logfile: Filename for log output.  If False, no log file is written.
    :type logfile: str
    """

    console = logging.StreamHandler(sys.stdout)
    if (
        "NO_COLOR" not in os.environ
        and sys.stdout.isatty()
        and ColoredFormatter is not None
    ):
        color_formatter = ColoredFormatter(
            "%(log_color)s%(message)s",
            log_colors={
                "DEBUG": "blue",
                "WARNING": "yellow",
                "ERROR": "red",
                "CRITICAL": "red,bg_white",
            },
        )
        console.setFormatter(color_formatter)

    dials_logger = logging.getLogger("dials")
    dials_logger.addHandler(console)

    logging.captureWarnings(True)
    warning_logger = logging.getLogger("py.warnings")
    warning_logger.addHandler(console)

    if verbosity > 1:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    if logfile:
        fh = logging.FileHandler(filename=logfile, mode="w", encoding="utf-8")
        fh.setLevel(loglevel)
        fh.setFormatter(DialsLogfileFormatter(timed=verbosity))
        dials_logger.addHandler(fh)
        warning_logger.addHandler(fh)

    dials_logger.setLevel(loglevel)
    #   logging.getLogger("dxtbx").setLevel(logging.DEBUG)
    console.setLevel(loglevel)


def log_out():
    config(logfile="foo.log")

    logger = logging.getLogger("dials")
    logger.info("Hello")
    logger.warning("Watch out!")  # Colour
    logger.info("Å σ")


if __name__ == "__main__":

    # This works:
    log_out()

    # This does not work on Windows, but does on Linux:
    procrunner.run(
        ["dials.python", "-c", "from logging_on_Windows import log_out; log_out()"]
    )
