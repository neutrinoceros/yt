"""
Logging facility for yt
Will initialize everything, and associate one with each module



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import logging
import sys
from yt.config import ytcfg
from rich.logging import RichHandler
from rich.console import Console
# This next bit is grabbed from:
# http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored

def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[0].levelno
        if(levelno >= 50):
            color = '\x1b[31m'  # red
        elif(levelno >= 40):
            color = '\x1b[31m'  # red
        elif(levelno >= 30):
            color = '\x1b[33m'  # yellow
        elif(levelno >= 20):
            color = '\x1b[32m'  # green
        elif(levelno >= 10):
            color = '\x1b[35m'  # pink
        else:
            color = '\x1b[0m'  # normal
        ln = color + args[0].levelname + '\x1b[0m'
        args[0].levelname = ln
        return fn(*args)
    return new

level = min(max(ytcfg.getint("yt", "loglevel"), 0), 50)
ufstring = "%(name)-3s: [%(levelname)-9s] %(asctime)s %(message)s"
cfstring = "%(name)-3s: [%(levelname)-18s] %(asctime)s %(message)s"

if ytcfg.getboolean("yt", "stdoutStreamLogging"):
    stream = sys.stdout
else:
    stream = sys.stderr

ytLogger = logging.getLogger("yt")

def disable_stream_logging():
    if len(ytLogger.handlers) > 0:
        ytLogger.removeHandler(ytLogger.handlers[0])
    h = logging.NullHandler()
    ytLogger.addHandler(h)

if ytcfg.getboolean("yt", "suppressStreamLogging"):
    disable_stream_logging()
else:
    handler = RichHandler(console=Console(file=stream, width=120))

    # add the handler to the logger
    ytLogger.addHandler(handler)
    ytLogger.setLevel(level)
    ytLogger.propagate = False

    f = logging.Formatter("%(name)-3s: %(message)s")
    ytLogger.handlers[0].setFormatter(f)

ytLogger.debug("Set log level to %s", level)
