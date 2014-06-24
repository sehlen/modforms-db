import logging
try:
    import colorlog
    from colorlog import ColoredFormatter
    LOG_LEVEL = logging.DEBUG
    logging.root.setLevel(LOG_LEVEL)
    LOGFORMAT = "  %(log_color)s%(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(LOGFORMAT)
    stream = logging.StreamHandler()
    stream.setLevel(LOG_LEVEL)
    stream.setFormatter(formatter)
    clogger = logging.getLogger(__name__)
    if not clogger.handlers:
        clogger.addHandler(stream)
except:
    clogger = logging.getLogger(__name__)

clogger.setLevel(logging.DEBUG)

import character_conversions
import filesdb
import compute
import mongodb
from compute import ComputeMFData
from filesdb import FilenamesMFDB, FilenamesMFDBLoading
from mongodb import CompMF
