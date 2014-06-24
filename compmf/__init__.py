import logging
LOG_LEVEL = logging.DEBUG
logging.root.setLevel(LOG_LEVEL)
try:
    import colorlog
    from colorlog import ColoredFormatter
    LOGFORMAT = "  %(log_color)s%(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(LOGFORMAT)
except:
    LOGFORMAT = "  %(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = logging.Formatter(LOGFORMAT)
stream = logging.StreamHandler()
stream.setLevel(LOG_LEVEL)
stream.setFormatter(formatter)
clogger = logging.getLogger(__name__)
if not clogger.handlers:
    clogger.addHandler(stream)
    
clogger.setLevel(logging.DEBUG)

import character_conversions
import filesdb
import compute
import mongodb
from compute import ComputeMFData
from filesdb import FilenamesMFDB, FilenamesMFDBLoading
from mongodb import CompMF
