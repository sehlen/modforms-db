from lmfdb.utils import make_logger
import logging
wmf_logger = make_logger('wdb')
try:
    import colorlog
    from colorlog import ColoredFormatter
    LOGFORMAT = "  %(log_color)s%(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(LOGFORMAT)
except:
    LOGFORMAT = "  %(levelname)-10s%(filename)s:%(lineno)d | %(message)s"
    formatter = logging.Formatter(LOGFORMAT)
stream = logging.StreamHandler()
stream.setLevel(LOG_LEVEL)
stream.setFormatter(formatter)
if len(wmf_logger.handlers)>0:
    wmf_logger.handlers = []
wmf_logger.handlers.append(stream)
wmf_logger.propagate=False
from web_newforms_computing import WebNewForm_computing
from web_modform_space_computing import WebModFormSpace_computing

from generate_web_data import generate_web_modform_spaces,dimension_from_db,web_modformspace_in_db,update_dimension_tables,generate_dimension_tables

from utils import orbit_label,orbit_index_from_label
