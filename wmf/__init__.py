from lmfdb.utils import make_logger
wmf_logger = make_logger('wdb')
try:
    import colorlog
    from colorlog import ColoredFormatter
    LOGFORMAT = "  %(log_color)s%(levelname)-10s%(filename)s:%(lineno)d%(reset)s | %(log_color)s%(message)s%(reset)s"
    formatter = ColoredFormatter(LOGFORMAT)
except:
    LOGFORMAT = "  %(levelname)-10s%(filename)s:%(lineno)d | %(message)s"
    formatter = logging.Formatter(LOGFORMAT)
wmf_logger.handlers[0].setFormatter(formatter)
from web_newforms_computing import WebNewForm_computing
from web_modform_space_computing import WebModFormSpace_computing

from generate_web_data import generate_web_modform_spaces,dimension_from_db,web_modformspace_in_db,update_dimension_tables,generate_dimension_tables
