# scoreprs/__init__.py
__version__ = "1.1.0"
__author__ = "Seth A. Sharp"
__email__ = "ssharp@stanford.edu"

from .core import score_bcf
from .extensions import gen_dm
from .extensions.SQL import get_snp_db, get_dm_sql

__all__ = ["score_bcf", "gen_dm", "get_snp_db", "get_dm_sql"]
