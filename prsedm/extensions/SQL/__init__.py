"""
SQL utilities for PRSedm — includes helper functions for locating/downloading the
variants.db file used for diabetes PRS scoring.
"""

from .get_dm_sql import get_snp_db, get_dm_sql

__all__ = ["get_snp_db", "get_dm_sql"]
