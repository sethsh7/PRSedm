# scoreprs/__init__.py
__version__ = '1.0.1a'
__author__ = 'Seth A. Sharp'
__email__ = 'ssharp@stanford.edu'

from .core import score_bcf, get_snp_db
from .extensions import gen_dm

__all__ = ['score_bcf', 'gen_dm', 'get_snp_db']
