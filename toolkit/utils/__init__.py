__author__ = "Zhenyu Wei"
__maintainer__ = "Zhenyu Wei" 
__email__ = "zhenyuwei99@gmail.com"
__copyright__ = "Copyright 2021-2021, Southeast University and Zhenyu Wei"
__license__ = "GPLv3"

from .gpuInfo import num_gpus
from .check import check_path, check_quantity

__all__ = [
    'num_gpus',
    'check_path', 'check_quantity'
]