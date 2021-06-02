import nvgpu

from ..exceptions import *

def num_gpus():
    gpu_list = nvgpu.gpu_info()
    return len(gpu_list)