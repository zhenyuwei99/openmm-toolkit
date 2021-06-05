#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: exceptions.py
created time : 2021/06/02
last edit time : 2021/06/02
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

class NotFoundCUDAError(Exception):
    '''
    NotFoundCUDAError related to error when CUDA driver is invalid 
    Used in:
    - utils.gpuinfo
    '''    
    pass

class UndefinedSimulationInstanceError(Exception):
    ''' 
    UndefinedSimulationInstanceError related to error when self._simulation is not defined.
    Used in:
    - equilibrator.equilibrator
    '''
    pass