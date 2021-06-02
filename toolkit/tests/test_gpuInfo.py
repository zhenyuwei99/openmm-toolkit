#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: test_gpuInfo.py
created time : 2021/06/02
last edit time : 2021/06/02
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

from ..utils.gpuInfo import *

def test_num_gpu():
    assert num_gpus() == 1