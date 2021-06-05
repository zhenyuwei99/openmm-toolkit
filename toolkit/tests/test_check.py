#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: test_check.py
created time : 2021/06/05
last edit time : 2021/06/05
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import simtk.unit as unit
from ..utils import *

def test_check_quantity():
    assert check_quantity(1, unit.kelvin) == 1 * unit.kelvin
    assert check_quantity(1 * unit.kelvin, unit.kelvin) == 1 * unit.kelvin