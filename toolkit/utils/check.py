#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: check.py
created time : 2021/06/04
last edit time : 2021/06/05
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import os
import simtk.unit as unit

from ..exceptions import *


def check_path(path):
    if not os.path.exists(path):
        raise InvalidPathError(
            '%s is an invalid path, please check!'
        )
    else:
        return 0

def check_quantity(quantity, target_unit: unit.Unit):
    if unit.is_quantity(quantity):
        return quantity.in_units_of(target_unit)
    else:
        return quantity * target_unit