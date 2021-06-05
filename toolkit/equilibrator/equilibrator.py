#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: equilibrator.py
created time : 2021/06/02
last edit time : 2021/06/05
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import os, sys, datetime
import simtk.openmm.app as app
import simtk.openmm.openmm as openmm
import simtk.unit as unit
from .. import AVAILABLE_PLATFORM
from ..utils import check_quantity
from ..exceptions import *

class Equilibrator:
    def __init__(self, out_dir, cut_off, pdb_file, out_prefix, platform) -> None:
        # Read input
        if pdb_file == '':
            self._pdb = None
        else:
            self._pdb = app.PDBFile(pdb_file)
        self._cut_off = check_quantity(cut_off, unit.angstrom)
        if not platform in AVAILABLE_PLATFORM:
            raise InvalidPlatformError(
                '%s is an invalid platform name. Support list:\n%s' 
                %(AVAILABLE_PLATFORM)
            )
        else:
            self._platform = openmm.Platform_getPlatformByName(platform)

        # Check and define path
        self._out_log_dir = os.path.join(out_dir, 'log_files')
        self._out_pdb_dir = os.path.join(out_dir, 'pdb_files')

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
            os.mkdir(self._out_log_dir)
            os.mkdir(self._out_pdb_dir)
        else:
            
            if not os.path.exists(self._out_log_dir):
                os.mkdir(self._out_log_dir)
            if not os.path.exists(self._out_pdb_dir):
                os.mkdir(self._out_pdb_dir)

        self._out_prefix = out_prefix
        self._out_log_file_path = os.path.join(self._out_log_dir, out_prefix + '.log')
        self._out_pdb_file_path = os.path.join(self._out_pdb_dir, out_prefix + '.pdb')
        self._out_pdb_restart_file_path = os.path.join(self._out_pdb_dir, out_prefix + '_restart.pdb')
        
        # Predefined attribute
        self._simulation = None

    def equilibrate(self):
        start_time = datetime.datetime.now().replace(microsecond=0)
        self._log_file = open(self._out_log_file_path, 'w')
        self._setup()
        self._execute()
        self._teardown()
        end_time = datetime.datetime.now().replace(microsecond=0)
        print('Total running time:', end='\t', file=self._log_file)
        print(end_time - start_time, file=self._log_file)    
        self._log_file.close()
    
    def _setup(self):
        pass

    def _execute(self):
        pass

    def _teardown(self):
        pass

    def _check_simulation(self):
        if self._simulation == None:
            raise UndefinedSimulationInstanceError(
                'self._simulation is not defined, check self.setup method'
            )
        else:
            return 0

    @property
    def pdb_restart_file(self):
        return self._out_pdb_restart_file_path

    @property
    def out_prefix(self):
        return self._out_prefix

    @out_prefix.setter
    def out_prefix(self, out_prefix):
        self._out_prefix = out_prefix
        self._out_log_file_path = os.path.join(self._out_log_dir, out_prefix + '.log')
        self._out_pdb_file_path = os.path.join(self._out_pdb_dir, out_prefix + '.pdb')
        self._out_pdb_restart_file_path = os.path.join(self._out_pdb_dir, out_prefix + '_restart.pdb')
