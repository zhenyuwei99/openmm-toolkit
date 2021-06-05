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

from ..exceptions import *

class Equilibrator:
    def __init__(self, pdb_file, output_dir, out_prefix='equilibrator') -> None:
        # Read input
        self._pdb = app.PDBFile(pdb_file)

        # Check and define path
        out_log_dir = os.path.join(output_dir, 'log_files')
        out_pdb_dir = os.path.join(output_dir, 'pdb_files')

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
            os.mkdir(out_log_dir)
            os.mkdir(out_pdb_dir)
        else:
            
            if not os.path.exists(out_log_dir):
                os.mkdir(out_log_dir)
            if not os.path.exists(out_pdb_dir):
                os.mkdir(out_pdb_dir)

        self._out_log_file_path = os.path.join(out_log_dir, out_prefix + '.log')
        self._log_file = open(self._out_log_file_path, 'w')
        sys.stdout = self._log_file

        self._out_pdb_file_path = os.path.join(out_pdb_dir, out_prefix + '.pdb')
        self._out_pdb__restart_file_path = os.path.join(out_pdb_dir, out_prefix + '_restart.pdb')
        
        # Predefined attribute
        self._simulation = None

    def equilibrate(self):
        start_time = datetime.datetime.now().replace(microsecond=0)
        self._setup()
        self._execute()
        self._teardown()
        end_time = datetime.datetime.now().replace(microsecond=0)
        print('Total running time:', end='\t')
        print(end_time - start_time)    
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