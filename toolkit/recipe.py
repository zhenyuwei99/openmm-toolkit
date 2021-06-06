#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: recipe.py
created time : 2021/06/05
last edit time : 2021/06/05
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import simtk.openmm.app as app
from toolkit.equilibrator.equilibrator import Equilibrator
from .equilibrator import Equilibrator
from .sampler import Sampler

class Recipe:
    def __init__(self, initial_pdb, jobs) -> None:
        self._jobs = jobs
        self._initial_pdb = app.PDBFile(initial_pdb)

    def appendJobs(self, *jobs):
        for job in jobs:
            self._jobs.append(job)

    def run(self):
        for index, job in enumerate(self._jobs):
            out_index = index + 1
            job.out_prefix = (str(out_index) if out_index > 9 else '0' + str(out_index)) + '_' + job.out_prefix
            if index != 0:
                job._pdb = app.PDBFile(self._jobs[index-1].pdb_restart_file)
            else:
                job._pdb = self._initial_pdb
            if isinstance(job, Equilibrator):
                job.equilibrate()
            elif isinstance(job, Sampler):
                job.sample()