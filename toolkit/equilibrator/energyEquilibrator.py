#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: energyEquilibrator.py
created time : 2021/06/05
last edit time : 2021/06/05
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import simtk.openmm.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as unit
from . import Equilibrator
from ..utils import *

class EnergyEquilibrator(Equilibrator):
    def __init__(
        self, out_dir, cut_off=12, pdb_file='', out_prefix='energy_equilibrator', platform='CUDA'
    ) -> None:
        super().__init__(out_dir, cut_off, pdb_file, out_prefix, platform)

    def _setup(self):
        # Forcefield
        force_field = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml')
        system = force_field.createSystem(
            self._pdb.topology, nonbondedMethod=app.PME, 
            nonbondedCutoff=self._cut_off, constraints=app.HBonds
        )

        # Integrator
        integrator = openmm.LangevinIntegrator(
            0 * unit.kelvin, 0.001/unit.femtosecond, 1 * unit.femtosecond
        )

        # Simulation
        self._simulation = app.Simulation(
            self._pdb.topology, system, integrator, self._platform
        )
        self._simulation.context.setPositions(self._pdb.positions)
        init_state = self._simulation.context.getState(getEnergy=True)
        
        # Output
        print('Start to minimize system energy', file=self._log_file)
        print(
            'Initial potential energy:\t%.2f\t kj/mol' 
            %(init_state.getPotentialEnergy()/unit.kilojoule_per_mole), file=self._log_file
        )

    def _execute(self):
        self._simulation.minimizeEnergy(tolerance=1e-4 * unit.kilocalorie_per_mole, maxIterations=1000)

    def _teardown(self):
        final_state = self._simulation.context.getState(
            getPositions=True, getParameters=True, 
            enforcePeriodicBox=True, getEnergy=True
        )
        
        # Output
        print(
            'Final potential energy:\t\t%.2f\t kj/mol' 
            %(final_state.getPotentialEnergy()/unit.kilojoule_per_mole), file=self._log_file
        )

        # Write restart file
        file_restart = open(self._out_pdb_restart_file_path, 'w')
        app.PDBFile.writeFile(
            self._simulation.topology, final_state.getPositions(), file_restart
        )