#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: pressureEquilibrator.py
created time : 2021/06/02
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

class PressureEquilibrator(Equilibrator):
    def __init__(
        self, press_target, temp_target, time_sim, time_step, out_freq, out_dir, 
        cut_off=12, pdb_file='', out_prefix='pressure_equilibrator', platform='CUDA'
    ) -> None:
        super().__init__(out_dir, cut_off, pdb_file, out_prefix, platform)

        # Read input
        self._press_target = check_quantity(press_target, unit.bar)
        self._temp_target = check_quantity(temp_target, unit.kelvin)
        self._time_sim = check_quantity(time_sim, unit.femtosecond)
        self._time_step = check_quantity(time_step, unit.femtosecond)
        self._out_freq = out_freq

        # Deduce attribute
        self._num_sim_steps = round(self._time_sim / self._time_step)

    def _setup(self):
        # Force Field
        self._force_field = app.ForceField(
            'amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml'
        )

        # Log reporter
        self._log_reporter = app.StateDataReporter(
            self._log_file, self._out_freq, step=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, speed=True, density=True,
            totalSteps=self._num_sim_steps, remainingTime=True,separator='\t'
        )

        # PDB reporter
        self._pdb_reporter = app.PDBReporter(
            self._out_pdb_file_path, self._out_freq, enforcePeriodicBox=True
        )

        # System
        self._system = self._force_field.createSystem(
            self._pdb.topology, nonbondedMethod=app.PME,nonbondedCutoff=self._cut_off, 
            constraints=app.HBonds, ewaldErrorTolerance=0.0005
        )
        barostat = openmm.MonteCarloBarostat(self._press_target, self._temp_target, 1000)
        self._system.addForce(barostat)

        # Integrator
        self._integrator = openmm.LangevinIntegrator(self._temp_target, 0.001/unit.femtosecond, self._time_step)

        # Simulation
        self._simulation = app.Simulation(self._pdb.topology, self._system, self._integrator, self._platform)     
        self._simulation.context.setPositions(self._pdb.positions)
        self._simulation.context.setVelocitiesToTemperature(self._temp_target)
        self._simulation.reporters.append(self._log_reporter)
        self._simulation.reporters.append(self._pdb_reporter)   

    def _execute(self):
        self._simulation.step(self._num_sim_steps)

    def _teardown(self):
        # Write restart file
        pdb_restart_file = open(self._out_pdb_restart_file_path, 'w')

        final_state = self._simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True)
        self._simulation.topology.setPeriodicBoxVectors(final_state.getPeriodicBoxVectors())

        app.PDBFile.writeFile(self._simulation.topology, final_state.getPositions(), pdb_restart_file)

        pdb_restart_file.close()