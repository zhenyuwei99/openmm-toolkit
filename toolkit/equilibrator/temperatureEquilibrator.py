#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: temperatureEquilibrator.py
created time : 2021/06/02
last edit time : 2021/06/05
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import sys, datetime
import numpy as np
import simtk.openmm.openmm as openmm
import simtk.openmm.app as app
import simtk.unit as unit
from . import Equilibrator
from ..utils import *

class TemperatureEquilibrator(Equilibrator):
    def __init__(
        self, temp_target, temp_step, time_sim, time_step, output_freq, output_dir,
        cut_off=12, pdb_file='', out_prefix='temperature_equilibrator', platform='CUDA'
    ) -> None:
        super().__init__(output_dir, cut_off, pdb_file, out_prefix, platform)

        # Read input
        self._temp_target = check_quantity(temp_target, unit.kelvin)
        self._temp_step = check_quantity(temp_step, unit.kelvin)  
        self._time_sim = check_quantity(time_sim, unit.femtosecond)
        self._time_step = check_quantity(time_step, unit.femtosecond)
        self._output_freq = output_freq

        # Deduce attribute
        self._num_sim_steps = round(self._time_sim / self._time_step)
        self._num_episodes = round((self._temp_target - 1 * unit.kelvin) / self._temp_step) # Heating from 1K
        self._num_sim_steps_per_episode = round(self._num_sim_steps / self._num_episodes)

    def _setup(self):
        # Output equilibrator info
        print(
            'Equilibrator heats the system from 1 K to %.2f K\n' 
            %(self._temp_target / unit.kelvin)
        )
        print(
            '%d %d-steps (%.2f ps) episodes will be performed.' 
            %(self._num_episodes, self._num_sim_steps_per_episode, self._num_sim_steps_per_episode * self._time_step / unit.picosecond)
        )

        # Force Field
        self._force_field = app.ForceField(
            'amber14/protein.ff14SB.xml', 'amber14/tip3pfb.xml'
        )

        # Log reporter
        self._log_reporter = app.StateDataReporter(
            sys.stdout, self._output_freq, step=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, speed=True, density=True,
            totalSteps=self._num_sim_steps, remainingTime=True,separator='\t'
        )

        # PDB reporter
        self._pdb_reporter = app.PDBReporter(
            self._out_pdb_file_path, self._output_freq, enforcePeriodicBox=True
        )

    def _execute(self):
        for temp in np.linspace(1, self._temp_target / unit.kelvin, self._num_episodes):
            # Output
            print('\nHeating system at %.2f K. \n' %temp) 
            
            # System
            system = self._force_field.createSystem(
                self._pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=self._cut_off, 
                constraints=app.HBonds, ewaldErrorTolerance=0.0005
            )

            # Integrator
            integrator = openmm.LangevinIntegrator(temp*unit.kelvin, 0.001/unit.femtosecond, self._time_step)  

            # Simulation
            simulation = app.Simulation(self._pdb.topology, system, integrator, self._platform)
            simulation.context.setPositions(self._pdb.positions)
            simulation.context.setVelocitiesToTemperature(temp*unit.kelvin)
            simulation.reporters.append(self._log_reporter)
            simulation.reporters.append(self._pdb_reporter)

            simulation.step(self._num_sim_steps_per_episode)

            # Write restart file
            final_state = simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True)
            file_restart = open(self._out_pdb_restart_file_path, 'w')
            app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), file_restart)
            file_restart.close()
            self._pdb = app.PDBFile(self._out_pdb_restart_file_path)

    def _teardown(self):
        pass