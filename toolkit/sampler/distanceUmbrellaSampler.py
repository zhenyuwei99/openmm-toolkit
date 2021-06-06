#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: distanceUmbrellaSampler.py
created time : 2021/06/06
last edit time : 2021/06/06
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import os, shutil
import numpy as np
import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit as unit
from . import NVTSampler
from ..utils import check_quantity

class DistanceUmbrellaSampler(NVTSampler):
    def __init__(
        self, atom1_id: int, atom2_id: int, sample_freq,
        distance_min, distance_max, num_replicas: int, harmonic_constant, 
        time_sim, time_step, out_dir, 
        temp_target=300, cut_off=12, pdb_file='', out_freq=10000, 
        out_prefix='distance_umbrella_sampler', platform='CUDA'
    ) -> None:
        super().__init__(temp_target, time_sim, time_step, out_dir, cut_off=cut_off, pdb_file=pdb_file, out_freq=out_freq, out_prefix=out_prefix, platform=platform)

        # Read input
        self._atom1_id = atom1_id
        self._atom2_id = atom2_id
        self._sample_freq = sample_freq
        self._distance_min = check_quantity(distance_min, unit.angstrom)
        self._distance_max = check_quantity(distance_max, unit.angstrom)
        self._num_replicas = num_replicas
        self._harmonic_constant = check_quantity(harmonic_constant, unit.kilojoule_per_mole / unit.angstrom**2)

        # Deduce attributes
        self._origin_list = np.linspace(self._distance_min / unit.angstrom, self._distance_max / unit.angstrom, num_replicas)
        self._num_sim_steps_per_replica = round(self._num_sim_steps / self._num_replicas)
        self._num_samples_per_replica = round(self._num_sim_steps_per_replica / self._sample_freq)

        # Check and define path
        self._out_cv_dir = os.path.join(self._out_sample_dir, out_prefix)
        if not os.path.exists(self._out_cv_dir):
            os.mkdir(self._out_cv_dir)

    @staticmethod
    def radius(coord1, coord2, pbc):
        pbc = np.array(pbc / unit.angstrom)
        pbc_inv = np.linalg.inv(pbc)
        coord1 = pbc_inv.dot(np.array(coord1 / unit.angstrom))
        coord2 = pbc_inv.dot(np.array(coord2 / unit.angstrom))
        
        diff = coord1-coord2
        diff = diff - np.round(diff)
        diff = pbc.dot(diff)
        
        return np.sqrt( (diff**2).sum() )

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

        # Bias potential
        self._bias_potential = openmm.CustomBondForce('0.5 * k * (r - r0)^2')
        self._bias_potential.addGlobalParameter('k', self._harmonic_constant)
        self._bias_potential.addPerBondParameter('r0')
        self._bias_potential.setUsesPeriodicBoundaryConditions(True)
        self._bias_potential.addBond(self._atom1_id, self._atom2_id, [0])
        self._bias_potential.setForceGroup(30)

        # System
        self._system = self._force_field.createSystem(
            self._pdb.topology, nonbondedMethod=app.PME,nonbondedCutoff=self._cut_off, 
            constraints=app.HBonds, ewaldErrorTolerance=0.0005
        )
        self._system.addForce(self._bias_potential)

        # Integrator
        self._integrator = openmm.LangevinIntegrator(self._temp_target, 0.001/unit.femtosecond, self._time_step)

        # Simulation
        self._simulation = app.Simulation(self._pdb.topology, self._system, self._integrator, self._platform)     
        self._simulation.context.setPositions(self._pdb.positions)
        self._simulation.context.setVelocitiesToTemperature(self._temp_target)
        self._simulation.reporters.append(self._log_reporter)
        self._simulation.reporters.append(self._pdb_reporter) 

    def _execute(self):
        # Conformation refinement
        self._bias_potential.setBondParameters(0, self._atom1_id, self._atom2_id, [self._origin_list[0]*unit.angstrom])
        self._bias_potential.updateParametersInContext(self._simulation.context)
        self._simulation.step(5000)

        for index, origin in enumerate(self._origin_list):
            # Initialization
            sample_res = []
            self._bias_potential.setBondParameters(0, self._atom1_id, self._atom2_id, [origin*unit.angstrom])
            self._bias_potential.updateParametersInContext(self._simulation.context)
            cv_file = os.path.join(self._out_cv_dir, 'replica_%d.npz' %index)

            # Log output
            print(
                '%d steps in one run centered at %.2f A, giving %d samples' 
                %(self._num_sim_steps_per_replica, origin, self._num_samples_per_replica), file=self._log_file
            )

            # Sample
            for sample in range(self._num_samples_per_replica):
                self._simulation.step(self._sample_freq)
                state = self._simulation.context.getState(getPositions=True)
                position = state.getPositions()
                pbc = state.getPeriodicBoxVectors()
                sample_res.append(self.radius(position[self._atom1_id], position[self._atom2_id], pbc))
                if sample % 200 == 0:
                    np.savez(
                        cv_file, origin=origin, sample_res=sample_res,
                        harmonic_constant=self._harmonic_constant*unit.angstrom**2/unit.kilojoule_per_mole
                    )

            # Save sample result
            np.savez(
                cv_file, origin=origin, sample_res=sample_res,
                harmonic_constant=self._harmonic_constant*unit.angstrom**2/unit.kilojoule_per_mole
            )
            

    def _teardown(self):
        # Write restart file
        pdb_restart_file = open(self._out_pdb_restart_file_path, 'w')

        final_state = self._simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True)
        self._simulation.topology.setPeriodicBoxVectors(final_state.getPeriodicBoxVectors())

        app.PDBFile.writeFile(self._simulation.topology, final_state.getPositions(), pdb_restart_file)

        pdb_restart_file.close()

    @property
    def out_prefix(self):
        return self._out_prefix

    @out_prefix.setter
    def out_prefix(self, out_prefix):
        self._out_prefix = out_prefix
        self._out_log_file_path = os.path.join(self._out_log_dir, out_prefix + '.log')
        self._out_pdb_file_path = os.path.join(self._out_pdb_dir, out_prefix + '.pdb')
        self._out_pdb_restart_file_path = os.path.join(self._out_pdb_dir, out_prefix + '_restart.pdb')
        shutil.rmtree(self._out_cv_dir)
        self._out_cv_dir = os.path.join(self._out_sample_dir, out_prefix)
        os.mkdir(self._out_cv_dir)