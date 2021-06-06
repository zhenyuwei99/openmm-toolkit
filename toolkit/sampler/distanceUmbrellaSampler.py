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

import numpy as np
import simtk.openmm as openmm
import simtk.unit as unit
from . import UmbrellaSampler
from ..utils import check_quantity

class DistanceUmbrellaSampler(UmbrellaSampler):
    def __init__(
        self, atom1_id: int, atom2_id: int, sample_freq,
        distance_min, distance_max, num_replicas: int, harmonic_constant, 
        time_sim, time_step, out_dir, 
        temp_target=300, cut_off=12, pdb_file='', out_freq=10000, 
        out_prefix='distance_umbrella_sampler', platform='CUDA'
    ) -> None:
        super().__init__(
            sample_freq, num_replicas, time_sim, time_step, out_dir, 
            temp_target=temp_target, cut_off=cut_off, pdb_file=pdb_file, 
            out_freq=out_freq, out_prefix=out_prefix, platform=platform
        )

        # Read input
        self._atom1_id = atom1_id
        self._atom2_id = atom2_id
        self._distance_min = check_quantity(distance_min, unit.angstrom)
        self._distance_max = check_quantity(distance_max, unit.angstrom)
        self._harmonic_constant = check_quantity(harmonic_constant, unit.kilojoule_per_mole / unit.angstrom**2)

    def createOriginList(self):
        self._origin_list = np.linspace(self._distance_min / unit.angstrom, self._distance_max / unit.angstrom, self._num_replicas) * unit.angstrom

    def createBiasPotential(self):
        self._bias_potential = openmm.CustomBondForce('0.5 * k * (r - r0)^2')
        self._bias_potential.addGlobalParameter('k', self._harmonic_constant)
        self._bias_potential.addPerBondParameter('r0')
        self._bias_potential.setUsesPeriodicBoundaryConditions(True)
        self._bias_potential.addBond(self._atom1_id, self._atom2_id, [0])
        self._bias_potential.setForceGroup(31)

    def setBiasPotential(self, origin):
        self._bias_potential.setBondParameters(0, self._atom1_id, self._atom2_id, [origin])
        self._bias_potential.updateParametersInContext(self._simulation.context)

    @staticmethod
    def _radius(coord1, coord2, pbc):
        pbc = np.array(pbc / unit.angstrom)
        pbc_inv = np.linalg.inv(pbc)
        coord1 = pbc_inv.dot(np.array(coord1 / unit.angstrom))
        coord2 = pbc_inv.dot(np.array(coord2 / unit.angstrom))
        
        diff = coord1 - coord2
        diff = diff - np.round(diff)
        diff = pbc.dot(diff)
        
        return np.sqrt( (diff**2).sum())

    def getCV(self):
        state = self._simulation.context.getState(getPositions=True)
        position = state.getPositions()
        pbc = state.getPeriodicBoxVectors()
        return self._radius(position[self._atom1_id], position[self._atom2_id], pbc)