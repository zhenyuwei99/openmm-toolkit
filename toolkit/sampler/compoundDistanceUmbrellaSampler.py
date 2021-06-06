#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: compoundDistanceUmbrellaSampler.py
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

class CompoundDistanceUmbrellaSampler(UmbrellaSampler):
    def __init__(
        self, group1_ids, group2_ids, sample_freq,
        distance_min, distance_max, num_replicas: int, 
        harmonic_constant, time_sim, time_step, out_dir, 
        temp_target=300, cut_off=12, pdb_file='', out_freq=10000, 
        out_prefix='distance_umbrella_sampler', platform='CUDA'
    ) -> None:
        super().__init__(
            sample_freq, num_replicas, time_sim, time_step, out_dir, 
            temp_target=temp_target, cut_off=cut_off, pdb_file=pdb_file, 
            out_freq=out_freq, out_prefix=out_prefix, platform=platform
        )

        # Read input
        self._group1_ids = group1_ids
        self._group2_ids = group2_ids
        self._distance_min = check_quantity(distance_min, unit.angstrom)
        self._distance_max = check_quantity(distance_max, unit.angstrom)
        self._harmonic_constant = check_quantity(harmonic_constant, unit.kilojoule_per_mole / unit.angstrom**2)

    def createOriginList(self):
        self._origin_list = np.linspace(self._distance_min / unit.angstrom, self._distance_max / unit.angstrom, self._num_replicas) * unit.angstrom

    def createBiasPotential(self):
        self._bias_potential = openmm.CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-r0)^2')
        self._bias_potential.addGlobalParameter('k', self._harmonic_constant)
        self._bias_potential.addPerBondParameter('r0')
        self._bias_potential.addGroup(self._group1_ids, [])
        self._bias_potential.addGroup(self._group2_ids, [])
        self._bias_potential.setUsesPeriodicBoundaryConditions(True)
        self._bias_potential.addBond([0, 1], [0])
        self._bias_potential.setForceGroup(31)

    def setBiasPotential(self, origin):
        self._bias_potential.setBondParameters(0, [0, 1], [origin])
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

        group1_info = [[position[i] / unit.angstrom, self._system.getParticleMass(i) / unit.dalton] for i in self._group1_ids]
        group2_info = [[position[i] / unit.angstrom, self._system.getParticleMass(i) / unit.dalton] for i in self._group2_ids]

        coord1 = np.zeros(3)
        mass1 = 0
        for info in group1_info:
            coord1 += info[0] * info[1]
            mass1 += info[1]
        coord1 = coord1 / mass1

        coord2 = np.zeros(3)
        mass2 = 0
        for info in group2_info:
            coord2 += info[0] * info[1]
            mass2 += info[1]
        coord2 = coord2 / mass2

        return self._radius(coord1, coord2, pbc)
