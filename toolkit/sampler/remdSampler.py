#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: remdSampler.py
created time : 2021/06/02
last edit time : 2021/06/02
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import simtk.openmm.openmm as openmm
import simtk.openmm.app as app

class REMDSampler:
    def __init__(self, t_max, t_min, num_replicas: int, switch_freq=1) -> None:
        """
        Parameters
        ----------
        t_max : int or unit
            Maximum temperature of replicas
        t_min : int or unit
            Minimum temperature of replicas
        num_replicas : int
            The number of   
        switch_freq : int, optional
            Switch freqency of REMD in unit of simulation step, by default 1
        """        
        # Read input
        self._t_max = t_max
        self._t_min = t_min
        self._num_replicas = num_replicas
        self._switch_freq = switch_freq
        
        # Get hardware information

        # Calculate attributes
        

        # Set default attributes
        self._topology = None
        self._system = None
        self._integrator = None

    def addTopology(self, topology: app.Topology):
        self._topology = topology

    def addSystem(self, system: openmm.System):
        self._system = system

    def addIntegrator(self, integrator: openmm.Integrator):
        self._integrator = integrator

    def step(self, num_step: int):
        _ = self._checkAttribute()

    def _checkAttribute(self):
        if self._topology == None:
            raise AttributeError('Topology has not been specified')
        elif self._system == None:
            raise AttributeError('System has not been specified')
        elif self._integrator == None:
            raise AttributeError('Integrator has not been specified')
        return 0