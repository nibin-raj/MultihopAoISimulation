import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import os
from itertools import combinations
# from scipy.special import comb

from scipy.io import loadmat
from scipy import linalg as la

import logging
import itertools
from math import comb

from .constants import *

class Source:
    def __init__(self, hs, ps_val, max_tau=10000):
        self.hs = hs
        self.aoi = 0 #1e-20  # Initial AoI
        self.age_track = []
        self.h = self.hs
        self.s = 0
        self.h_track = []
        self.s_track = [] 
        self.max_tau =max_tau
        self.ps_val = ps_val
        
        self.rho_cache = self._precompute_rho(self.ps_val)
        
    def _precompute_rho(self, pl):
        rho_vals = np.zeros(self.max_tau)
        if self.hs == 1:
            rho_vals[1] = 1.0
        else:
            for a in range(self.hs, self.max_tau):
                coeff = comb(a - 2, self.hs - 2)
                rho_vals[a] = coeff * (pl ** (self.hs - 1)) * ((1 - pl) ** (a - self.hs))
        return rho_vals
    
    def step(self, scheduled, p_l):
        if scheduled:
            if np.random.rand() < p_l:
                if self.h > 1:
                    self.h -= 1
                    self.s += 1
                    self.aoi += 1
                elif self.h == 1:
                    self.h = self.hs
                    self.aoi = self.s + 1
                    self.s = 0
            else:
                if self.h != self.hs:
                    self.s += 1
                self.aoi += 1  
        else:
            if self.h != self.hs:
                self.s += 1
            self.aoi += 1  
        self.age_track.append(self.aoi)

    
    def compute_D_tau(self, tau_s, pl):  
        firstterm = self.hs / pl
        sum_term = sum(self.rho_cache[a] * (tau_s - a) for a in range(self.hs, tau_s))
        total_cycle_time = firstterm + sum_term
        return firstterm / total_cycle_time

    def compute_A_tau(self, tau_s, pl):
        numerator_term1 = (self.hs * (3 * self.hs - 1)) / (2 * pl**2)
        numerator_sum = 0.0
        denominator_sum = 0.0

        for a in range(self.hs, tau_s):
            rho_a = self.rho_cache[a]
            delta = tau_s - a
            term = (delta / 2) * (a + tau_s + (2 * self.hs / pl) - 1)
            numerator_sum += rho_a * term
            denominator_sum += rho_a * delta

        numerator = numerator_term1 + numerator_sum
        denominator = (self.hs / pl) + denominator_sum
        return numerator / denominator



    def compute_index(self, ps_val):
        tau = self.aoi
        D_tau = self.compute_D_tau(tau, ps_val)
        A_tau = self.compute_A_tau(tau, ps_val)
        D_tau_plus1 = self.compute_D_tau(tau + 1, ps_val)
        A_tau_plus1 = self.compute_A_tau(tau + 1, ps_val)

        epsilon1 = 1e-6
        if abs(D_tau - D_tau_plus1) < epsilon1:
            return np.Inf
        return (A_tau_plus1 - A_tau) / (D_tau - D_tau_plus1)

