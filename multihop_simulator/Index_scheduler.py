import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import os
from itertools import combinations
from scipy.special import comb

from scipy.io import loadmat
from scipy import linalg as la

import logging
import itertools

from .constants import *

class Source:
    E_T_cache = {}
    E_T2_cache = {}

    def __init__(self, hs, ps_val, max_tau=10000):
        self.hs = hs
        self.aoi = 1e-20  # Initial AoI
        self.age_track = []
        self.h = self.hs
        self.s = 0
        self.h_track = []
        self.s_track = []
        self.max_tau = max_tau
        
        # Precompute expected times
        self.E_T = self.compute_expected_times_cached(ps_val, self.hs)[self.hs - 1]
        self.E_T2 = self.compute_second_moment_times_cached(ps_val, self.hs)[self.hs - 1]
        

        # Precompute P_Tminus values
        self.P_Tminus_cache = self.precompute_P_Tminus(max_tau, ps_val)

    def step(self, scheduled, p_s):
        if scheduled:
            if np.random.rand() < p_s:
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

    def transition_matrix(self, ps_val):
        num_states = self.hs + 1
        T = np.zeros((num_states, num_states))
        for i in range(1, num_states):
            T[i, i] = 1 - ps_val  
            T[i, i - 1] = ps_val  
        T[0, 0] = 1  
        return T[1:, 1:]

    def compute_expected_times_cached(self, ps, hs):
        if (ps, hs) not in self.E_T_cache:
            self.E_T_cache[(ps, hs)] = self.compute_expected_times(ps)
#             print('ET', self.E_T_cache)
        return self.E_T_cache[(ps, hs)]

    def compute_second_moment_times_cached(self, ps, hs):
        if (ps, hs) not in self.E_T2_cache:
            self.E_T2_cache[(ps, hs)] = self.compute_second_moment_times(ps)
        return self.E_T2_cache[(ps, hs)]


    def compute_expected_times(self, ps):
        T = self.transition_matrix(ps)
        I_minus_T = np.eye(self.hs) - T
        ones = np.ones((self.hs, 1))
        return np.linalg.solve(I_minus_T, ones)

    def compute_second_moment_times(self, ps):
        T = self.transition_matrix(ps)
        I_minus_T = np.eye(self.hs) - T
        ones = np.ones((self.hs, 1))
        X = np.linalg.solve(I_minus_T, ones)
        term = ones + 2 * T @ X
        return np.linalg.solve(I_minus_T, term)

    def precompute_P_Tminus(self, max_tau, p_s):
        P_Tminus_values = np.zeros(max_tau + 1)
        for n in range(self.hs, max_tau + 1):
            P_Tminus_values[n] = comb(n - 1, n - self.hs) * (p_s ** self.hs) * ((1 - p_s) ** (n - self.hs))
        return P_Tminus_values

    def compute_D_tau(self, tau, ps):
        tau = int(tau)
        if tau > self.hs:
            a_values = np.arange(self.hs, self.max_tau + 1)
            den = np.sum(self.P_Tminus_cache[:tau] * (tau - np.arange(tau))) + np.sum(self.P_Tminus_cache[a_values] * self.E_T)
            return self.E_T / den
        else:
            return 1

    def compute_A_tau(self, tau, ps):
        tau = int(tau)
        if tau > self.hs:
            a_values = np.arange(self.hs, self.max_tau + 1)
            num = np.sum(self.P_Tminus_cache[a_values] * (a_values * self.E_T)) + \
                  (self.E_T2 - self.E_T) / 2 + \
                  np.sum(self.P_Tminus_cache[:tau] * ( np.arange(tau) * (tau - np.arange(tau)) + ((tau - np.arange(tau)) ** 2 +  (tau - np.arange(tau)) * (2 * self.E_T - 1)) / 2))

            den = np.sum(self.P_Tminus_cache[:tau] * (tau - np.arange(tau))) + np.sum(self.P_Tminus_cache[a_values] * self.E_T)
            return num / den
        else:
            return self.E_T + ((self.E_T2 - self.E_T) * 0.5 / self.E_T)

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
