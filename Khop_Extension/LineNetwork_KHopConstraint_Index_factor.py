#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
from scipy.io import loadmat
from scipy import linalg as la
from itertools import combinations
import logging
import itertools
from scipy.special import comb


def update_age(src, scheduled_links, nodes_links, nodes_ages, ps):
    next_ages = [0] * src
    for l in nodes_links[src]:
        if l in scheduled_links:
            if np.random.rand() <= ps:
                if l < src:
                    next_ages[l - 1] = nodes_ages[src][l] + 1
                else:
                    next_ages[l - 1] = 1
            else:
                next_ages[l - 1] = nodes_ages[src][l - 1] + 1
        else:
            next_ages[l - 1] = nodes_ages[src][l - 1] + 1
    return next_ages


class Source:
    def __init__(self, hs, ps_val, max_tau=100000, max_aoi = 5000):
        self.hs = hs
        self.aoi = 1e-20  # Initial AoI
        self.age_track = []
        self.E_T_cache = {}
        self.E_T2_cache = {}

        self.h = self.hs
        self.s = 0
        self.h_track = []
        self.s_track = []
        self.max_tau=max_tau
        self.E_T = self.compute_expected_times_cached(ps_val, self.hs)[self.hs - 1]
        self.E_T2 = self.compute_second_moment_times_cached(ps_val, self.hs)[self.hs - 1]
        self.P_Tminus_cache = self.precompute_P_Tminus(max_tau, ps_val)

        self.bufferpacket_exists = True
        self.ps_val = ps_val
        self.max_aoi = max_aoi
        self.index_cache = self.precompute_index(self.max_aoi)
        
        self.cumulative_age = 0

    def step(self, scheduled):
        self.cumulative_age += self.aoi
        # self.age_track.append(self.aoi)
        if scheduled:
            if np.random.rand() < self.ps_val:
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

    def compute_link_index(self):
        ind = self.compute_index(self.ps_val)
        link_index = {}
        for i in range(1, self.hs + 1):
            if i == self.h:
                if not self.h == self.hs:
                    link_index[i] = ind
                else:
                    if self.bufferpacket_exists:
                        link_index[i] = ind
                    else:
                        link_index[i] = 0
            else:
                link_index[i] = 0
        return link_index
    
    def precompute_index(self, max_aoi):
        index_cache = []
        for a in range(max_aoi + 1):
            ind = self.compute_index_aoi(a, self.ps_val)
            index_cache.append(ind)
        return index_cache
            
    def compute_index_aoi(self, tau, ps_val):
        D_tau = self.compute_D_tau(tau, ps_val)
        A_tau = self.compute_A_tau(tau, ps_val)
        D_tau_plus1 = self.compute_D_tau(tau + 1, ps_val)
        A_tau_plus1 = self.compute_A_tau(tau + 1, ps_val)

        epsilon1 = 1e-6
        if abs(D_tau - D_tau_plus1) < epsilon1:
            return 0

        ind = (A_tau_plus1 - A_tau) / (D_tau - D_tau_plus1)
        return max(ind[0], 0)

        
    def compute_index(self, ps_val):
        tau = self.aoi
        if tau <= self.max_aoi:
            return self.index_cache[int(tau)]
        else:
            D_tau = self.compute_D_tau(tau, ps_val)
            A_tau = self.compute_A_tau(tau, ps_val)
            D_tau_plus1 = self.compute_D_tau(tau + 1, ps_val)
            A_tau_plus1 = self.compute_A_tau(tau + 1, ps_val)
    
            epsilon1 = 1e-6
            if abs(D_tau - D_tau_plus1) < epsilon1:
                return 0
    
            ind = (A_tau_plus1 - A_tau) / (D_tau - D_tau_plus1)
            return max(ind[0], 0)

class Source_ExternalPacket:
    def __init__(self, hs, ps_val, max_tau=100000, max_aoi = 5000):
        self.hs = hs
        self.aoi = 1e-20  # Initial AoI
        self.age_track = []
        self.E_T_cache = {}
        self.E_T2_cache = {}

        self.h = self.hs
        self.s = 0
        self.h_track = []
        self.s_track = []
        self.max_tau=max_tau
        self.E_T = self.compute_expected_times_cached(ps_val, self.hs)[self.hs - 1]
        self.E_T2 = self.compute_second_moment_times_cached(ps_val, self.hs)[self.hs - 1]
        self.P_Tminus_cache = self.precompute_P_Tminus(max_tau, ps_val)

        self.bufferpacket_exists = False
        self.bufferpacket_age = 0

        self.ps_val = ps_val
        self.max_aoi = max_aoi
        self.index_cache = self.precompute_index(self.max_aoi)
        self.cumulative_age = 0


    def add_packet(self, age):
        self.bufferpacket_exists = True
        self.bufferpacket_age = age

    def step(self, scheduled):
        self.cumulative_age += self.aoi
        # self.age_track.append(self.aoi)
        self.aoi += 1
        if scheduled:
            if np.random.rand() < self.ps_val:
                if self.h == self.hs:
                    if self.bufferpacket_exists:
                        self.h -= 1
                        self.s = self.bufferpacket_age + 1
                        self.bufferpacket_exists = False
                        self.bufferpacket_age = 0
                elif self.h < self.hs and self.h > 1:
                    self.h -= 1
                    self.s += 1
                elif self.h == 1:
                    self.h = self.hs
                    self.aoi = self.s + 1
                    self.s = 0
            else:
                if self.h != self.hs:
                    self.s += 1
        else:
            if self.h != self.hs:
                self.s += 1

            if self.bufferpacket_exists:
                self.bufferpacket_age += 1


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

    def compute_link_index(self):
        ind = self.compute_index(self.ps_val)
        link_index = {}
        for i in range(1, self.hs + 1):
            if i == self.h:
                if not self.h == self.hs:
                    link_index[i] = ind
                else:
                    if self.bufferpacket_exists:
                        link_index[i] = ind
                    else:
                        link_index[i] = 0
            else:
                link_index[i] = 0
        return link_index
    
    def precompute_index(self, max_aoi):
        index_cache = []
        for a in range(max_aoi + 1):
            ind = self.compute_index_aoi(a, self.ps_val)
            index_cache.append(ind)
        return index_cache

    def compute_index_aoi(self, tau, ps_val):
        D_tau = self.compute_D_tau(tau, ps_val)
        A_tau = self.compute_A_tau(tau, ps_val)
        D_tau_plus1 = self.compute_D_tau(tau + 1, ps_val)
        A_tau_plus1 = self.compute_A_tau(tau + 1, ps_val)
    
        epsilon1 = 1e-6
        if abs(D_tau - D_tau_plus1) < epsilon1:
            return 0
    
        ind = (A_tau_plus1 - A_tau) / (D_tau - D_tau_plus1)
        return max(ind[0], 0)
    
    def compute_index(self, ps_val):
        if self.bufferpacket_exists:
            tau = max(self.aoi - self.bufferpacket_age, 0)
        else:
            tau = self.aoi
            
        if tau <= self.max_aoi:
            return self.index_cache[int(tau)]
        else:
            D_tau = self.compute_D_tau(tau, ps_val)
            A_tau = self.compute_A_tau(tau, ps_val)
            D_tau_plus1 = self.compute_D_tau(tau + 1, ps_val)
            A_tau_plus1 = self.compute_A_tau(tau + 1, ps_val)
    
            epsilon1 = 1e-6
            if abs(D_tau - D_tau_plus1) < epsilon1:
                return 0
    
            ind = (A_tau_plus1 - A_tau) / (D_tau - D_tau_plus1)
            return max(ind[0], 0)

def simulate_indexheuristic(iterations, K, M, ps):
    # Initialize node parameters
    ffactor = 1/(max(M - K, 0) + 1)
    activation_sets = []
    for i in range(1, (K + 1) + 1):
        activation_set = range(i, M + 1, K + 1)
        if len(activation_set) > 0:
            activation_sets.append(activation_set)
    # for a in activation_sets:
    #     print([j for j in a])
    
    nodes_links = {}
    nodes_ages = {}
    for i in range(1, M + 1):
        nodes_links[i] = range(1, i + 1) # link (0,1) is 1, (1,2) is 2 etc.
        nodes_ages[i] = [i - x for x in range(i)] # here note that 0th index is the AoI_0 and AoI_src is 0
    
    sources_forindex = {}
    for s in range(1, M + 1):
        if s <= K + 1:
            hs = s
            src = Source(hs, ps)
        else:
            hs = K + 1
            src = Source_ExternalPacket(hs, ps)
        sources_forindex[s] = src
        
    for t in range(iterations):
        ad_link_metric = {}
        ad_link_source = {}
    
        comint_linkmetrics = {}
        for s in range(1, M + 1):
            metric = sources_forindex[s].compute_link_index()
            comint_linkmetrics[s] = metric
    
        
        for l in range(1, K + 2): # this is the initial set of links where complete interference is assumed
            max_ad = -np.inf
            sl = 1
            for s in range(l, M + 1):
                metric = comint_linkmetrics[s][l]
                if metric > max_ad:
                    max_ad = metric
                    sl = s
            ad_link_metric[l] = max_ad
            ad_link_source[l] = sl
            
        for l in range(K + 2, M + 1): # this is the rest of the links where index is computed using AgeDifference
            max_ad = -np.inf
            sl = M
            for s in range(l, M + 1):
                if l == s:
                    agediff_s = nodes_ages[s][l - 1]
                else:
                    agediff_s = np.max([nodes_ages[s][l - 1] - nodes_ages[s][l], 0])
                if agediff_s > max_ad:
                    max_ad = agediff_s
                    sl = s
    
            ad_link_metric[l] = max_ad
            ad_link_source[l] = sl
    
        max_admetric = -np.inf
        actset_index = 0
        for ai, a in enumerate(activation_sets):
            actset_metric = 0
            for l in a:
                f = 1
                if l <= K + 1:
                    f = ffactor
                actset_metric += f * ad_link_metric[l]
                
            if actset_metric >= max_admetric:
                max_admetric = actset_metric
                actset_index = ai
    
        # print("Activation set, source")
        # if t > 4500:
            # print("-" * 80)
            # print("Iteration", t)
            # print("AoI")
            # print([sources_forindex[s].aoi for s in range(1, M + 1)])
            # print("Node AoI")
            # print(nodes_ages)
            # print("Link metrics")
            # print(comint_linkmetrics)
            # print("AD link metric")
            # print(ad_link_metric)
            # print("AD selected source")
            # print(ad_link_source)
            # print(t, [(a, ad_link_source[a]) for a in activation_sets[actset_index]])
    
        # We have obtained the activation set index here
        scheduled_sources_links = {}
        for s in range(1, M + 1):
            scheduled_sources_links[s] = []
    
        for l in activation_sets[actset_index]:
            if l <= K + 1:
                for s in range(1, M + 1):
                    if s == ad_link_source[l]:
                        sources_forindex[s].step(True)
                    else:
                        sources_forindex[s].step(False)
            else:
                scheduled_sources_links[ad_link_source[l]].append(l)
    
        for s in range(1, M + 1):
            t = update_age(s, scheduled_sources_links[s], nodes_links, nodes_ages, ps)
            nodes_ages[s] = t
    
        if K + 2 in activation_sets[actset_index]:
            src = ad_link_source[K + 2]
            sources_forindex[src].add_packet(nodes_ages[src][K + 1])
    
    # scheduler_aoi = [np.mean(sources_forindex[source].age_track) for source in sources_forindex]
    scheduler_aoi = [sources_forindex[source].cumulative_age / iterations for source in sources_forindex]      
    return (np.mean(scheduler_aoi))

if __name__ == "__main__":
    K = 5
    M = 50
    ps = 0.8
    iterations = 20000
    
    t = simulate_indexheuristic(iterations, K, M, ps)
    print(t)

