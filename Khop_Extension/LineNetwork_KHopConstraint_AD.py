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

def simulate_agedifference(iterations, K, M, ps):
    # Initialize node parameters
    activation_sets = []
    for i in range(1, (K + 1) + 1):
        activation_sets.append(range(i, M + 1, K + 1))
        
    nodes_links = {}
    nodes_ages = {}
    cumulative_ages = [0] * M
    for i in range(1, M + 1):
        nodes_links[i] = range(1, i + 1) # link (0,1) is 1, (1,2) is 2 etc.
        nodes_ages[i] = [i - x for x in range(i)] # here note that 0th index is the AoI_0 and AoI_src is 0
    
    for t in range(iterations):
        ad_link_metric = {}
        ad_link_source = {}
            
        for l in range(1, M + 1):
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
                actset_metric += ad_link_metric[l]
                
            if actset_metric > max_admetric:
                max_admetric = actset_metric
                actset_index = ai
        # if t > 4000:
        #   print(t, [(a, ad_link_source[a]) for a in activation_sets[actset_index]])

        scheduled_sources_links = {}
        for s in range(1, M + 1):
            scheduled_sources_links[s] = []

        for l in activation_sets[actset_index]:
            scheduled_sources_links[ad_link_source[l]].append(l)
        
        for s in range(1, M + 1):
            cumulative_ages[s - 1] = cumulative_ages[s - 1] + nodes_ages[s][0]
            t = update_age(s, scheduled_sources_links[s], nodes_links, nodes_ages, ps)
            nodes_ages[s] = t
            
    AoI_0 = []
    for s in range(1, M + 1):
        AoI_0.append(cumulative_ages[s - 1])
    return np.mean(AoI_0)/iterations

if __name__ == "__main__":
    b = simulate_agedifference(20000, 4, 10, 0.5)
    print(b)