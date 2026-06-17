#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os

try:
    os.mkdir('Results_EXT_M')
except FileExistsError:
    pass
# Experiment compare between the policies as a function of the K-hop constraint

import LineNetwork_KHopConstraint_AD as AD
import LineNetwork_KHopConstraint_Index_factor as Index

Ms = [5, 10]#, 15, 20, 25, 30, 50, 100]
ps = 0.5
K = 5
iterations = 2500

OUTPUT_FILE = "Results_EXT_M/Line_ps05_k5_Mvary.csv"

outfile = open(OUTPUT_FILE, "w")
for M in Ms:
    tad = AD.simulate_agedifference(iterations, K, M, ps)
    tind = Index.simulate_indexheuristic(iterations, K, M, ps)
    print(M, tad, tind)
    outfile.write("%u,%f, %f\n" % (M, tad, tind))
outfile.close()
    
    
    
