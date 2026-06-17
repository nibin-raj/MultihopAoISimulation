#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Experiment compare between the policies as a function of the K-hop constraint
import LineNetwork_KHopConstraint_AD as AD
import LineNetwork_KHopConstraint_Index_factor as Index

import os

try:
    os.mkdir('Results_EXT_K')
except FileExistsError:
    pass

Ks = [1,2,3,4]#,5,6,7,8,10,15]
ps = 0.8
M = 20
iterations = 2500

OUTPUT_FILE = "Results_EXT_K/Line_ps08_M20_kvary.csv"

outfile = open(OUTPUT_FILE, "w")
for K in Ks:
    tad = AD.simulate_agedifference(iterations, K, M, ps)
    tind = Index.simulate_indexheuristic(iterations, K, M, ps)
    print(K, tad, tind)
    outfile.write("%u,%f, %f\n" % (K, tad, tind))
outfile.close()
    
