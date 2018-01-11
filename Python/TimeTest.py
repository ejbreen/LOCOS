# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 18:08:20 2017

@author: Evan
"""

import numpy as np
import pandas as pd

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC

C_pop_full, T_pop_full = FNC.Import_DataSets(1)
# model options are "Pd" - build using pandas dataframes
#                   "TD" - build using gurobi touple dicts
model = "Pd"

D_sets   = [2, 3]
T_n_sets =          [   4,   8,  16,  32, 64, 128, 256]
rounds  = pd.Series([1000, 500, 500, 100, 50,  50,  50], index = T_n_sets)
Matches = pd.Series([   5,   5,   5,   5,  3,   2,   2], index = T_n_sets)
weights = np.ones(len(T_pop_full.columns))

Timing_Cols = ['model','data set', 'T_n', 'matches', 'Setup Time',
               'c1_t','c2_t','c3_t','c4_t', 'reps', 'Solve Time']
Timing_Data_Pd = pd.DataFrame(0, index = list(range(Matches.sum()*len(D_sets))),
                              columns = Timing_Cols, dtype = np.float)
dataSet = D_sets[0]
writer  = FNC.build_writer("timings_%s_%i-%i_treatments"%
                           ("both_models", T_n_sets[0],T_n_sets[len(T_n_sets)-1]))

# timing data iterator
td = 0
for dataSet in D_sets:
    C_pop_full, T_pop_full = FNC.Import_DataSets(dataSet) 
    for T_n in T_n_sets:
        C_pop, T_pop = FNC.Shrink_pop(C_pop_full, T_pop_full, T_n)
        for matches in list(range(1, Matches[T_n]+1)):
            Timing_Data_Pd.loc[td][0:4] = [0, dataSet, T_n, matches]
            
            #set up the model and return the timings of the different element's
            m, assign, setup_Time = FNC.Build_Model(model, C_pop, T_pop,
                                                    matches, weights)
            Timing_Data_Pd.loc[td][4:9] = setup_Time
            
            FNC.printMessage("Start Presolve")
            m.presolve()
            
            runTimes = pd.Series(range(rounds[T_n]), dtype = np.float)
            for rnd in list(range(rounds[T_n])):
                m.optimize()
                runTimes[rnd] = m.runtime
            Timing_Data_Pd.loc[td][9:11] = [rounds[T_n], round(runTimes.mean(),5)]
            
            td = td + 1
        print "Data set %i, treatment set size %i complete"%(dataSet, T_n)
FNC.write_set(Timing_Data_Pd, 0, writer, model)

model = "TD"
print ""
print ""
print "starting TD testing"
print ""
print ""
Timing_Data_TD = pd.DataFrame(0, index = list(range(Matches.sum()*len(D_sets))),
                              columns = Timing_Cols, dtype = np.float)
td = 0
for dataSet in D_sets:
    C_pop_full, T_pop_full = FNC.Import_DataSets(dataSet)
    for T_n in T_n_sets:
        C_pop, T_pop = FNC.Shrink_pop(C_pop_full, T_pop_full, T_n)
        for matches in list(range(1, Matches[T_n]+1)):
            Timing_Data_TD.loc[td][0:4] = [1, dataSet, T_n, matches]
            
            #set up the model and return the timings of the different element's
            m, assign, setup_Time = FNC.Build_Model(model, C_pop, T_pop,
                                                    matches, weights)
            Timing_Data_TD.loc[td][4:9] = setup_Time
            
            FNC.printMessage("Start Presolve")
            m.presolve()
            
            runTimes = pd.Series(range(rounds[T_n]), dtype = np.float)
            for rnd in list(range(rounds[T_n])):
                m.optimize()
                runTimes[rnd] = m.runtime
            Timing_Data_TD.loc[td][9:11] = [rounds[T_n], round(runTimes.mean(),5)]
            
            td = td + 1
        print "Data set %i, treatment set size %i complete"%(dataSet, T_n)
FNC.write_set(Timing_Data_TD, 0, writer, model)

Timing_Data = Timing_Data_Pd.append(Timing_Data_TD)
FNC.write_set(Timing_Data, 0, writer, "combined")
FNC.write_out(writer)
