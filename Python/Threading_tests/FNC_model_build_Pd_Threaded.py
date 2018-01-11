# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 17:08:39 2017

@author: Evan
"""
import pandas as pd
import numpy as np
from gurobipy import *
from gurobipy import Model, GRB, quicksum

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC

from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

def constr_1(m, assign, matches, Treat):
    m.addConstrs(matches <= assign.T.sum()[t] for t in Treat)
def constr_2(m, assign, Ctrl):
    m.addConstrs(1 >= assign.sum()[c] for c in Ctrl)
def constr_3(m, assign, z, C_pop, matches, T_n, mean_T_pop, Treat, Covar):
    m.addConstrs(z[i] >= quicksum(((assign.T[t]*C_pop[i])/(matches*T_n)).sum()
            for t in Treat) + mean_T_pop[i] for i in Covar)
def constr_4(m, assign, z, C_pop, matches, T_n, mean_T_pop, Treat, Covar):
    m.addConstrs(z[i] >= -quicksum(((assign.T[t]*C_pop[i])/(matches*T_n)).sum()
            for t in Treat) - mean_T_pop[i] for i in Covar) 
    
def threaded_constr_2(m, assign, Ctrl, threads):
    
    def f(sublist):
        FNC.printMessage("a thread has started at C_pop index %i  ||"%(sublist[0]))
        constr_2(m, assign, sublist)
        FNC.printMessage("the thread started at C_pop index %i has finished ||"%(sublist[0]))
    
    length = len(Ctrl)
    seg = int(round(length/threads, 0))
    sublists = []
    for i in reversed(list(range(threads))):
        temp = list(range(int(round(length-((i+1)*seg),0)),
                          int(round(length-(i*seg),0))))
        sublists.append(temp)
    
    
    pool = ThreadPool(8)
    pool.map(f, sublists)
    pool.close()
    pool.join()
    pool.terminate()
    


C_pop_full, T_pop_full = FNC.Import_DataSets(1)
#########################################################
C_pop, T_pop = FNC.Shrink_pop(C_pop_full, T_pop_full, 50)
#########################################################
weights = np.ones(len(T_pop_full.columns))




def Build(C_pop, T_pop, matches, weights):
    
    #start the setup timer
    setup_time = FNC.timerStart()
    
    #set weights for covariates to their min
    weight_base, mean_T_pop,  dist = FNC.Pop_Calcuations(C_pop, T_pop) 
    weights = weights*weight_base
    
    #start creating model elements
    Ctrl  = list(range(len(C_pop)))
    C_pop.index = Ctrl
    Treat = list(range(len(T_pop)))
    T_pop.index = Treat
    T_n = len(T_pop)
    
    dist = pd.DataFrame(dist, index = Treat, columns = Ctrl)
    
    Covar = list(T_pop)
    
    weights = pd.Series(weights, index = Covar)
    
    #define the model
    m = Model('match')
    m.Params.LogFile = ""
    m.Params.LogToConsole = 0
    
    #create variables
    FNC.printMessage("Creating Gurobi Variables")
    assign = [[m.addVar(vtype = GRB.BINARY, name = "%i, %i" % (t, c)) 
                for c in Ctrl] for t in Treat]
    assign = pd.DataFrame(assign, index = Treat, columns = Ctrl)
    z      = m.addVars(Covar, vtype = GRB.CONTINUOUS, name = "z")
    
    m.update()
    
    #objective fuction
    FNC.printMessage("Creating Gurobi Objective Function")
    m.setObjective((dist*assign).sum().sum() +  
            quicksum(weights[i]*z[i] for i in Covar), 
            sense = GRB.MINIMIZE)
    
    #define the constraints
    FNC.printMessage("Creating Gurobi Constraints")
    c1_t = FNC.timerStart()
    constr_1(m, assign, matches, Treat)
    FNC.printMessage("Constraint 1 done")
    c1_t = FNC.timerStop(c1_t, 3)
    
    c2_t = FNC.timerStart()
    #####################################################
    
    threaded_constr_2(m, assign, Ctrl, 30)
    
    #####################################################
    FNC.printMessage("Constraint 2 done")
    c2_t = FNC.timerStop(c2_t, 3)

    c3_t = FNC.timerStart()
    constr_3(m, assign, z, C_pop, matches, T_n, mean_T_pop, Treat, Covar)
    FNC.printMessage("Constraint 3 done")
    c3_t = FNC.timerStop(c3_t, 3)

    c4_t = FNC.timerStart()
    constr_4(m, assign, z, C_pop, matches, T_n, mean_T_pop, Treat, Covar)   
    FNC.printMessage("Constraint 4 done")
    c4_t = FNC.timerStop(c4_t, 3)
    
    m.update()
    
    setup_time = FNC.timerStop(setup_time, 3)
    
    Timings = [setup_time, c1_t, c2_t, c3_t, c4_t]
    Timings = pd.Series(Timings, index = ["Setup Time","C1","C2","C3","C4"])
    
    return m, assign, Timings



m, assign, Timings = Build(C_pop, T_pop, 5, weights)
m.optimize()

Means, Assignment, C_matched = FNC.compare_PD_output(C_pop, T_pop, assign, 5)

print 'the time it took to create C_2 was %f'%(Timings[2])













