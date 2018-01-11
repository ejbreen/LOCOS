# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 17:07:50 2017

@author: Evan
"""
import pandas as pd
from gurobipy import *
from gurobipy import Model
from gurobipy import GRB
from gurobipy import quicksum

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC

#@profile
def Build(C_pop, T_pop, matches, weights):
    
    #start the setup timer
    t_setup = FNC.timerStart()
    
    #set weights for covariates to their min
    weight_base, mean_T_pop,  dist = FNC.Pop_Calcuations(C_pop, T_pop) 
    weights = weights*weight_base
    
    #start creating model elements
    Ctrl  = list(range(len(C_pop)))
    C_pop.index = Ctrl
    Treat = list(range(len(T_pop)))
    T_pop.index = Treat
    T_n = len(T_pop)
    
    Covar = list(T_pop)
    
    #build touple dict structures
    distance = { (t,c) : dist[m][n] for n,c in enumerate(Ctrl)
                                    for m,t in enumerate(Treat)}
    
    Ctrl_pop = { (c,i) : C_pop.loc[m].values[n] for n,i in enumerate(Covar)
                                                for m,c in enumerate(Ctrl)}
    
    T_avg = { i : mean_T_pop[n] for n,i in enumerate (Covar)}
    weight= { i : weights[n] for n,i in enumerate (Covar)}
    
    
    #define the model
    m = Model('match')
    m.Params.OutputFlag = 0
    
    #create variables
    FNC.printMessage("Creating Gurobi Variables")
    assign = m.addVars(distance.keys(), vtype = GRB.BINARY, name = "assign")
    z      = m.addVars(T_avg.keys(), vtype = GRB.CONTINUOUS, name = "z")
    
    
    m.update()
    
    #objective fuction
    FNC.printMessage("Creating Gurobi Objective Function")
    m.setObjective(quicksum(dist[t,c]*assign[t,c] for [t,c] in distance) + 
                   quicksum(weight[i]*z[i] for i in Covar),
                   sense = GRB.MINIMIZE)
    
    #define the constraints
    FNC.printMessage("Creating Gurobi Constraints")
    t_c1 = FNC.timerStart()
    m.addConstrs(matches <= quicksum(assign[t,c] for c in Ctrl) for t in Treat)
    FNC.printMessage("Constraint 1 done")
    t_c1 = FNC.timerStop(t_c1, 5)
    
    t_c2 = FNC.timerStart()
    m.addConstrs(1 >= quicksum(assign[t,c] for t in Treat) for c in Ctrl)
    FNC.printMessage("Constraint 2 done")
    t_c2 = FNC.timerStop(t_c2, 5)

    t_c3 = FNC.timerStart()
    m.addConstrs(z[i] >= quicksum((Ctrl_pop[c,i]*assign[t,c])/(matches*T_n) 
                for t,c in distance) + T_avg[i] for i in Covar)    
    FNC.printMessage("Constraint 3 done")
    t_c3 = FNC.timerStop(t_c3, 5)

    t_c4 = FNC.timerStart()
    m.addConstrs(z[i] >= -quicksum((Ctrl_pop[c,i]*assign[t,c])/(matches*T_n) 
                for t,c in distance) - T_avg[i] for i in Covar)    
    FNC.printMessage("Constraint 4 done")
    t_c4 = FNC.timerStop(t_c4, 5)
    
    m.update()
    
    t_setup = FNC.timerStop(t_setup, 4)
    
    Timings = [t_setup, t_c1, t_c2, t_c3, t_c4]
    Timings = pd.Series(Timings, index = ["Setup Time","C1","C2","C3","C4"])
    
    return m, assign, Timings
























