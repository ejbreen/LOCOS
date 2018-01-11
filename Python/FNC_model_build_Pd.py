# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 17:08:39 2017

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
    m.Params.OutputFlag = 0
    
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
    m.addConstrs(matches <= assign.T.sum()[t] for t in Treat)
    FNC.printMessage("Constraint 1 done")
    c1_t = FNC.timerStop(c1_t, 5)
    
    c2_t = FNC.timerStart()
    m.addConstrs(1 >= assign.sum()[c] for c in Ctrl)
    FNC.printMessage("Constraint 2 done")
    c2_t = FNC.timerStop(c2_t, 5)

    c3_t = FNC.timerStart()
    m.addConstrs(z[i] >= quicksum(((assign.T[t]*C_pop[i])/(matches*T_n)).sum()
            for t in Treat) + mean_T_pop[i] for i in Covar)
    FNC.printMessage("Constraint 3 done")
    c3_t = FNC.timerStop(c3_t, 5)

    c4_t = FNC.timerStart()
    m.addConstrs(z[i] >= -quicksum(((assign.T[t]*C_pop[i])/(matches*T_n)).sum()
            for t in Treat) - mean_T_pop[i] for i in Covar)    
    FNC.printMessage("Constraint 4 done")
    c4_t = FNC.timerStop(c4_t, 5)
    
    m.update()
    
    setup_time = FNC.timerStop(setup_time, 4)
    
    Timings = [setup_time, c1_t, c2_t, c3_t, c4_t]
    Timings = pd.Series(Timings, index = ["Setup Time","C1","C2","C3","C4"])
    
    return m, assign, Timings





















