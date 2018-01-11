# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:43:04 2017

@author: evan
"""

import pandas as pd
from gurobipy import *
from gurobipy import Model, GRB, quicksum

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC

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
    
##### Pandas inputs
    dist = pd.DataFrame(dist, index = Treat, columns = Ctrl)
    Covar = list(T_pop)
    weights = pd.Series(weights, index = Covar)
    
##### Touple Dict inputs
    distance = { (t,c) : dist[m][n] for n,c in enumerate(Ctrl)
                                    for m,t in enumerate(Treat)}
    Ctrl_pop = { (c,i) : C_pop.loc[m].values[n] for n,i in enumerate(Covar)
                                                for m,c in enumerate(Ctrl)}
    T_avg = { i : mean_T_pop[n] for n,i in enumerate (Covar)}
    weight= { i : weights[n] for n,i in enumerate (Covar)}
    
    #define the model
    m = Model('match')
    m.Params.OutputFlag = 0

##### model functions using pandas      
    #create variables
    assign = [[m.addVar(vtype = GRB.BINARY, name = "%i, %i" % (t, c)) 
                for c in Ctrl] for t in Treat]
    assign = pd.DataFrame(assign, index = Treat, columns = Ctrl)
    z      = m.addVars(Covar, vtype = GRB.CONTINUOUS, name = "z")
    m.update()
  
    #objective fuction
    m.setObjective((dist*assign).sum().sum() +  
            quicksum(weights[i]*z[i] for i in Covar), 
            sense = GRB.MINIMIZE)
    
    #define the constraints
    ## this one is a little slow
    m.addConstrs(matches <= assign.T.sum()[t] for t in Treat)
    
    ## this one is insainly slow
    m.addConstrs(1 >= assign.sum()[c] for c in Ctrl)
    
    #these runs much faster than with touple dicts
    m.addConstrs(z[i] >= quicksum(((assign.T[t]*C_pop[i])/(matches*T_n)).sum()
            for t in Treat) + mean_T_pop[i] for i in Covar)
    m.addConstrs(z[i] >= -quicksum(((assign.T[t]*C_pop[i])/(matches*T_n)).sum()
            for t in Treat) - mean_T_pop[i] for i in Covar)    
    
##### model fuctions using touple dicts
    #create variables
    assign = m.addVars(distance.keys(), vtype = GRB.BINARY, name = "assign")
    z      = m.addVars(T_avg.keys(), vtype = GRB.CONTINUOUS, name = "z")
    m.update()
    
    #objective fuction
    m.setObjective(quicksum(dist[t,c]*assign[t,c] for [t,c] in distance) + 
                   quicksum(weight[i]*z[i] for i in Covar),
                   sense = GRB.MINIMIZE)
    
    #define the constraints
    m.addConstrs(matches <= quicksum(assign[t,c] for c in Ctrl) for t in Treat)
    
    #this is just miles and miles faster
    m.addConstrs(1 >= quicksum(assign[t,c] for t in Treat) for c in Ctrl)

    # these preform worse than their pandas counterpart
    m.addConstrs(z[i] >= quicksum((Ctrl_pop[c,i]*assign[t,c])/(matches*T_n) 
                for t,c in distance) + T_avg[i] for i in Covar)    
    m.addConstrs(z[i] >= -quicksum((Ctrl_pop[c,i]*assign[t,c])/(matches*T_n) 
                for t,c in distance) - T_avg[i] for i in Covar)    
                
    
    m.update()
    
    setup_time = FNC.timerStop(setup_time, 3)
    
    Timings = [setup_time, c1_t, c2_t, c3_t, c4_t]
    Timings = pd.Series(Timings, index = ["Setup Time","C1","C2","C3","C4"])
    
    return m, assign, Timings