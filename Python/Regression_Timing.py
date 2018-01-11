#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 18:41:31 2017

@author: evan
"""

import pandas as pd
import statsmodels.api as sm

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC

FNC.printMessage("starting timing regression")

Timings = pd.read_excel("Data/timings_both_models_4-128_treatments.xlsx", 2)

def timing_data_regression(Timings):
    Timings = Timings.assign(sqrt_T_n = lambda x: x.T_n**.5)
    Timings = Timings.assign(T_n_sqrd = lambda x: x.T_n**2)
    Timings = Timings.assign(T_n_cubd = lambda x: x.T_n**3)
    
    Timings = Timings.rename(index = str, columns = {'Setup Time':'setup_time'})
    Timings = Timings.rename(index = str, columns = {'Solve Time':'solve_time'})
    
    parameters = 'C(model) + T_n + matches + sqrt_T_n + T_n_sqrd + T_n_cubd'
    
    m_setup=sm.OLS.from_formula('setup_time ~ %s'%(parameters), Timings)
    m_c1 = sm.OLS.from_formula('c1_t ~ %s'%(parameters), Timings)
    m_c2 = sm.OLS.from_formula('c2_t ~ %s'%(parameters), Timings)
    m_c3 = sm.OLS.from_formula('c3_t ~ %s'%(parameters), Timings)
    m_c4 = sm.OLS.from_formula('c4_t ~ %s'%(parameters), Timings)
    m_solve=sm.OLS.from_formula('solve_time ~ %s'%(parameters), Timings)
    
    models = [m_setup, m_c1, m_c2, m_c3, m_c4, m_solve]
    
    r = m_setup.fit()
    results = [r, r, r, r, r, r]
    p = r.params
    params = pd.DataFrame({0:p, 1:p})
    pval = r.pvalues
    pvalues = pd.DataFrame({0:pval, 1:pval})
    
    i=0
    for m in models:
        results[i]=m.fit()
        params[i]=m.fit().params
        pvalues[i] = m.fit().pvalues
        i=i+1
    
    
    params = params.rename(index=str, columns={0:'setup', 1:'c1', 2:'c2', 3:'c3', 4:'c4', 5:'solve'})
    pvalues = pvalues.rename(index=str, columns={0:'setup', 1:'c1', 2:'c2', 3:'c3', 4:'c4', 5:'solve'})
    
    return params, pvalues, results, models

params, pvalues, results, models = timing_data_regression(Timings)
params_Pd, pvalues_Pd, y, z = timing_data_regression(Timings.loc[Timings['model']==0])
params_TD, pvalues_TD, y, z = timing_data_regression(Timings.loc[Timings['model']==1])






