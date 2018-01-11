# -*- coding: utf-8 -*-
"""
Created on Tue Dec 05 14:48:06 2017

@author: Evan
"""

import pandas as pd
from sklearn import linear_model as linReg
import statsmodels.api as sm

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC

def compare_MLR(C_pop, T_pop, C_matched):
    
    FNC.printMessage("starting multiple linear regression")
    
    MLR_raw, MLR_sub = linReg.LinearRegression(), linReg.LinearRegression()
    
    Raw = T_pop.append(C_pop)
    Raw_feat = Raw.drop('AGE', axis=1)
    Raw_resp = Raw['AGE']
    
    Sub = T_pop.append(C_matched)
    Sub_feat = Sub.drop('AGE', axis=1)
    Sub_resp = Sub['AGE']
    
    MLR_raw.fit(Raw_feat, Raw_resp)
    MLR_sub.fit(Sub_feat, Sub_resp)
    
    temp = pd.Series(MLR_raw.coef_, index = Raw_feat.columns)
    temp = temp.append(pd.Series({"const":MLR_raw.intercept_}))
    
    temp2 = pd.Series(MLR_sub.coef_, index = Sub_feat.columns)
    temp2 = temp2.append(pd.Series({"const":MLR_sub.intercept_}))
    
    Models = pd.DataFrame({"raw data MRL":temp, "subset MRL":temp2})
    
    Scores = pd.Series({"raw data MRL":MLR_raw.score(Raw_feat, Raw_resp),
                        "subset MRL":MLR_sub.score(Sub_feat, Sub_resp)})
    Scores.name = "R^2 Score"
    
    Models = Models.append(Scores)
        
    return Models, MLR_raw, MLR_sub


def timing_data_regression(Timings):
    Timings = Timings.assign(sqrt_T_n = lambda x: x.T_n**.5)
    Timings = Timings.assign(T_n_sqrd = lambda x: x.T_n**2)
    Timings = Timings.assign(T_n_cubd = lambda x: x.T_n**3)
    
    Timings = Timings.rename(index = str, columns = {'Setup Time':'setup_time'})
    Timings = Timings.rename(index = str, columns = {'Solve Time':'solve_time'})
    
    parameters = 'C(model) + T_n + matches + sqrt_T_n + T_n_sqrd + T_n_cubd'
    
    m_setup = sm.OLS.from_formula('setup_time ~ %s'%(parameters), Timings)
    m_c1 = sm.OLS.from_formula('c1_t ~ %s'%(parameters), Timings)
    m_c2 = sm.OLS.from_formula('c2_t ~ %s'%(parameters), Timings)
    m_c3 = sm.OLS.from_formula('c3_t ~ %s'%(parameters), Timings)
    m_c4 = sm.OLS.from_formula('c4_t ~ %s'%(parameters), Timings)
    m_solve = sm.OLS.from_formula('solve_time ~ %s'%(parameters), Timings)
    
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
