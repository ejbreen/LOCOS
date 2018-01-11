# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 22:22:03 2017

@author: evan
"""
import pandas as pd
import numpy as np
import scipy.spatial
import xlsxwriter
import platform
import time

import sys
import os
sys.path.append(os.path.abspath("Python/"))
import FNC_model_build_Pd
import FNC_model_build_TD
import FNC_regression

#sets the working directory for the os that we are in
#this way once the directory is known for an os, we can
#let it set itself up automatically
#----This probably needs to be updated for the pixelbook and Flux
def Set_WD():
    def winset():
        os.chdir(r'D:/Programing/TKA_T-C_Matching')
    def linset():
        linux_login = pd.Series({"Flux":"ejbreen", "pixel":"evan"})
        # change this to pixelbook when running there but change back when done
        os.chdir(r'/home/%s/TKA_T-C_Matching'%linux_login["Flux"])
    setdir = {'Windows' : winset,
              'Linux'   : linset}
    setdir[platform.system()]()

#Import the specified dataset (1, 2, or 3) from the data folder into 
#C_pop_full and T_pop_full
def Import_DataSets(dataSet): 
    Set_WD()
    filePath = "Data/"
    fileName_C = "C_pop "
    fileName_T = "T_pop "
    fileExt = " .csv"
    C_pop_full = pd.read_csv("%s%s%i%s"%(filePath, fileName_C, 
                                         dataSet, fileExt), index_col = 0)
    T_pop_full = pd.read_csv("%s%s%i%s"%(filePath, fileName_T, 
                                         dataSet, fileExt), index_col = 0)
    return C_pop_full, T_pop_full

#shrink the full datasets to down to T_n for T_pop and T_n*30 for C_pop
#for the shrink, it takes the top of the list
def Shrink_pop(C_pop_full, T_pop_full, T_n):
    C_pop = C_pop_full.head(T_n*30)
    T_pop = T_pop_full.head(T_n)
    return C_pop, T_pop

#calculate the min weights, mean_T_pop, and bray curtis distance
def Pop_Calculations(C_pop, T_pop):
    weights = np.tile(len(T_pop)*.1,len(T_pop.columns))
    weights = pd.Series(weights, index = list(T_pop))
    mean_T_pop = T_pop.mean()
    dist = scipy.spatial.distance.cdist(T_pop, C_pop, 'braycurtis')
    return weights, mean_T_pop, dist

#start time of some timer
def timerStart():
    timer = time.clock()
    return timer
#end time of some timer
def timerStop(timer, sig):
    timer = time.clock()-timer
    timer = round(timer, sig)
    return timer
#prints a message with time header
def printMessage(message):
    tm = time.localtime()
    timeStr = '%i:%i:%i'%(tm.tm_hour, tm.tm_min, tm.tm_sec)
    print ('%s %s' %(timeStr, message))
    
#export the timing data to an excel file
def build_writer(fileName):
    filePath = "Data/"
    fileName_out = "/%s"%(fileName)
    writer  = pd.ExcelWriter("%s%s%s"%(filePath, fileName_out, '.xlsx'),
                             engine = 'xlsxwriter')
    return writer
def write_set(TimingData, dataSet, writer, modelType):
    TimingData.to_excel(writer,
                        "Data Set %i Timings (%s)"%(dataSet, modelType))
def write_df(Data, message, writer, modelType):
    Data.to_excel(writer, "%s (%s)"%(message, modelType))
def write_out(writer):
    writer.save()
    print "data exported to excel"
    
#a simple writer for CSVs
def csv_write(fileName, Data):
    import csv
    with open('%s.csv'%(fileName), 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(Data)
    print ('%s written as csv'%(fileName))

# models are defined in their own files, this is just meant as a pass 
# through to FNC so FNC can be a one stop shop
def Build_Pd_Model(C_pop, T_pop, matches, weights): 
    return FNC_model_build_Pd.Build(C_pop, T_pop, matches, weights)

def Build_TD_Model(C_pop, T_pop, matches, weights):
    return FNC_model_build_TD.Build(C_pop, T_pop, matches, weights)

#allows you to select the which model you want to create 
def Build_Model(model_base, C_pop, T_pop, matches, weights):
    selectM = {'Pd': Build_Pd_Model,
               'TD': Build_TD_Model}
    m, assign, Timings = selectM[model_base](C_pop, T_pop, matches, weights)
    return m, assign, Timings


# outputs the average value of C_pop, T_pop, and the elements of C_pop 
# selected in by the optimization
def compare_Pd_output(C_pop, T_pop, assign, matches):
    Assignment = pd.DataFrame(0, index = C_pop.index, columns = T_pop.index)
    for i in T_pop.index:
        for j in C_pop.index:
            Assignment[i][j] = assign[j][i].X
    
    C_matched = pd.DataFrame(0, index = range(matches*len(T_pop)), 
                                 columns = C_pop.columns)
    j=0
    for i in C_pop.index:
        if Assignment.loc[i].sum() >=1:
            C_matched.loc[j] = C_pop.loc[i]
            j = j+1
            
    mean_C_matched = C_matched.mean()
    mean_C_pop = C_pop.mean()
    mean_T_pop = T_pop.mean()
    means = pd.concat([mean_T_pop, mean_C_matched, mean_C_pop], axis=1)
    means.columns = ['Treatment', 'Matched Ctrl', 'Control']
    
    return means, Assignment, C_matched

def compare_TD_output(C_pop, T_pop, assign, matches):
    Assignment = pd.DataFrame(0, index = C_pop.index, columns = T_pop.index)
    for i in T_pop.index:
        for j in C_pop.index:
            Assignment[i][j] = assign[i,j].X
    
    C_matched = pd.DataFrame(0, index = range(matches*len(T_pop)), 
                                 columns = C_pop.columns)
    j=0
    for i in C_pop.index:
        if Assignment.loc[i].sum() >=1:
            C_matched.loc[j] = C_pop.loc[i]
            j = j+1
            
    mean_C_matched = C_matched.mean()
    mean_C_pop = C_pop.mean()
    mean_T_pop = T_pop.mean()
    means = pd.concat([mean_T_pop, mean_C_matched, mean_C_pop], axis=1)
    means.columns = ['Treatment', 'Matched Ctrl', 'Control']
    
    return means, Assignment, C_matched

# compare the multiple linear regression models of the raw data 
# and the selected subset
def compare_MLR(C_pop, T_pop, C_matched):
    return FNC_regression.compare_MLR(C_pop, T_pop, C_matched)

def timing_data_regression(Timings):
    return FNC_regression.timing_data_regression(Timings)

#get the files needed in place to run julia within the Julia folder
def setup_julia(dataSet, T_n):
    C, T = Import_DataSets(dataSet)
    T = T.head(T_n)
    C = C.head(T_n*30)
    a, b, dist = Pop_Calculations(C, T)
    T.to_csv("Julia/T_pop_%i.csv"%(dataSet), index = False)
    C.to_csv('Julia/C_pop_%i.csv'%(dataSet), index = False)
    csv_write('Julia/dist_%i'%(dataSet), dist)




















