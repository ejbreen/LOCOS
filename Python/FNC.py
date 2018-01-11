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

def Set_WD():
    """
    Sets the working directory for the os that we are in
    this way once the directory is known for an os, we can
    let it set itself up automatically.
    
    Inputs:     none
    Outputs:    none
    
    This probably needs to be updated for the pixelbook and Flux
    """
    def winset():
        os.chdir(r'D:/Programing/LOCOS')
    def linset():
        linux_login = pd.Series({"Flux":"ejbreen", "pixel":"evan"})
        # change this to pixelbook when running there but change back when done
        os.chdir(r'/home/%s/LOCOS'%linux_login["Flux"])
    setdir = {'Windows' : winset,
              'Linux'   : linset}
    setdir[platform.system()]()


def Import_DataSets(dataSet): 
    """
    Imports the#Import the specified dataset (1, 2, or 3) CSV from the data 
    folder into C_pop_full and T_pop_full CSV that was created in R.
    
    Inputs:     
        dataSet(int) - The version of the dataset you with to import
    Outputs:    
        C_pop_full - Dataframe imported from Data/C_pop(dataSet).csv
        T_pop_full - Dataframe imported from Data/T_pop(dataSet).csv

    """
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

def Shrink_pop(C_pop_full, T_pop_full, T_n):
    """
    Shrink the full datasets to down to T_n for T_pop and T_n*30 for C_pop
    for the shrink, it takes the top of the list.
    Inputs:     
        C_pop_full  - The full C_pop 
        T_pop_full  - The full T_pop
        T_n(int)    - The number of observations in the output T_pop
    Outputs:    
        C_pop - The shrunken C_pop
        T_pop - The shrunken T_pop
    """
    C_pop = C_pop_full.head(T_n*30)
    T_pop = T_pop_full.head(T_n)
    return C_pop, T_pop

def Pop_Calculations(C_pop, T_pop):
    """
    Preform calculations on C_pop and T_pop to get the inputs requrired
    in the linear programs.
    
    Inputs:     
        C_pop - The population of potential controls
        T_pop - The population of potential treatments
    Outputs:    
        weights - The base weights for the linear program
        mean_T_pop - The mean value of each covariate in the T_pop
        dist    - The Bray Curtis distance between each element of the 
                    treatment population and each element of the potential
                    control population
    """
    weights = np.tile(len(T_pop)*.1,len(T_pop.columns))
    weights = pd.Series(weights, index = list(T_pop))
    mean_T_pop = T_pop.mean()
    dist = scipy.spatial.distance.cdist(T_pop, C_pop, 'braycurtis')
    return weights, mean_T_pop, dist

class Timer:
    """
    A simple timer object
    
    Attributes:
        t       - The time recorded when stop() is called, -1 if 
                    stop not called
        t_start - The clock time recorded at timer creation or after a reset
        sig     - The number of significant figures used
    """
    def __init__(self, t=-1, t_start=time.clock(), sig=5):    
        self.t   = t
        self.t_start = t_start
        self.sig = sig
    
    def reset_t(self):
        self.t_start = time.clock()
    def stop(self):
        self.t = time.clock() - self.t
        self.t = round(self.t, self.sig)
    def sig_set(self, sigfig):
        self.sig = sigfig

##testing the Timer object
t = Timer()
t.sig
t.reset_t()
t.t_start
t.stop()
t.t
#works

#start time of some timer
def timerStart():
    """
    Inputs:     
    Outputs:    
    """
    timer = time.clock()
    return timer
#end time of some timer
def timerStop(timer, sig):
    """
    Inputs:     
    Outputs:    
    """
    timer = time.clock()-timer
    timer = round(timer, sig)
    return timer
#prints a message with time header
def printMessage(message):
    """
    Inputs:     
    Outputs:    
    """
    tm = time.localtime()
    timeStr = '%i:%i:%i'%(tm.tm_hour, tm.tm_min, tm.tm_sec)
    print ('%s %s' %(timeStr, message))
    
#export the timing data to an excel file
def build_writer(fileName):
    """
    Inputs:     
    Outputs:    
    """
    filePath = "Data/"
    fileName_out = "/%s"%(fileName)
    writer  = pd.ExcelWriter("%s%s%s"%(filePath, fileName_out, '.xlsx'),
                             engine = 'xlsxwriter')
    return writer
def write_set(TimingData, dataSet, writer, modelType):
    """
    Inputs:     
    Outputs:    
    """
    TimingData.to_excel(writer,
                        "Data Set %i Timings (%s)"%(dataSet, modelType))
def write_df(Data, message, writer, modelType):
    """
    Inputs:     
    Outputs:    
    """
    Data.to_excel(writer, "%s (%s)"%(message, modelType))
def write_out(writer):
    """
    Inputs:     
    Outputs:    
    """
    writer.save()
    print "data exported to excel"
    
#a simple writer for CSVs
def csv_write(fileName, Data):
    """
    Inputs:     
    Outputs:    
    """
    import csv
    with open('%s.csv'%(fileName), 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(Data)
    print ('%s written as csv'%(fileName))

# models are defined in their own files, this is just meant as a pass 
# through to FNC so FNC can be a one stop shop
def Build_Pd_Model(C_pop, T_pop, matches, weights): 
    import FNC_model_build_Pd
    """
    Inputs:     
    Outputs:    
    """
    return FNC_model_build_Pd.Build(C_pop, T_pop, matches, weights)

def Build_TD_Model(C_pop, T_pop, matches, weights):
    import FNC_model_build_TD
    """
    Inputs:     
    Outputs:    
    """
    return FNC_model_build_TD.Build(C_pop, T_pop, matches, weights)

#allows you to select the which model you want to create 
def Build_Model(model_base, C_pop, T_pop, matches, weights):
    """
    Inputs:     
    Outputs:    
    """
    selectM = {'Pd': Build_Pd_Model,
               'TD': Build_TD_Model}
    m, assign, Timings = selectM[model_base](C_pop, T_pop, matches, weights)
    return m, assign, Timings


# outputs the average value of C_pop, T_pop, and the elements of C_pop 
# selected in by the optimization
def compare_Pd_output(C_pop, T_pop, assign, matches):
    """
    Inputs:     
    Outputs:    
    """
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
    """
    Inputs:     
    Outputs:    
    """
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
    import FNC_regression
    """
    Inputs:     
    Outputs:    
    """
    return FNC_regression.compare_MLR(C_pop, T_pop, C_matched)

def timing_data_regression(Timings):
    import FNC_regression
    """
    Inputs:     
    Outputs:    
    """
    return FNC_regression.timing_data_regression(Timings)

#get the files needed in place to run julia within the Julia folder
def setup_julia(dataSet, T_n):
    """
    Inputs:     
    Outputs:    
    """
    C, T = Import_DataSets(dataSet)
    T = T.head(T_n)
    C = C.head(T_n*30)
    a, b, dist = Pop_Calculations(C, T)
    T.to_csv("Julia/T_pop_%i.csv"%(dataSet), index = False)
    C.to_csv('Julia/C_pop_%i.csv'%(dataSet), index = False)
    csv_write('Julia/dist_%i'%(dataSet), dist)




















