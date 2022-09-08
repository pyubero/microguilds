# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:28:47 2022

@author: logslab
"""

import os
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt

# Function to calculate Roman values
def intToRoman(num):
  
    # Storing roman values of digits from 0-9
    # when placed at different places
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D",
         "DC", "DCC", "DCCC", "CM "]
    x = ["", "X", "XX", "XXX", "XL", "L",
         "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V",
         "VI", "VII", "VIII", "IX"]
  
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
  
    ans = (thousands + hundreds +
           tens + ones)
  
    return ans
    

def nanToZero(array):
    idx = np.isnan(array)
    array[idx]=0
    return array
    
    
WIDTH = 0.25  
SLABELS= ['Epipelagic','Mesopelagic','Bathypelagic']
COLORS = plt.cm.viridis(np.linspace(0, 1, 5))[-3::-1,:]
R_STEP = 1.0
SHOW_EMPTY_BARS = False
EXPORT_FIGURES  = True

FILENAME_LIST = ['amoA_rawlogcounts.csv',
              'amt_guildvector.csv',
              'hzsA_rawlogcounts.csv',
               'nosZ_guildvector.csv',
              'nrt_guildvector.csv',
              'nxr_guildvector.csv']

# Uncomment this to process all files in the current folder
# ... that are .csv
#files = os.listdir()
#FILENAME_LIST = [file  for file in files if file.split('.')[-1] =='csv']




for FILENAME in FILENAME_LIST:
    DATA = np.loadtxt(FILENAME, delimiter=',', encoding="utf-8-sig")

    nrows, ncols = DATA.shape
    labels = [ 'F. '+intToRoman(jj+1) for jj in range(ncols) ]
    
    
    
    # Compute pie slices
    theta = np.linspace(0.0, 2 * np.pi, ncols, endpoint=False)
    rlims = np.arange(0, DATA.max(), R_STEP)
    
    plt.figure( figsize=(16,6), dpi=300)
    
    for jj in range(nrows):
        ax = plt.subplot(1,3,jj+1,projection='polar')
        
        ax.bar(theta, DATA[jj,:], width=WIDTH, color= COLORS[jj,:])
        

        if SHOW_EMPTY_BARS:
            if jj!=0:
                ax.bar(theta, DATA[0,:], width=WIDTH, edgecolor= COLORS[0,:], fill=False, linewidth=1.5)
            if jj!=1:
                ax.bar(theta, DATA[1,:], width=WIDTH, edgecolor= COLORS[1,:], fill=False, linewidth=1.5)
            if jj!=2:
                ax.bar(theta, DATA[2,:], width=WIDTH, edgecolor= COLORS[2,:], fill=False, linewidth=1.5)
            
        ax.set_rgrids( rlims)
        ax.set_thetagrids( theta*57.3, labels=labels)
        plt.title( SLABELS[jj])
        ax.xaxis.grid(linewidth=0)
        ax.yaxis.grid(linewidth=0.5)
    
        
    plt.suptitle(FILENAME)
    if EXPORT_FIGURES:
        plt.savefig( FILENAME.split('.')[0]+'.png', facecolor='w', dpi=600)
