import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys
import pathlib
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import scipy.stats as st
from IPython.display import HTML, display
from IPython.display import set_matplotlib_formats





#from DataManagers import *
def get_long_lineages(this_cells, minFrames):
    complete_cells=[]
    for this_cell in this_cells:
        try:
            if len(this_cell['roiFrames'])>minFrames:
                complete_cells.append(this_cell)
        except TypeError:
            continue
    return complete_cells

def get_analysis_lineages(this_cells,an_start,an_end,frame_experiment_start):
    
    complete_cells=[]
    inif=an_start
    endf=an_end
    #print(inif,endf)
    for this_cell in this_cells:
        this_frames=this_cell['roiFrames']
        #print(inif,endf,this_frames)
        if(inif in this_frames and endf in this_frames) :
            complete_cells.append(this_cell)
    
    return complete_cells





def fromFileData(fileName):
    #d = np.loadtxt(fileName, delimiter="\t")
    
    f = open(fileName,'r')

    data = []
    numCells=0
    for line in f.readlines():
        if numCells>0:
            data.append(line.replace('\n','').split('\t'))
        numCells=numCells+1
    f.close()

    return data

#Returns a list with index of data_frame in array(t)
def filterFrames(data_frame, t):
    i=0
    ret=[]
    for this_data in data_frame:
        if this_data in t:
            ret.append(i)
        i=i+1
    return ret

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


def toFileData(fileName, row_data):
    f = open(fileName, 'w')
    f.write("id\timageGFP\timageDsRed\tpos\tframe\troi_label\tGFP\tDsRed\tNorm_GFP\tNorm_Red")
    
    for row in row_data:
        str_row="\t".join(map(str,row))
        f.write("\n%s"%str_row)
    f.close()
    


#Returns a list with index of data_frame in array(pos)
def filterPos(data_pos, pos):
    i=0
    ret=[]
    for this_data in data_pos:
        if this_data==pos:
            ret.append(i)
        i=i+1
    return ret

print("Analysis Functions loaded")