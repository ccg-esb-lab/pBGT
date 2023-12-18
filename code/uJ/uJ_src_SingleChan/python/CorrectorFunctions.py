import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.patches as patches
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure


from shapely.geometry import Point,MultiPoint,MultiPolygon,Polygon,box
from shapely import wkt
from shapely import affinity
from shapely.geometry import LineString
from descartes.patch import PolygonPatch

import random
import pandas as pd

import re

from shapely import geometry
import pickle
import time

from copy import deepcopy
import skimage.measure as measure
from magicgui import magicgui
import magicgui.widgets as mwidgets
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout



import warnings
warnings.simplefilter('ignore', np.RankWarning)



from DataManagers import *
from DataStructs import *

print("CorrectorFunctions... ",end='')



def get_trap(sel,this_trap_list):
    #global current_trap
    import __main__
    if(sel=="All"):
        __main__.trap_list=this_trap_list
    else:
        tindex=this_trap_list.index(sel)
        __main__.trap_list=this_trap_list[tindex:tindex+1]
    __main__.current_trap=sel
    current_trap=sel
    
    print("Selected trap:",current_trap,end="\r")
    return current_trap

######################### Correction erros fixers

def fix_cells_with_double_lineages(df_lineages_geomsl):
    my_df=df_lineages_geomsl.copy()
    fill_cell_ids=my_df.cellID.unique()
    for this_cell in fill_cell_ids:
        celldf=my_df[my_df["cellID"]==this_cell]
        cell_lids=celldf.lineageID.unique()
        if(len(cell_lids)>1):
            last_cell='0.0'
            selected_id=-1
            for this_lid in cell_lids:
                lineage_cells=list(my_df[my_df["lineageID"]==this_lid]["cellID"])
                lineage_cells.sort()
                this_last_cell=lineage_cells[-1]
                if(this_last_cell>last_cell):
                    last_cell=this_last_cell
                    selected_id=this_lid
            print("conflicting lineages",cell_lids)
            print("Selecting",selected_id)
            for this_lid in cell_lids:
                lineage_df=my_df[my_df["lineageID"]==this_lid]
                inxs=lineage_df.index
                my_df.loc[inxs,"lineageID"]=selected_id
    
    return my_df


def dataframe_error_reporter(df_lineages_geomsl):
    my_df=df_lineages_geomsl.copy()
    track_list=my_df.trackID.unique()
    wrong_track_names=[]
    wrong_track_structure=[]
    for tt in track_list:
        track_df=my_df[my_df["trackID"]==tt]
        cells_list=list(track_df["cellID"])
        cells_listS=cells_list.copy()
        cells_listS.sort()
        tname0,tname1=tt.split("-")
        c0=cells_list[0]
        c1=cells_list[-1]
        if((tname0!=c0)or (tname1!=c1)):
            print("track not named ok",tt,cells_list)
            wrong_track_names.append(tt)
            #rename
        if(cells_listS!=cells_list):
            print("track not made ok",tt,cells_list)
            wrong_track_structure.append(tt)
            #droping
    print("Wrong structure tracks",wrong_track_structure)
    print("Wrong named tracks",wrong_track_names)
    return wrong_track_names,wrong_track_names

def dataframe_error_name_corrrector(df_lineages_geomsl,wrong_track_names,confirmation=False):
    
    my_df=df_lineages_geomsl.copy()
    track_list=wrong_track_names
    wrong_track_names=[]
    wrong_track_structure=[]
    if(confirmation):
        print("wirtes \"y\" or \"n\" to change track name...")
    else:
        print("Automatic changing wrog track names")
    for wrong_name in track_list:
        track_df=my_df[my_df["trackID"]==wrong_name]
        tname0,tname1=wrong_name.split("-")
        cells_list=list(track_df["cellID"])
        c0=cells_list[0]
        c1=cells_list[-1]
        new_track_name=c0+"-"+c1
        print("Changing track name %s to %s"%(wrong_name,new_track_name))
        if(confirmation):
            display(track_df)
            print("proceed? (y/n)")
            choise=input()
            if(choise=="n"):
                continue
            elif(choise=="y"):
                inxs=track_df.index
                my_df.loc[inxs,"trackID"]=new_track_name
        else:
            inxs=track_df.index
            my_df.loc[inxs,"trackID"]=new_track_name
    return my_df
    
def fix_splited_tracks_and_divisions(df_lineages_geomsl):
    my_df=df_lineages_geomsl.copy()
    lineages_list=[int(x) for x in list(my_df.lineageID.unique())]
    lineages_list.sort()
    print("checking lineage:")
    for this_lineage in lineages_list:
        print(this_lineage,end=", ")
        trackID_list1=list(my_df[my_df['lineageID']==this_lineage]["trackID"])
        trackID_list1.sort()
        for i,this_trackID in enumerate(trackID_list1):
            sub2_1=my_df[my_df['trackID']==this_trackID]
            if(sub2_1.shape[0]==0):
                continue
            this_cells=list(sub2_1["cellID"])
            this_first_cell=this_cells[0]
            this_last_cell=this_cells[-1]
            if(this_first_cell==this_last_cell):
                sub2_1=my_df[my_df['trackID']==this_trackID]
                my_df=my_df.drop(sub2_1.index)
                sub2_1=my_df[my_df['cellID']==this_first_cell]
                my_df.loc[sub2_1.index,"division"]=0

            for j,next_trackID in enumerate(trackID_list1[i+1:]):
                #next_trackID=second_trackid_list[i+1]
                next_first_cell=next_trackID.split("-")[0]
                next_last_cell=next_trackID.split("-")[1]
                if(this_last_cell==next_first_cell):
                    sub2_1=my_df[my_df['trackID']==this_trackID]
                    sub2_2=my_df[my_df['trackID']==next_trackID]
                    if(sub2_1.shape[0]==0 or sub2_2.shape[0]==0):
                        continue
                    print("\nJoining track %s with track %s"%(this_trackID,next_trackID))
                    sub_joint=pd.concat([sub2_1,sub2_2])
                    my_df=my_df.drop(sub_joint.index)
                    repeated_lines=sub_joint[sub_joint["cellID"]==this_last_cell]
                    sub_joint.loc[repeated_lines.index[1],"division"]=0
                    sub_joint=sub_joint.drop(repeated_lines.index[0])
                    sub_joint["trackID"]=this_first_cell+"-"+next_last_cell
                    #display(sub_joint)
                    my_df=pd.concat([my_df,sub_joint])

    my_df.reset_index(drop=True,inplace=True)
    return my_df

def check_wrong_divisions(df_lineages_geoms,lineageIDx):
    print("Checking wrong divs for lineage",lineageIDx)
    my_df=df_lineages_geoms.copy()
    lienageDf=my_df[my_df["lineageID"]==lineageIDx]
    cellids=lienageDf.cellID.unique()
    for thisc in cellids:
        celldf=my_df[my_df["cellID"]==thisc]
        cell_tracks=celldf.trackID.unique()
        cell_tracks.sort()
        track0=cell_tracks[0]
        cells0=list(my_df[my_df["trackID"]==track0]["cellID"])
        for this_track in cell_tracks:
            fisrt_cell,last_cell=this_track.split("-")
            if(thisc==fisrt_cell):
                #print("division event on tracks")
                if((celldf.shape[0]>sum(celldf.division)) &(celldf.shape[0]>>1)):
                    cells1=list(my_df[my_df["trackID"]==this_track]["cellID"])
                    arein=[True for c1 in cells1 if c1 in cells0]
                    if(sum(arein)==len(cells1)):
                        print("----------------------------------------------------->%s Not droped properly"%this_track)
                    else:
                        print("----------------------------------------------------->%s Division not made properly"%this_track)
    #return my_df
    return None

def fix_wrong_divisions(df_lineages_geoms,lineageIDx):
    my_df=df_lineages_geoms.copy()
    print("Checking wrong divs for lineage",lineageIDx)
    lienageDf=my_df[my_df["lineageID"]==lineageIDx]
    cellids=lienageDf.cellID.unique()
    for thisc in cellids:
        celldf=my_df[my_df["cellID"]==thisc]
        cell_tracks=celldf.trackID.unique()
        cell_tracks.sort()
        track0=cell_tracks[0]
        cells0=list(my_df[my_df["trackID"]==track0]["cellID"])
        for this_track in cell_tracks:
            fisrt_cell,last_cell=this_track.split("-")
            if(thisc==fisrt_cell):
                #print("division event on tracks")
                if((celldf.shape[0]>sum(celldf.division)) &(celldf.shape[0]>>1)):
                    track_df=my_df[my_df["trackID"]==this_track]
                    cells1=list(track_df["cellID"])
                    arein=[True for c1 in cells1 if c1 in cells0]
                    if(sum(arein)==len(cells1)):
                        print("%s Not droped properly"%this_track)
                        inxs=track_df.index
                        my_df=my_df.drop(inxs)
                    else:
                        print("%s Division not made properly"%this_track)
                        inxs=celldf.index
                        my_df.loc[inxs,"division"]=1
                        
                        
                        
    return my_df


####################Correction assosiated

def get_this_cell_df(mcell,nlid): 
    global frame_experiment_start,tracked_cells_corrected,expeLabel,current_trap
    get_main_vars()
    f,n=mcell.split('.')
    this_cell=tracked_cells_corrected[int(f)-frame_experiment_start][int(n)]
    ttrack=mcell+"-"+mcell
    tdivision=0
    tmother=mcell #this_cell['motherID'][0]
    tscore=0
    if(len(this_cell['tracking_score'])>0):
        tscore=this_cell['tracking_score'][0]
    if(len(this_cell['trackedBy_next_frame'])>1):
        tdivision=1
    mcell_df=pd.DataFrame.from_records([{"lineageID":nlid,"trackID":ttrack,"cellID":this_cell['cellID'],"motherID":tmother,"frame":int(f),"roiID":this_cell['roiID'],
                         "length":this_cell['axis'].length,"division":tdivision,"state":this_cell['state'],"tracking_score":tscore,"GFP":this_cell['GFP'],"DsRed":this_cell['DsRed'],
                      "cellColor":this_cell['cellColor'],"center":this_cell['center'],"roiPoly":this_cell['roiPoly'],"axis":this_cell['axis']
                     }])
    mcell_df=mcell_df.set_index(pd.Series([-1]))
    return mcell_df

def change_subtree_lineage(my_dfl,cID2,tree_tracksl,new_lineageID):
    
    print("Changing sub-tree of track %s..."%tree_tracksl)
    
    my_dft=my_dfl.copy()
    while(tree_tracksl!=[]):
        this_track=tree_tracksl[0]
        this_track_df=my_dft[(my_dft['trackID']==this_track)]
        div_cells=list(this_track_df[this_track_df["division"]==1]["cellID"])
        div_cells=[x for x in div_cells if x >=cID2]
        print("this track",this_track,"div_cells",div_cells,end="\t")
        for this_div_cell in div_cells:
            p_next_track_df=my_dft[(my_dft['cellID']==this_div_cell)&(my_dft['trackID']!=this_track)&(my_dft['lineageID']!=new_lineageID)]
            if(p_next_track_df.shape[0]==0):
                continue
            next_track_list=list(p_next_track_df.trackID.unique())
            for next_track in next_track_list:
                next_track_df=my_dft[my_dft['trackID']==next_track]
                print(next_track,end="\t")
                tree_tracksl.append(next_track)
                my_dft.loc[next_track_df.index,"lineageID"]=new_lineageID
        #print("\nPoping ",this_track)
        tree_tracksl.pop(tree_tracksl.index(this_track))
        #tree_tracksl.sort()
    print("Pending tracks:",tree_tracksl)
    print("<****")
    return my_dft

def lineages_df_rearrage_disassociation(cID1,cID2,df_lineages_geoms):
    print("Disassociating dataframe for %s and %s"%(cID1,cID2))
    #display(df_lineages_geoms.head(1))
    my_df=df_lineages_geoms.copy()
    nlineages=max(my_df.lineageID.unique())
    c_df1=my_df[my_df['cellID']==cID1]
    c_df2=my_df[my_df['cellID']==cID2]
    if((c_df1.shape[0]==0) and (c_df2.shape[0]==0)):
        c_df1=get_this_cell_df(cID1,nlineages+1)
        c_df2=get_this_cell_df(cID2,nlineages+1)
        my_df=pd.concat([my_df,c_df1])
        my_df=pd.concat([my_df,c_df2])
        my_df.reset_index(drop=True,inplace=True)
    elif((c_df1.shape[0]==0) and (c_df2.shape[0]!=0)):
        c2lid=c_df2["lineageID"].unique()[0]
        c_df1=get_this_cell_df(cID1,c2lid)
        my_df=pd.concat([my_df,c_df1])
        my_df.reset_index(drop=True,inplace=True)
    elif((c_df1.shape[0]!=0) and (c_df2.shape[0]==0)):
        c1lid=c_df1["lineageID"].unique()[0]
        c_df2=get_this_cell_df(cID2,c1lid)
        my_df=pd.concat([my_df,c_df2])
        my_df.reset_index(drop=True,inplace=True)
    nlineages=max(my_df.lineageID.unique())
    
    #nlineages=len(my_df.lineageID.unique())
    
    
    trackID_list1=list(my_df[my_df['cellID']==cID1]["trackID"])
    trackID_list2=list(my_df[my_df['cellID']==cID2]["trackID"])
    lineageID_list1=list(my_df[my_df['cellID']==cID1]["lineageID"])
    lineageID_list2=list(my_df[my_df['cellID']==cID2]["lineageID"])
    
    common_lineageID=list(set([x for x in lineageID_list1 if x in lineageID_list2]))[0]
    print("common lid:",common_lineageID)
    common_trackID_list=[x for x in trackID_list1 if x in trackID_list2]
    if(len(common_trackID_list)==0):
        print("This cells were not assosiated in dataframe")
        my_df.reset_index(drop=True,inplace=True)
        return my_df
    common_trackID=common_trackID_list[0]
    print("nlin",nlineages,"common track",common_trackID,"tracks1",trackID_list1,"tracks2",trackID_list2)
    trackID_i=trackID_list1.index(common_trackID)
    trackid1=trackID_list1[trackID_i]
    common_first_cell=common_trackID.split("-")[0]
    
    if(common_first_cell!=cID1):
        common_trackid_df=my_df[my_df['trackID']==common_trackID].copy()
        my_df=my_df.drop(common_trackid_df.index)
        first_cell1=common_trackID.split("-")[0]
        last_cell1=common_trackID.split("-")[1]
        sub_before1=common_trackid_df[common_trackid_df["cellID"]<=cID1].copy()
        new_before_trackID=first_cell1+"-"+cID1
        sub_before1["trackID"]=new_before_trackID
        my_df=pd.concat([my_df,sub_before1])    
        
        sub_after1=common_trackid_df[common_trackid_df["cellID"]>cID1].copy()
        new_after_trackID=cID2+"-"+last_cell1
        nlineages=nlineages+1
        print("New after track",new_after_trackID,nlineages)
        sub_after1["trackID"]=new_after_trackID
        sub_after1["lineageID"]=nlineages
        my_df=pd.concat([my_df,sub_after1])    
        
        tree_tracks=[new_after_trackID]
        
        my_df=change_subtree_lineage(my_df,cID2,tree_tracks,nlineages)
        
        
        starting_c1_tracks=[]
        starting_c1_lasts_cells=[]
        trackID_list1=list(my_df[(my_df['cellID']==cID1)&(my_df['lineageID']==common_lineageID)]["trackID"])
        for this_track in trackID_list1:
            this_first_cell=this_track.split("-")[0]
            this_last_cell=this_track.split("-")[1]
            if(this_first_cell==cID1):
                starting_c1_tracks.append(this_track)
                starting_c1_lasts_cells.append(this_last_cell)
        if(len(starting_c1_tracks)==1):
            #join and div=0
            afert_t=starting_c1_tracks[0]
            sub2_1=my_df[my_df['trackID']==new_before_trackID]
            sub2_2=my_df[my_df['trackID']==afert_t]
            fist_cell1,last_cell1=new_before_trackID.split("-")
            fist_cell2,last_cell2=afert_t.split("-")
            print("\nJoining track %s with track %s"%(new_before_trackID,afert_t))
            sub_joint=pd.concat([sub2_1,sub2_2])
            my_df=my_df.drop(sub_joint.index)
            repeated_lines=sub_joint[sub_joint["cellID"]==last_cell1]
            sub_joint.loc[repeated_lines.index[1],"division"]=0
            sub_joint=sub_joint.drop(repeated_lines.index[0])
            sub_joint["trackID"]=fist_cell1+"-"+last_cell2
            #display(sub_joint)
            my_df=pd.concat([my_df,sub_joint])
            
            
        elif(len(starting_c1_tracks)>1):
            #chose longest 
            maxindex=starting_c1_lasts_cells.index(max(starting_c1_lasts_cells))
            afert_t=starting_c1_tracks[maxindex]
            #join
            sub2_1=my_df[my_df['trackID']==new_before_trackID]
            sub2_2=my_df[my_df['trackID']==afert_t]
            fist_cell1,last_cell1=new_before_trackID.split("-")
            fist_cell2,last_cell2=afert_t.split("-")
            print("\nJoining track %s with track %s"%(new_before_trackID,afert_t))
            sub_joint=pd.concat([sub2_1,sub2_2])
            my_df=my_df.drop(sub_joint.index)
            repeated_lines=sub_joint[sub_joint["cellID"]==last_cell1]
            sub_joint=sub_joint.drop(repeated_lines.index[0])
            sub_joint["trackID"]=fist_cell1+"-"+last_cell2
            #display(sub_joint)
            my_df=pd.concat([my_df,sub_joint])
        #
    elif(common_first_cell==cID1):
        ##
        common_trackid_df=my_df[my_df['trackID']==common_trackID].copy()
        my_df=my_df.drop(common_trackid_df.index)
        first_cell1=common_trackID.split("-")[0]
        last_cell1=common_trackID.split("-")[1]
        sub_before1=common_trackid_df[common_trackid_df["cellID"]<=cID1].copy()
        print("Indeed before is just cell1. Length df:",sub_before1.shape[0])
        new_before_trackID=first_cell1+"-"+cID1
        sub_before1["trackID"]=new_before_trackID
        my_df=pd.concat([my_df,sub_before1])    
        
        sub_after1=common_trackid_df[common_trackid_df["cellID"]>cID1].copy()
        new_after_trackID=cID2+"-"+last_cell1
        nlineages=nlineages+1
        print("New after track",new_after_trackID,nlineages)
        sub_after1["trackID"]=new_after_trackID
        sub_after1["lineageID"]=nlineages
        my_df=pd.concat([my_df,sub_after1])    
        
        tree_tracks=[new_after_trackID]
        my_df=change_subtree_lineage(my_df,cID2,tree_tracks,nlineages)
        
        starting_c1_tracks=[]
        starting_c1_lasts_cells=[]
        trackID_list1=list(my_df[(my_df['cellID']==cID1)&(my_df['lineageID']==common_lineageID)]["trackID"])
        for this_track in trackID_list1:
            this_first_cell=this_track.split("-")[0]
            this_last_cell=this_track.split("-")[1]
            if(this_first_cell==cID1):
                starting_c1_tracks.append(this_track)
                starting_c1_lasts_cells.append(this_last_cell)
        if(len(starting_c1_tracks)==0):
            # div=0
            trackID_list1=list(my_df[(my_df['cellID']==cID1)]["trackID"])
            print("Indeed only track of cell1 . Tracks:",trackID_list1)
            this_trackID=trackID_list1[0]
            sub2_1=my_df[my_df['trackID']==this_trackID].copy()
            
            div_line=sub2_1[sub2_1["cellID"]==cID1]
            my_df.loc[div_line.index,"division"]=0
            
            #my_df=pd.concat([my_df,sub2_1])
    
    my_df.reset_index(drop=True,inplace=True)
    check_wrong_divisions(my_df,nlineages)
    check_wrong_divisions(my_df,common_lineageID)
    print("about to return")
    return my_df

def lineages_df_rearrage_association(cID1,cID2,df_lineages_geoms):
    global frame_experiment_start,tracked_cells_corrected
    print("Associating dataframe for %s and %s"%(cID1,cID2))
    my_df=df_lineages_geoms.copy()
    nlineages=max(my_df.lineageID.unique())
    c_df1=my_df[my_df['cellID']==cID1]
    c_df2=my_df[my_df['cellID']==cID2]
    if((c_df1.shape[0]==0) and (c_df2.shape[0]==0)):
        c_df1=get_this_cell_df(cID1,nlineages+1)
        c_df2=get_this_cell_df(cID2,nlineages+1)
        my_df=pd.concat([my_df,c_df1])
        my_df=pd.concat([my_df,c_df2])
        my_df.reset_index(drop=True,inplace=True)
    elif((c_df1.shape[0]==0) and (c_df2.shape[0]!=0)):
        c2lid=c_df2["lineageID"].unique()[0]
        c_df1=get_this_cell_df(cID1,c2lid)
        my_df=pd.concat([my_df,c_df1])
        my_df.reset_index(drop=True,inplace=True)
    elif((c_df1.shape[0]!=0) and (c_df2.shape[0]==0)):
        c1lid=c_df1["lineageID"].unique()[0]
        c_df2=get_this_cell_df(cID2,c1lid)
        my_df=pd.concat([my_df,c_df2])
        my_df.reset_index(drop=True,inplace=True)

    nlineages=max(my_df.lineageID.unique())
    
    trackID_list1=list(my_df[my_df['cellID']==cID1]["trackID"])
    trackID_list1.sort()
    trackID_i=0
    ntracks=len(trackID_list1)
    if(ntracks>1):
        print("Cell %s is linked to %s cells"%(cID1,ntracks))
        lastscells=[x.split("-")[1] for x in trackID_list1]
        trackID_i=lastscells.index(max(lastscells))
        print("Using %s"%trackID_list1[trackID_i])

    #trackid1=trackID_list1[trackID_i]
    trackid1=trackID_list1[0]
    trackid_df1=my_df[my_df['trackID']==trackid1].copy()
    sub_before1=trackid_df1[trackid_df1["cellID"]<cID1].copy()
    sub_after1=trackid_df1[trackid_df1["cellID"]>=cID1].copy()
    lineageID1=trackid_df1.lineageID.unique()[0]
    lineageid_df1=my_df[my_df['lineageID']==lineageID1]
    first_cell1=trackid1.split("-")[0]
    last_cell1=trackid1.split("-")[1]
    print(1,"trackids:",trackID_list1)
    print(1,trackid1,last_cell1,lineageID1)
    
    trackID_list2=list(my_df[my_df['cellID']==cID2]["trackID"])
    trackID_list2.sort()
    print(2,"trackids:",trackID_list2)
    trackid2=trackID_list2[0]
##############################################
    trackid_df2=my_df[my_df['trackID']==trackid2].copy()
    sub_before2=trackid_df2[trackid_df2["cellID"]<cID2].copy()
    sub_after2=trackid_df2[trackid_df2["cellID"]>=cID2].copy()

    lineageID2=trackid_df2.lineageID.unique()[0]
    lineageid_df2=my_df[my_df['lineageID']==lineageID2]
    first_cell2=trackid2.split("-")[0]
    last_cell2=trackid2.split("-")[1]
    print(2,trackid2,last_cell2,lineageID2)

    if(last_cell1==cID1):   ### not a division event
        print("Simple linking")
        
        #trackid_df1
        #firstcells=[x.split("-")[0] for x in trackID_list]
        if(first_cell2==cID2):
            print("a")
            print("first is",cID2,trackid2)
            new_trackID=first_cell1+"-"+last_cell2
            new_df=pd.concat([trackid_df1,trackid_df2])
            new_df["trackID"]=new_trackID
            cdf=new_df[new_df["cellID"]==cID1]
            new_df.loc[cdf.index,"division"]=0
            new_df["lineageID"]=lineageID1
            my_df=my_df.drop(new_df.index)
            my_df=pd.concat([my_df,new_df])
        else: 
            print("b")
            print("to track",trackid2)
            prevcells=list(sub_before2["cellID"])
            prevcell2=prevcells[-1]
            prevcell2f=float(prevcell2)

            f1=int(prevcell2f)
            f01=f1-frame_experiment_start
            n1=int(round(prevcell2f%1,3)*1000)
            cell1=tracked_cells_corrected[f01][n1]
            this_list=cell1["trackedBy_next_frame"]
            if(cID2 in this_list):
                index=this_list.index(cID2)
                cell1["trackedBy_next_frame"].pop(index)
                cell1["tracking_score"].pop(index)
            
            new_before2_trackID=first_cell2+"-"+prevcell2
            print("new before 2",new_before2_trackID)
            sub_before2["trackID"]=new_before2_trackID
            my_df=my_df.drop(sub_before2.index)
            my_df=pd.concat([my_df,sub_before2])

            new_trackID=first_cell1+"-"+last_cell2
            new_df=pd.concat([trackid_df1,sub_after2])
            new_df["trackID"]=new_trackID
            cdf=new_df[new_df["cellID"]==cID1]
            new_df.loc[cdf.index,"division"]=0
            new_df["lineageID"]=lineageID1
            my_df=my_df.drop(new_df.index)
            my_df=pd.concat([my_df,new_df])
            #my_df.reset_index(drop=True,inplace=True)
            

        tree_tracks=[new_trackID]
        my_df=change_subtree_lineage(my_df,cID2,tree_tracks,lineageID1)
        
    else: ### this a division event
        print("mixing tracks")
        print("c")
        c1_df=trackid_df1[trackid_df1["cellID"]==cID1].copy()  
        if(c1_df.shape[0]>1):
            inxs=c1_df.index[0]
            c1_df=c1_df[c1_df.index==inxs]
        
        c2_df=trackid_df2[trackid_df2["cellID"]==cID2].copy()  
        if(c2_df.shape[0]>1):
            inxs=c2_df.index[0]
            c2_df=c2_df[c2_df.index==inxs]
        #display(c1_df)
        #display(c2_df)
         #### First check which lineage has the longest track after joining cells
        if(sub_after2.shape[0]>sub_after1.shape[0]):
            print("Second tail is larger fusioning tracks")
            c1_df["division"]=1    
            new_first_trackID=first_cell1+"-"+last_cell2
            print("New 1 joint track",new_first_trackID)
            newdf=pd.concat([sub_before1,c1_df,sub_after2])
            newdf["trackID"]=new_first_trackID
            newdf["lineageID"]=lineageID1
            my_df=my_df.drop(newdf.index)
            my_df=pd.concat([my_df,newdf])
            #my_df.reset_index(drop=True,inplace=True)
            tree_tracks=[new_first_trackID]
            
            div_index=c1_df.index
            sub_after1=sub_after1.drop(div_index)
            
            newdfDiv=pd.concat([c1_df,sub_after1])
            new_div_trackID=cID1+"-"+last_cell1
            print("New div track",new_div_trackID)
            newdfDiv["trackID"]=new_div_trackID
            newdfDiv["lineageID"]=lineageID1
            my_df=my_df.drop(sub_after1.index)
            my_df=pd.concat([my_df,newdfDiv])
            
            second_before_cell_list=list(sub_before2["cellID"])
            if(second_before_cell_list==[]):
                print("Second trackID starts with linking cell")
                #return my_df
            else:
                
                prevcell2=second_before_cell_list[-1]
                prevcell2f=float(prevcell2)
                f1=int(prevcell2f)
                f01=f1-frame_experiment_start
                n1=int(round(prevcell2f%1,3)*1000)
                cell1=tracked_cells_corrected[f01][n1]
                this_list=cell1["trackedBy_next_frame"]
                if(cID2 in this_list):
                    index=this_list.index(cID2)
                    cell1["trackedBy_next_frame"].pop(index)
                    cell1["tracking_score"].pop(index)
                new_second_last_cell=second_before_cell_list[-1]
                new_second_first_cell=second_before_cell_list[0]
                new_second_trackID=new_second_first_cell+"-"+new_second_last_cell
                print("New 2 before",new_second_trackID)
                sub_before2["trackID"]=new_second_trackID
                my_df=my_df.drop(sub_before2.index)
                my_df=pd.concat([my_df,sub_before2])
                #my_df.reset_index(drop=True,inplace=True)
                
            #display(sub_before2)
            my_df=change_subtree_lineage(my_df,cID2,tree_tracks,lineageID1)
            
        else:
            print("Firts track is longer, joining as a division")

            new_div_trackID=cID1+"-"+last_cell2
            print("New div ",new_div_trackID)
            line1=trackid_df1[trackid_df1["cellID"]==cID1].copy()  
            div_index=line1.index
            my_df.loc[div_index,"division"]=1
            line1["division"]=1    
            newdf=pd.concat([line1,sub_after2])
            newdf["trackID"]=new_div_trackID
            newdf["lineageID"]=lineageID1
            my_df=my_df.drop(sub_after2.index)
            my_df=pd.concat([my_df,newdf])
            
            tree_tracks=[new_div_trackID]
            
            second_before_cell_list=list(sub_before2["cellID"])
            if(second_before_cell_list==[]):
                print("Second trackID starts with linking cell")
                #return my_df
            else:
                prevcell2=second_before_cell_list[-1]
                prevcell2f=float(prevcell2)
                f1=int(prevcell2f)
                f01=f1-frame_experiment_start
                n1=int(round(prevcell2f%1,3)*1000)
                cell1=tracked_cells_corrected[f01][n1]
                this_list=cell1["trackedBy_next_frame"]
                if(cID2 in this_list):
                    index=this_list.index(cID2)
                    cell1["trackedBy_next_frame"].pop(index)
                    cell1["tracking_score"].pop(index)
                 
                new_second_last_cell=second_before_cell_list[-1]
                new_second_first_cell=second_before_cell_list[0]
                new_second_trackID=new_second_first_cell+"-"+new_second_last_cell
                print("New 2 before",new_second_trackID)
                sub_before2["trackID"]=new_second_trackID
                my_df=my_df.drop(sub_before2.index)
                my_df=pd.concat([my_df,sub_before2])
                #my_df.reset_index(drop=True,inplace=True)
            
            my_df=change_subtree_lineage(my_df,cID2,tree_tracks,lineageID1)
    
    

    #trackID_list1=list(my_df[my_df['cellID']==cID1]["trackID"])
    print("check for missjoined tracks")
    for this_lineage in [lineageID1,lineageID2]:
        trackID_list1=list(my_df[my_df['lineageID']==this_lineage]["trackID"].unique())
        trackID_list1.sort()
        #print(this_lineage,trackID_list1)
        for i,this_trackID in enumerate(trackID_list1):
            sub2_1=my_df[my_df['trackID']==this_trackID]
            if(sub2_1.shape[0]==0):
                #print("---------------------------->This should never happend")
                continue
            this_cells=list(sub2_1["cellID"])
            this_first_cell=this_cells[0]
            this_last_cell=this_cells[-1]
            
            if(this_first_cell==this_last_cell):
                print("Dropping track ---> ",this_trackID)
                sub2_1=my_df[my_df['trackID']==this_trackID]
                my_df=my_df.drop(sub2_1.index)
                
                last_cell_trackID_list=list(my_df[(my_df['cellID']==this_first_cell)&(my_df['lineageID']==this_lineage)]["trackID"].unique())
                starting_cn_tracks=[]
                starting_cn_lasts_cells=[]
                for this_last_track in last_cell_trackID_list:
                    this_last_first_cell=this_last_track.split("-")[0]
                    this_last_last_cell=this_last_track.split("-")[1]
                    if(this_last_first_cell==this_last_cell):
                        starting_cn_tracks.append(this_last_track)
                        starting_cn_lasts_cells.append(this_last_last_cell)
                                                      

                if(len(starting_cn_tracks)==0):
                    sub2_1=my_df[my_df['cellID']==this_first_cell]
                    my_df.loc[sub2_1.index,"division"]=0
                                                      
                continue
            starting_cn_tracks=[]
            starting_cn_lasts_cells=[]
            last_cell_trackID_list=list(my_df[(my_df['cellID']==this_last_cell)&(my_df['lineageID']==this_lineage)]["trackID"].unique())
            trackID_listX=list(my_df[my_df['cellID']==this_last_cell]["trackID"])
            
            if(last_cell_trackID_list!=trackID_listX):
                print("woooooppssie")
            
            for this_last_track in last_cell_trackID_list:
                this_last_first_cell=this_last_track.split("-")[0]
                this_last_last_cell=this_last_track.split("-")[1]
                if(this_last_first_cell==this_last_cell):
                    starting_cn_tracks.append(this_last_track)
                    starting_cn_lasts_cells.append(this_last_last_cell)
            
            if(len(starting_cn_tracks)==1):
                #join and div=0
                print("j1")
                afert_t=starting_cn_tracks[0]
                
                sub2_1=my_df[my_df['trackID']==this_trackID]
                sub2_2=my_df[my_df['trackID']==afert_t]
                fist_cell1,last_cell1=this_trackID.split("-")
                fist_cell2,last_cell2=afert_t.split("-")
                print("\nJoining track %s with track %s"%(this_trackID,afert_t))
                sub_joint=pd.concat([sub2_1,sub2_2])
                my_df=my_df.drop(sub_joint.index)
                repeated_lines=sub_joint[sub_joint["cellID"]==last_cell1]
                sub_joint.loc[repeated_lines.index[1],"division"]=0
                sub_joint=sub_joint.drop(repeated_lines.index[0])
                sub_joint["trackID"]=fist_cell1+"-"+last_cell2
                #display(sub_joint)
                my_df=pd.concat([my_df,sub_joint])


            elif(len(starting_cn_tracks)>1):
                print("j>1")
                #chose longest 
                maxindex=starting_cn_lasts_cells.index(max(starting_cn_lasts_cells))
                afert_t=starting_cn_tracks[maxindex]
                #join
                sub2_1=my_df[my_df['trackID']==this_trackID]
                sub2_2=my_df[my_df['trackID']==afert_t]
                fist_cell1,last_cell1=this_trackID.split("-")
                fist_cell2,last_cell2=afert_t.split("-")
                print("\nJoining track %s with track %s"%(this_trackID,afert_t))
                sub_joint=pd.concat([sub2_1,sub2_2])
                my_df=my_df.drop(sub_joint.index)
                repeated_lines=sub_joint[sub_joint["cellID"]==last_cell1]
                sub_joint=sub_joint.drop(repeated_lines.index[0])
                sub_joint["trackID"]=fist_cell1+"-"+last_cell2
                #display(sub_joint)
                my_df=pd.concat([my_df,sub_joint])

    my_df.reset_index(drop=True,inplace=True)
    check_wrong_divisions(my_df,lineageID1)
    check_wrong_divisions(my_df,lineageID2)
    print("about to return")
    return my_df



######key binded fuctions




def print_roiLabel_sel_info(cellViewer):
    global frame_experiment_start,tracked_cells_corrected,df_lineages_geoms
    get_main_vars()
    
    print("This cell info...")
    l_id=0
    c_id=0.
    for layer in cellViewer.layers:
        if layer.name=="RoiLabels":
            c_id=layer.selected_label
            pos=layer.position
            #print(layer.name)
        #if "Lineage" in layer.name:
         #   l_id=layer.selected_label
    if(type(c_id)==int):
        print("Something went wrong. Click it again")
        return None

    c_ids="%.3f"%c_id
    #print("CellID:",c_ids)
    [f1,n1]=c_ids.split(".")
    f01=int(f1)-frame_experiment_start
    n1=int(n1)
    this_cell=tracked_cells_corrected[f01][n1]
    #print(this_cell)
    cell_df=df_lineages_geoms[df_lineages_geoms["cellID"]==c_ids]
    print("cellID:",c_ids,"lineageIDs",cell_df.lineageID.unique(),"trackIDs:",cell_df.trackID.unique())    
    print("tracked by next:",this_cell["trackedBy_next_frame"])
    print("tracked by prevoius:",this_cell["trackedBy_previous_frame"])
    print("score:",this_cell["tracking_score"])
    print("division:",list(cell_df["division"]))
    print("status:",this_cell["state"])
    #display(cell_df)
    return c_id

def link_these(cellViewer):
    global df_lineages_geoms,indx1,indx2,frame_experiment_start,tracked_cells_corrected,current_lineageID,current_trackID,plotwidget
    get_main_vars()
    if(type(indx1)==int or type(indx2)==int ):
        print("One or two cells did not were captured correctly")
        return None
    if(indx1>indx2):
        t=indx1
        indx1=indx2
        indx2=t
    print("\nLinking Cells...",indx1,indx2)
    f1=int(indx1)
    f01=f1-frame_experiment_start
    n1=int(round(indx1%1,3)*1000)
    n1s=str(n1).zfill(3)
    cid1str="%s.%s"%(f1,n1s)
    print(indx1,f1,f01,n1,cid1str)

    f2=int(indx2)
    f02=f2-frame_experiment_start
    n2=int(round(indx2%1,3)*1000)
    n2s=str(n2).zfill(3)
    cid2str="%s.%s"%(f2,n2s)
    temp_list=tracked_cells_corrected[f01][n1]["trackedBy_next_frame"]
    if(cid2str not in temp_list):
        df_lineages_geoms=lineages_df_rearrage_association(cid1str,cid2str,df_lineages_geoms)
        tracked_cells_corrected[f01][n1]["trackedBy_next_frame"].append(cid2str)
        tracked_cells_corrected[f01][n1]["tracking_score"].append('1000')
        tracked_cells_corrected[f02][n2]["trackedBy_previous_frame"].append(cid1str)
        ##here we link dataframe
    else:
        print("Already associciated in tracked cells")
        return None
    sync_corrected_vars()
    plotwidget=refresh_rendering(cellViewer,this_lineageID=current_lineageID,track_id=current_trackID)
    return None



def unlink_these(cellViewer):
    global df_lineages_geoms,indx1,indx2,tracked_cells_corrected,current_lineageID,frame_experiment_start,plotwidget,current_trackID
    get_main_vars()
    if(type(indx1)==int or type(indx2)==int ):
        print("One or two cells did not were captured correctly")
    if(indx1>indx2):
        t=indx1
        indx1=indx2
        indx2=t
    print("\nUnlinking Cells...",indx1,indx2)
    f1=int(indx1)
    f01=f1-frame_experiment_start
    n1=int(round(indx1%1,3)*1000)
    n1s=str(n1).zfill(3)
    cid1str="%s.%s"%(f1,n1s)
    #print(indx1,f1,f01,n1,cid1str)

    f2=int(indx2)
    f02=f2-frame_experiment_start
    n2=int(round(indx2%1,3)*1000)
    n2s=str(n2).zfill(3)
    cid2str="%s.%s"%(f2,n2s)

    cell1=tracked_cells_corrected[f01][n1]
    #display(df_lineages_geoms[df_lineages_geoms["cellID"]==cid1str])
    this_list=cell1["trackedBy_next_frame"]
    if(cid2str in this_list):
        index=this_list.index(cid2str)
        cell1["trackedBy_next_frame"].pop(index)
        cell1["tracking_score"].pop(index)
        df_lineages_geoms=lineages_df_rearrage_disassociation(cid1str,cid2str,df_lineages_geoms)
    else:
        print("Not assosiated in tracked cells...")
        return None
    cell2=tracked_cells_corrected[f02][n2]
    this_list=cell2["trackedBy_next_frame"]
    if(cid1str in this_list):
        index=this_list.index(cid1str)
        cell2["trackedBy_previous_frame"].pop(index)
    #print("1*",cell1)
    #print("1**",tracked_cells_corrected[f01][n1])
    #display(df_lineages_geoms[df_lineages_geoms["cellID"]==cid1str])
    #print("2*",cell2)
    sync_corrected_vars()
    plotwidget=refresh_rendering(cellViewer,this_lineageID=current_lineageID,track_id=current_trackID)
    return None


def get_first(cellViewer):
    global indx1,df_lineages_geoms,tracked_cells_corrected,frame_experiment_start
    layer=cellViewer.layers["RoiLabels"]
    indx1=layer.selected_label
    
    if(type(indx1)==int):
        print("Got int type, must be float. Please click again",indx1)
        indx1=-1
        return None
    #print("Cell 1 to link:",indx1)
    c_ids="%.3f"%indx1
    [f1,n1]=c_ids.split(".")
    f01=int(f1)-frame_experiment_start
    n1=int(n1)
    this_cell=tracked_cells_corrected[f01][n1]
    cell_df=df_lineages_geoms[df_lineages_geoms["cellID"]==c_ids]
    print("Cell 1 to link:",c_ids,"\tlineageIDs",cell_df.lineageID.unique(),"tracked by:",this_cell["trackedBy_next_frame"])
    return indx1

def get_second(cellViewer):
    global indx2,df_lineages_geoms,tracked_cells_corrected,frame_experiment_start
    layer=cellViewer.layers["RoiLabels"]
    indx2=layer.selected_label
    #indx2=roiLabelLayer.selected_label
    if(type(indx2)==int):
        print("Got int type, must be float. Please click again",indx2)
        indx2=-1
        return None
    #print("Cell 2 to link:",indx2)
    c_ids="%.3f"%indx2
    [f1,n1]=c_ids.split(".")
    f01=int(f1)-frame_experiment_start
    n1=int(n1)
    this_cell=tracked_cells_corrected[f01][n1]
    cell_df=df_lineages_geoms[df_lineages_geoms["cellID"]==c_ids]
    print("Cell 2 to link:",c_ids,"\tlineageIDs",cell_df.lineageID.unique(),"tracked by:",this_cell["trackedBy_next_frame"])
    return indx2

def reload_session_data(cellViewer):
    global tracked_cells_corrected,current_lineageID,current_trackID,df_lineages_geoms,backup_df_lineages,backup_cells,plotwidget
    get_main_vars()
    print("Reloading this session data...",end="\t")
    tracked_cells_corrected=deepcopy(backup_cells)
    df_lineages_geoms=deepcopy(backup_df_lineages)
    sync_corrected_vars()
    plotwidget=refresh_rendering(cellViewer,this_lineageID=current_lineageID,track_id=current_trackID)
    print("Loaded")
    return None

def reload_last_save_data(cellViewer):
    global tracked_cells_corrected,df_lineages_geoms,last_save_cells,last_save_df_lineages
    global current_lineageID,current_trackID,plotwidget
    get_main_saved_vars()
    get_main_vars()
    print("Reloading last save data...")
    tracked_cells_corrected=deepcopy(last_save_cells)
    df_lineages_geoms=deepcopy(last_save_df_lineages)
    sync_corrected_vars()
    plotwidget=refresh_rendering(cellViewer,this_lineageID=current_lineageID,track_id=current_trackID)
    print("Loaded")
    return None


def save_progress(cellViewer):
    global tracked_cells_corrected,df_lineages_geoms,last_save_cells,last_save_df_lineages
    get_main_saved_vars()
    get_main_vars()
    
    print("Saving corrected tracked cells:")
    dirNameTRACKEDCELLS=rootDir+'data_cells_tracked/'+current_trap+"/"+data_type+"/"
    cells_tracked=tracked_cells_corrected
    dirNameTRACKEDCELLS_Corrrected=rootDir+'data_cells_tracked_corrected/'+current_trap+"/"+data_type+"/"
    print(dirNameTRACKEDCELLS)
    print(dirNameTRACKEDCELLS_Corrrected)
    if not os.path.exists(dirNameTRACKEDCELLS_Corrrected):
        print("making dir lineages corrected...")
        os.makedirs(dirNameTRACKEDCELLS_Corrrected)

    file_cells= list(f for f in os.listdir(dirNameTRACKEDCELLS) if f.endswith('.pkl'))
    file_cells.sort()
    for frame_i, this_file in enumerate(file_cells):
        this_cells=cells_tracked[frame_i] 
        fileNameTrackedCells="%s%s.pkl"%(dirNameTRACKEDCELLS_Corrrected,os.path.splitext(os.path.basename(this_file))[0])    
        #print(fileNameTrackedCells)
        save_cells(this_cells, fileNameTrackedCells)

    print("Saving corrected lineages dataframe...")
    
    last_save_cells=deepcopy(tracked_cells_corrected)
    last_save_df_lineages=deepcopy(df_lineages_geoms)
    sync_save_vars()
    dirNameDataLineagesCorrected=rootDir+'data/lineages_corrected/'
    filename=dirNameDataLineagesCorrected+expeLabel+"_"+current_trap+"_lineages_all.csv"
    #this_df=df_lineages_geoms.copy()
    this_df=df_lineages_geoms[['lineageID','trackID','cellID','motherID','frame','roiID','length','division','state','tracking_score','GFP','DsRed']].copy()
    #display(this_df.head(1))
    print("Saving: ",filename)
    this_df.to_csv(filename,index=False)
    return None

def change_cell_status(cellViewer):
    global df_lineages_geoms
    get_main_vars()
    status_options=[-1,0,1,2]
    l_id=0
    c_id=0.
    for layer in cellViewer.layers:
        if layer.name=="RoiLabels":
            c_id=layer.selected_label
            pos=layer.position
        if "Lineage" in layer.name:
            l_id=layer.selected_label

    if(type(c_id)==int):
        print("Something went wrong. Click it again")
        return None

    c_ids="%.3f"%c_id
    print(c_ids,l_id)
    [f1,n1]=c_ids.split(".")

    f01=int(f1)-frame_experiment_start
    n1=int(n1)
    cell_df=df_lineages_geoms[df_lineages_geoms["cellID"]==c_ids]
    cell_index=cell_df.index
    current_status=list(cell_df['state'])[0]
    status_index=status_options.index(current_status)
    status_index=status_index+1
    if(status_index>=len(status_options)):
        status_index=0
    new_status=status_options[status_index]
    tracked_cells_corrected[f01][n1]['state']=new_status
    df_lineages_geoms.loc[cell_index,'state']=new_status
    print("Cell %s status changed from %s to %s"%(c_ids,current_status,new_status))
    
    sync_corrected_vars()
    return None

def select_thisCell_lineage(cellViewer):
    global frame_experiment_start,tracked_cells_corrected,df_lineages_geoms,current_lineageID,current_trackID,lineage_slider
    get_main_vars()
    from __main__ import lineage_slider 
    print("Selecting this cell lineage...")
    l_id=0
    c_id=0.
    for layer in cellViewer.layers:
        if layer.name=="RoiLabels":
            c_id=layer.selected_label
            pos=layer.position
        
    if(type(c_id)==int):
        print("Something went wrong. Click it again")
        return None

    c_ids="%.3f"%c_id
    cell_df=df_lineages_geoms[df_lineages_geoms["cellID"]==c_ids]
    print("cellID:",c_ids,"lineageIDs",cell_df.lineageID.unique(),"trackIDs:",cell_df.trackID.unique())    
    this_lineage=cell_df.lineageID.unique()[0]
    print(this_lineage)
    lineage_slider.lineage_IDw2.value=int(this_lineage)
    lineage_slider.lineage_IDw3.value=this_lineage
    current_lineageID=this_lineage
    sync_lineage_slider()
    #lineage_slider()  #uncomment this line to automatic render when choosing new lineage
    return None

########### wigdets  and rendering

def get_cmap_Lineages(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


def clear_my_layers(cellViewer,shapeM):
    data_null=np.empty(shapeM)
    tnames=['RoiLabels','Lineage','Track']
    for i in range(len(cellViewer.layers)):
        for il,layer in enumerate(cellViewer.layers):
            for tname in tnames:
                if tname in layer.name:
                    cellViewer.layers[layer.name].data=data_null
                    cellViewer.layers[layer.name].color=dict()
                    cellViewer.layers[layer.name].refresh()
                    cellViewer.layers.pop(il)


def get_lineage_mask_lid(empty_mask,df_lineages_geoms,frame_experiment_start,lid):
    matrix_shape = empty_mask.shape[1:]
    print()
    
    lineage_mask=deepcopy(empty_mask)
    lineage_df=df_lineages_geoms[df_lineages_geoms["lineageID"]==lid]
    for row in lineage_df.itertuples():#(index=False):
        this_poly=row.roiPoly
        napari_coords=np.array([[int(y),int(x)] for (x,y) in this_poly.exterior.coords])
        mask = measure.grid_points_in_poly(matrix_shape, napari_coords)
        frame=row.frame-frame_experiment_start
        lineage_mask[frame ][mask] = lid
    return lineage_mask

def get_lineage_mask_cid(empty_mask,df_lineages_geoms,frame_experiment_start,lid):
    matrix_shape = empty_mask.shape[1:]
    print()
    
    lineage_mask=deepcopy(empty_mask)
    lineage_df=df_lineages_geoms[df_lineages_geoms["lineageID"]==lid]
    for row in lineage_df.itertuples():#(index=False):
        this_poly=row.roiPoly
        cell_id=float(row.cellID)
        napari_coords=np.array([[int(y),int(x)] for (x,y) in this_poly.exterior.coords])
        mask = measure.grid_points_in_poly(matrix_shape, napari_coords)
        frame=row.frame-frame_experiment_start
        lineage_mask[frame ][mask] = cell_id
    return lineage_mask



def get_track_mask(empty_mask,df_lineages_geoms,frame_experiment_start,tid):
    matrix_shape = empty_mask.shape[1:]
    #df_lineages_geoms=trap_df_lineages_geoms[i]
    track_mask=deepcopy(empty_mask)
    track_df=df_lineages_geoms[df_lineages_geoms["trackID"]==tid]
    #display(track_df)
    for row in track_df.itertuples():#(index=False):
        this_poly=row.roiPoly
        cell_id=float(row.cellID)
        #print()
        napari_coords=np.array([[int(y),int(x)] for (x,y) in this_poly.exterior.coords])
        mask = measure.grid_points_in_poly(matrix_shape, napari_coords)
        frame=row.frame-frame_experiment_start
        track_mask[frame ][mask] = cell_id
    return track_mask

# def widget_plot_this_lineage_ID(this_lineageID,df_lineages_geoms):
#     ncolors=50
#     cmapL = cm.get_cmap("Reds", ncolors)
#     maxred=120
#     cmapN=cm.get_cmap("Greens", ncolors)
#     maxgreen=700
    
#     lineage_df=df_lineages_geoms[df_lineages_geoms['lineageID']==this_lineageID]
#     trackID_list=lineage_df.trackID.unique()
#     lineageG=nx.Graph()
#     lineage_nodeSizes=dict()
    
#     for this_trackID in trackID_list:
#         track_df=lineage_df[lineage_df['trackID']==this_trackID]
#         cellID_list=track_df.cellID.tolist()
#         #print("track:",this_trackID,cellID_list)
#         for this_cellID in cellID_list:
#             lineageG.add_node(this_cellID)
#             this_rfp=track_df[track_df["cellID"]==this_cellID]["DsRed"]
#             this_gfp=track_df[track_df["cellID"]==this_cellID]["GFP"]
#             this_color=track_df[track_df["cellID"]==this_cellID]["cellColor"]
#             lineage_nodeSizes[this_cellID]=this_rfp
#             lineageG.add_node(this_cellID,gfp=this_gfp,cellcolor=this_color)
#         for i in range(len(cellID_list)-1):
#             node1=cellID_list[i]
#             node2=cellID_list[i+1]
#             this_rfp=track_df[track_df["cellID"]==cellID_list[i+1]]["DsRed"]
#             lineageG.add_edge(node1,node2,rfp=this_rfp)
    
#     #node_colorsGFP=cmapN([(n[1]["gfp"].item())/maxgreen for n in lineageG.nodes(data=True)])
#     node_colors=[(n[1]["cellcolor"].item()) for n in lineageG.nodes(data=True)]
#     node_colors=mpl.colors.to_rgba_array(node_colors)
#     ecolors =cmapL([(lineageG[u][v]['rfp'].item())/maxred for u,v in lineageG.edges()])
#     positions=graphviz_layout(lineageG, prog='dot')
#     with plt.style.context('dark_background'):
#         #plt.rcParams['figure.dpi'] = 100
#         mpl_widget = FigureCanvas(Figure(figsize=(2, 6)))
#         static_ax = mpl_widget.figure.subplots()
#         nx.draw(lineageG,positions,ax=static_ax,nodelist=lineageG.nodes(), node_size=[x*1 for x in lineage_nodeSizes.values()],node_color=node_colors,edge_color=ecolors, alpha=0.75,  with_labels=True)
#         mpl_widget.figure.tight_layout()
#     return mpl_widget
def widget_plot_this_lineage_ID(this_lineageID,df_lineages_geoms):
    rs=df_lineages_geoms.DsRed.unique()
    gs=df_lineages_geoms.GFP.unique()
    
    ncolors=50
    #cmapL = cm.get_cmap("Reds", ncolors)
    cmapL = cm.get_cmap("RdBu_r", ncolors)
    #maxred=120
    maxred=max(rs)
    minred=min(rs)
    cmapN=cm.get_cmap("Greens", ncolors)
    #maxgreen=700
    maxgreen=max(gs)
    mingreen=min(gs)
    
    lineage_df=df_lineages_geoms[df_lineages_geoms['lineageID']==this_lineageID]
    trackID_list=lineage_df.trackID.unique()
    lineageG=nx.DiGraph()
    lineage_nodeSizes=dict()
    
    for this_trackID in trackID_list:
        track_df=lineage_df[lineage_df['trackID']==this_trackID]
        cellID_list=track_df.cellID.tolist()
        #print("track:",this_trackID,cellID_list)
        for this_cellID in cellID_list:
            lineageG.add_node(this_cellID)
            this_rfp=track_df[track_df["cellID"]==this_cellID]["DsRed"].item()
            this_gfp=track_df[track_df["cellID"]==this_cellID]["GFP"].item()
            this_size=track_df[track_df["cellID"]==this_cellID]["length"].item()
            this_color=track_df[track_df["cellID"]==this_cellID]["cellColor"].item()
            lineage_nodeSizes[this_cellID]=this_size
            lineageG.add_node(this_cellID,gfp=this_gfp,cellcolor=this_color)
        for i in range(len(cellID_list)-1):
            node1=cellID_list[i]
            node2=cellID_list[i+1]
            this_rfp=track_df[track_df["cellID"]==cellID_list[i+1]]["DsRed"].item()
            lineageG.add_edge(node1,node2,rfp=this_rfp)
    
    #node_colorsGFP=cmapN([(n[1]["gfp"])/maxgreen for n in lineageG.nodes(data=True)])
    #node_colors=node_colorsGFP
    node_colors=[(n[1]["cellcolor"]) for n in lineageG.nodes(data=True)]
    #print(lineageG.nodes())
    ecolors =cmapL([(r['rfp']-minred)/(maxred-minred) for n1,n2,r in lineageG.edges(data=True)])
    esizes=[3*(r['rfp'])/maxred for n1,n2,r in lineageG.edges(data=True)]
    positions=graphviz_layout(lineageG, prog='dot') #neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps, sccmap, tred, sfdp, unflatten.
    with plt.style.context('dark_background'):
        #plt.rcParams['figure.dpi'] = 100
        mpl_widget = FigureCanvas(Figure(figsize=(2, 6)))
        static_ax = mpl_widget.figure.subplots()
        #nx.draw(lineageG,positions,ax=static_ax,nodelist=lineageG.nodes(), node_size=[x*1 for x in lineage_nodeSizes.values()],node_color=node_colors,edge_color=ecolors, alpha=0.75,  with_labels=True)
        nx.draw(lineageG,positions,ax=static_ax,nodelist=lineageG.nodes(), node_size=[x*2 for x in lineage_nodeSizes.values()],node_color=node_colors,edge_color=ecolors,width=esizes,arrows=False, alpha=1,  with_labels=True)
        mpl_widget.figure.tight_layout()
    return mpl_widget

def refresh_rendering(cellViewer,this_lineageID=1,track_id="None"):
    global plotwidget,df_lineages_geoms
    get_main_vars()
    clear_my_layers(cellViewer,empty_mask.shape)
    plotwidget=refresh_plot_widget(this_lineageID,df_lineages_geoms,cellViewer)
    sync_plotwidget()
    roi_mask0=roi_mask.copy()
    #lineage_mask=get_lineage_mask_lid(empty_mask,df_lineages_geoms,frame_experiment_start,this_lineageID)
    lineage_mask=get_lineage_mask_cid(empty_mask,df_lineages_geoms,frame_experiment_start,this_lineageID)
    
    if(track_id=="None"):
        0
    elif(track_id=="all"):
        trackids=df_lineages_geoms[df_lineages_geoms["lineageID"]==this_lineageID].trackID.unique()
        for this_trackID in trackids:
            tname="Track_%s"%this_trackID
            track_mask=get_track_mask(empty_mask,df_lineages_geoms,frame_experiment_start,this_trackID)
            trackLayer=cellViewer.add_labels(data=track_mask,color=colorsCid,name=tname,blending="additive",opacity=.95)
            trackLayer.selected=False   
    elif(track_id):
        this_trackID=track_id
        tname="Track_%s"%this_trackID
        track_mask=get_track_mask(empty_mask,df_lineages_geoms,frame_experiment_start,this_trackID)
        trackLayer=cellViewer.add_labels(data=track_mask,color=colorsCid,name=tname,blending="additive",opacity=.95)
        
    lname="Lineage_%s"%this_lineageID
    #lineageLayer=cellViewer.add_labels(data=lineage_mask,color=colorsLid,name=lname,blending="additive",opacity=.95)
    #lineageLayer=cellViewer.add_labels(data=lineage_mask,color=colorsCid,name=lname,blending="additive",opacity=1.)
    #roiLabelLayer=cellViewer.add_labels(data=roi_mask,color=colorsCid,name="RoiLabels",blending="translucent",opacity=.3)
    lineageLayer=cellViewer.add_labels(data=lineage_mask,color=colorsCid,name=lname,blending="opaque",opacity=1.)
    roiLabelLayer=cellViewer.add_labels(data=roi_mask0,color=colorsCid,name="RoiLabels",blending="additive",opacity=.2)
    
    ## Needed for correct label/cellID picking
    lineageLayer.selected=False   
    
    roiLabelLayer.selected=True
    roiLabelLayer.mode="pick"   

    
    cellViewer.window.qt_viewer.console.execute()
    
    cellViewer.window.activate()
    ndims_pos=cellViewer.dims.current_step
    cellViewer.dims.set_current_step(0,ndims_pos[0])
    cellViewer.dims._increment_dims_left()
    cellViewer.dims._increment_dims_right()
    return plotwidget

    
def refresh_plot_widget(lineage_ID,df_lineages_geoms,cellViewer):
    global plotwidget
    
    if(plotwidget):
        cellViewer.window.remove_dock_widget(plotwidget)
    mpl_widget=widget_plot_this_lineage_ID(lineage_ID,df_lineages_geoms)
    plotwidget = cellViewer.window.add_dock_widget(mpl_widget,name='Graph plot',area='right')
    return plotwidget


############## data and vars management

def sync_corrected_vars():
    global tracked_cells_corrected,df_lineages_geoms
    import __main__
    __main__.tracked_cells_corrected=deepcopy(tracked_cells_corrected)
    __main__.df_lineages_geoms=deepcopy(df_lineages_geoms)
    return None

def sync_save_vars():
    global last_save_cells,last_save_df_lineages
    import __main__
    __main__.last_save_cells=deepcopy(last_save_cells)
    __main__.last_save_df_lineages=deepcopy(last_save_df_lineages)
    return None


def sync_plotwidget():
    global plotwidget
    import __main__
    __main__.plotwidget=plotwidget
    return None

def sync_lineage_slider():
    global slider_default_val,current_lineageID,current_trackID,lineage_slider
    import __main__
    __main__.current_lineageID=current_lineageID
    __main__.slider_default_val=slider_default_val
    __main__.lineage_slider=lineage_slider
    __main__.current_trackID=current_trackID
    
    return None





def get_main_saved_vars():
    global last_save_cells,last_save_df_lineages
    from __main__ import last_save_cells,last_save_df_lineages
    
    return None

def get_main_vars():
    global cellViewer, plotwidget,lineage_slider 
    global df_lineages_geoms, tracked_cells_corrected,current_lineageID,current_trackID,current_trap
    global roi_mask, empty_mask, frame_experiment_start, colorsCid, colorsLid, colorsCid, colorsLid
    global backup_df_lineages,backup_cells,rootDir,expeLabel,data_type
    from __main__ import cellViewer, plotwidget#,lineage_slider 
    from __main__ import df_lineages_geoms, tracked_cells_corrected,current_lineageID,current_trackID,current_trap
    from __main__ import roi_mask, empty_mask, frame_experiment_start, colorsCid, colorsLid, colorsCid, colorsLid
    from __main__ import backup_df_lineages,backup_cells,rootDir,expeLabel,data_type
    
    return None



######################## widgets ##########################
slider_default_val=-1
current_lineageID=-1
current_trackID=0
plotwidget=False
lineage_slider=False



def add_lineage_slider_widget(cellViewer):
    global df_lineages_geoms,slider_default_val,current_lineageID,lineage_slider
    get_main_vars()
    slider_default_val=1
    current_lineageID=1
    lineages_list=[int(x) for x in list(df_lineages_geoms.lineageID.unique())]
    lineages_list.sort()
    n_lineages=max(lineages_list)
    print("Current number of lineages",n_lineages)
    tracks_list=list(df_lineages_geoms[df_lineages_geoms["lineageID"]==1].trackID.unique())
    tracks_list=["None"]+tracks_list #["None","all"]+tracklist
    
    @magicgui(
    auto_call=False,   ##### autocall = True gets double clicks
    call_button="render!",
    #lineage_IDw={"widget_type": "Slider", "min":1, "max": n_lineages},
    lineage_IDw2={"choices": lineages_list},
    lineage_IDw3={"widget_type": "SpinBox", 'max':n_lineages},  
    trackID_IDw={"choices": tracks_list},
    )
    def lineage_slider(lineage_IDw2=slider_default_val,lineage_IDw3=slider_default_val,trackID_IDw=tracks_list[0]):
        global slider_default_val,current_lineageID,current_trackID,plotwidget,df_lineages_geoms,lineages_list

        if(lineage_IDw2!=slider_default_val):
            trackID_IDw="None"
            lineage_IDw3=lineage_IDw2
            lineage_slider.lineage_IDw3.value=lineage_IDw2
        elif(lineage_IDw3!=slider_default_val):
            trackID_IDw="None"
            lineage_IDw2=lineage_IDw3
            lineage_slider.lineage_IDw2.value=lineage_IDw3

        print("Selecting lineage %s     "%lineage_IDw2,end="\n")

        tracklist=list(df_lineages_geoms[df_lineages_geoms["lineageID"]==lineage_IDw2].trackID.unique())
        tracklist=["None"]+tracks_list #["None","all"]+tracklist
        lineage_slider.trackID_IDw.choices=tracklist
        lineages_list=[int(x) for x in list(df_lineages_geoms.lineageID.unique())]
        lineages_list.sort()
        lineage_slider.lineage_IDw2.choices=lineages_list   ###### this line does not always take effect
        n_lineages=max(lineages_list)
        lineage_slider.lineage_IDw3.max=n_lineages


        #print(type(plotwidget),default_val,lineage_IDw,lineage_IDw2,lineage_IDw3)
        if(lineage_IDw2 not in lineages_list):
            print("Lineage does not exists: %s"%lineage_IDw2)

            slider_default_val=lineage_IDw3+1
            current_lineageID=lineage_IDw3+1
            #print(current_lineageID,lineage_IDw2,slider_default_val)
            return None
        slider_default_val=lineage_IDw2
        current_lineageID=lineage_IDw2
        current_trackID=trackID_IDw
        #print(slider_default_val,current_lineageID,current_trackID,lineage_IDw2)
        sync_lineage_slider()
        plotwidget=refresh_rendering(cellViewer,this_lineageID=current_lineageID,track_id=current_trackID)
        return None
    mywd=cellViewer.window.add_dock_widget(lineage_slider,name="Lineage Selector",area='bottom')
    lineage_slider.trackID_IDw.label="Track ID"
    lineage_slider.lineage_IDw2.label="Lineage list"
    lineage_slider.lineage_IDw3.label="Lineage #"
    
    sync_lineage_slider()
    print("*")
    return None



# def add_save_buttons(cellViewer):
#     @magicgui(call_button='Save progress (s)')
#     def save_widget():
#         save_progress(cellViewer)
#         return 
#     @magicgui(call_button='restart data (q)')
#     def reload_session_widget():
#         reload_session_data(cellViewer)
#         return 
#     @magicgui(call_button='reload from save (w)')
#     def reload_last_save_widget():
#         reload_last_save_data(cellViewer)
#         return 


#     cellViewer.window.add_dock_widget(save_widget, area='left',name="save button")
#     cellViewer.window.add_dock_widget(reload_last_save_widget, area='left',name="reload button")
#     cellViewer.window.add_dock_widget(reload_session_widget, area='left',name="restart data")
#     return None

def add_save_buttons(cellViewer):
    @magicgui(call_button='Save progress (s)')
    def save_widget():
        save_progress(cellViewer)
        return None
    @magicgui(call_button='restart data (q)')
    def reload_session_widget():
        reload_session_data(cellViewer)
        return None
    @magicgui(call_button='reload from save (w)')
    def reload_last_save_widget():
        reload_last_save_data(cellViewer)
        return None
    
    @magicgui(call_button='help!')
    def print_help():
        text="Key binded fuctions:\n\
        To link cells  choose them indepenedenlty with \"r\" and \"y\" keys, then press \"t\" to associate\n\
        To unlink cells  choose them indepenedenlty with \"r\" and \"y\" keys, the press \"g\" to disassociate\n\
        Show this cell info \"h\"\n\
        Change cell status  \"k\"\n\
        select this cell lineage \"o\"\n\
        save progress \"s\"\n\
        reload last saved data \"w\"\n\
        reload whole session data \"q\"\n\
        \n Hints:\n\
        Hold space bar key to transiently  activate zoom & pan\n\
        Ctrl + mouse scroll chage dims-slider. Or left right arrows\n\
        To delete corrections, delete  tracked_cells_corrected/\"trap\" folder\n\
        "
        print(text)
        return None
    
    data_saving_widget=mwidgets.Container(name='Data saving', label=None, tooltip="Careful!", visible=None, enabled=True,   layout='vertical', widgets=(save_widget,reload_last_save_widget,reload_session_widget,print_help), labels=False)
    #cellViewer.window.add_dock_widget(save_widget, area='left',name="save button")
    #cellViewer.window.add_dock_widget(reload_last_save_widget, area='left',name="reload button")
    #cellViewer.window.add_dock_widget(reload_session_widget, area='left',name="restart data")
    
    cellViewer.window.add_dock_widget(data_saving_widget, area='left',name="save button")
    return None


def LoadKeyBindedFuctions(cellViewer):
    
    cellViewer.bind_key('k',change_cell_status,overwrite=True)
    cellViewer.bind_key('q',reload_session_data,overwrite=True)
    cellViewer.bind_key('w',reload_last_save_data,overwrite=True)
    cellViewer.bind_key('s',save_progress,overwrite=True)

    cellViewer.bind_key('g',unlink_these,overwrite=True)
    cellViewer.bind_key('t',link_these,overwrite=True)
    cellViewer.bind_key('r',get_first,overwrite=True)
    cellViewer.bind_key('y',get_second,overwrite=True)
    
    cellViewer.bind_key('o',select_thisCell_lineage,overwrite=True)
    cellViewer.bind_key('h',print_roiLabel_sel_info,overwrite=True)
    
    return None


print("loaded!")


