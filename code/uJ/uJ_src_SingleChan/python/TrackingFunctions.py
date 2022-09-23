import sys
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from shapely.geometry import Point,MultiPoint,MultiPolygon,Polygon,box
from shapely import wkt
from shapely import affinity
from shapely.geometry import LineString
from descartes.patch import PolygonPatch
import ipywidgets as widgets
import ipywidgets.widgets.interaction
#from ipywidgets import interactive, fixed
from ipywidgets import *
import random
import pandas as pd

import re

from shapely import geometry
import pickle
import time

from ipyparallel import Client

import warnings
warnings.simplefilter('ignore', np.RankWarning)


#sys.path.insert(0, './lib/')
#from readroi import read_roi_zip
#from readroi import *
from DataManagers import *
from DataStructs import *

print("TrackingFunctions... 1",end='')




############### ROI PROCESSING FUNCTIONS (PARALLEL) ##################

def par_process_rois():
    
    import numpy as np
    
    def get_max_x_angle(poly):
        rangle=0
        angles= list(range(0,185,5))
        (bminx, bminy, bmaxx, bmaxy)=poly.bounds
        xl=bmaxx-bminx
        yl=bmaxy-bminy
        center=Point((xl/2)+bminx,(yl/2)+bminy)

        max_xrange=xl
        for i,angle in enumerate(angles):
            #print(i,angle)
            poly_t=affinity.rotate(poly,angle,center)
            (bminx, bminy, bmaxx, bmaxy)=poly_t.bounds
            xl=bmaxx-bminx

            if(xl>max_xrange):
                max_xrange=xl
                rangle=angle

        return rangle,center
    
    def fix_poly(p):
        buff=0
        while(p.is_valid==False or p.geom_type is not "Polygon"):
            p=p.buffer(buff)
            p=wkt.loads(wkt.dumps(p, rounding_precision=2)).simplify(0)
            #print(buff)
            buff+=.001

        return p


    def get_cell_axis(poly):
      #  t=time.time()
        line=LineString()
        pre_cloud=[]
        angle,center=get_max_x_angle(poly)
        poly=affinity.rotate(poly, angle,center )
        poly=fix_poly(poly)
        ext=poly.exterior.coords

        decs=0
        for i,seed_point in enumerate(ext):
            pre_cloud_i=[]
            seed_point=Point(seed_point)
            i=.5;
            intersection=poly.exterior.intersection(seed_point)
            sphere=seed_point.buffer(i)
            intersection=poly.exterior.intersection(sphere.exterior)
            #while(len(list(intersection))>0):
            while(intersection.is_empty==False):
                arc=Point()
                sphere=seed_point.buffer(i)
                intersection=poly.exterior.intersection(sphere.exterior)
                i+=3
                #i+=.5
    #             if(len(list(intersection))<=1):
    #                 continue
                if(intersection.is_empty):
                    continue

                arc=poly.intersection(sphere.exterior)

                if(arc.geom_type=="MultiLineString"):

                    for this_geom in arc.geoms:
                        arc_xy=list(this_geom.coords)

                        for pi in arc_xy:
                            pi=Point(pi)
                            pre_cloud_i.append((round(pi.x,decs),round(pi.y,decs)))

                elif(arc.geom_type=="LineString"):
                    for pi in list(arc.coords):
                        pi=Point(pi) 
                        pre_cloud_i.append((round(pi.x,decs),round(pi.y,decs)))

            pre_cloud.append(pre_cloud_i)

        cloud=MultiPoint(pre_cloud[0])

        for line in pre_cloud:

            line=MultiPoint(line)
            cloud=cloud.union(line)

        ###make poly fit to points cloud
        xs=[]
        ys=[]
        for pi in cloud:
            pi=Point(pi)

            xs.append(pi.x)
            ys.append(pi.y)

        pf=np.poly1d(np.polyfit(xs,ys,25))
        xr=np.linspace(np.min(xs),np.max(xs),100)

        line=[]
        for x in xr:
            pi=Point(x,pf(x))
            line.append(pi)

        line=LineString(line)
        linexy=np.array(line)
        linex=linexy[:,0]
        liney=linexy[:,1]

        line=poly.intersection(line)
        liner=affinity.rotate(line,-angle,center)
        if(liner.geom_type=="MultiLineString"):
            line=LineString()
            for this_geom in liner.geoms:
                if(this_geom.length>line.length):
                    line=this_geom
            liner=line

        linexy=np.array(liner)
        linex=linexy[:,0]
        liney=linexy[:,1]

        x=(np.max(xs)-np.min(xs))/2+np.min(xs)
        pi=Point(x,pf(x))

        center=affinity.rotate(pi,-angle,center)

        return liner,center
    
    frame_axis=[]
    frame_centers=[]
    
    for j,this_roi in enumerate(rois):
        
        this_axis,this_center=get_cell_axis(this_roi)
            
        frame_axis.append(this_axis)
        frame_centers.append(this_center)
    ret=[rois, roiIDs, frame_axis, frame_centers]
    return ret




#########################################################3

def get_max_x_angle(poly):
        rangle=0
        angles= list(range(0,185,5))
        (bminx, bminy, bmaxx, bmaxy)=poly.bounds
        xl=bmaxx-bminx
        yl=bmaxy-bminy
        center=Point((xl/2)+bminx,(yl/2)+bminy)

        max_xrange=xl
        for i,angle in enumerate(angles):
            #print(i,angle)
            poly_t=affinity.rotate(poly,angle,center)
            (bminx, bminy, bmaxx, bmaxy)=poly_t.bounds
            xl=bmaxx-bminx

            if(xl>max_xrange):
                max_xrange=xl
                rangle=angle

        return rangle,center





############### TRACKING FUNCTIONS ##################
def track_this_indexs(local_cells_index, local_cells_frames, local_cells, weight_vec):
    n_cells=len(local_cells_index)
    for i, this_index in enumerate(local_cells_index):
        to_frame=0
        from_frame=local_cells_frames[i]
        
        print("%s/%s"%(i+1,n_cells),end="\r")
        trackID=local_cells[local_cells_frames[i]][this_index]['cellID']
        local_cells[local_cells_frames[i]][this_index]['trackID'].append(trackID)    #change to something clever
        local_cells[local_cells_frames[i]][this_index]['motherID']=trackID
        
        local_cells=track_cells_weight_r(this_index, from_frame, to_frame, local_cells,weight_vec)
    
    return local_cells

def track_cells_weight_r( tracked_cell_index, from_frame, to_frame, local_cells, weight_vec):
    
    for this_frame in range(from_frame,to_frame,-1):
        this_cell=local_cells[this_frame][tracked_cell_index]
        to_track_cells=local_cells[this_frame-1]

        this_cell_trackID=local_cells[this_frame][tracked_cell_index]['trackID']
        this_cell_ID=local_cells[this_frame][tracked_cell_index]['cellID']
        this_cell_color=local_cells[this_frame][tracked_cell_index]['cellColor']
        
        [tracked_cell_index, tracked_info]=track_cell_weight(this_cell, to_track_cells,weight_vec)
        tscoreV,tscore=tracked_info
        tscore=float(tscore)
        if(tracked_cell_index==-1):
            print("** Cell %s not tracked",this_cell['cellID'])
            return local_cells
        
        this_tracked_trackID=local_cells[this_frame-1][tracked_cell_index]['trackID']
        
        if(len(this_tracked_trackID)>=1):
            #print('.',end='')
            local_cells[this_frame-1][tracked_cell_index]['motherID']=this_cell_trackID
            local_cells[this_frame-1][tracked_cell_index]['trackedBy_next_frame'].append(this_cell_ID)
            temp_list=local_cells[this_frame-1][tracked_cell_index]['trackedBy_next_frame']
            temp_set=list(set(temp_list))
            temp_indexs=[temp_list.index(x) for x in temp_set]
            local_cells[this_frame-1][tracked_cell_index]['trackedBy_next_frame']=temp_list
            local_cells[this_frame-1][tracked_cell_index]['tracking_score'].append(tscore)
            temp_list=local_cells[this_frame-1][tracked_cell_index]['tracking_score']
            new_list=[temp_list[x] for x in temp_indexs]
            local_cells[this_frame-1][tracked_cell_index]['tracking_score']=new_list
            #local_cells[this_frame-1][tracked_cell_index]['tracking_score']=list(set(local_cells[this_frame-1][tracked_cell_index]['tracking_score']))
            
            return local_cells
        
        else:
            local_cells[this_frame-1][tracked_cell_index]['trackID'].extend(this_cell_trackID)
            local_cells[this_frame-1][tracked_cell_index]['trackID']=list(set(local_cells[this_frame-1][tracked_cell_index]['trackID']))
            local_cells[this_frame-1][tracked_cell_index]['cellColor']=this_cell_color
    
        local_cells[this_frame-1][tracked_cell_index]['motherID']=this_cell_trackID
        local_cells[this_frame-1][tracked_cell_index]['trackedBy_next_frame'].append(this_cell_ID)
        temp_list=local_cells[this_frame-1][tracked_cell_index]['trackedBy_next_frame']
        temp_set=list(set(temp_list))
        temp_indexs=[temp_list.index(x) for x in temp_set]
        local_cells[this_frame-1][tracked_cell_index]['trackedBy_next_frame']=temp_list
        local_cells[this_frame-1][tracked_cell_index]['tracking_score'].append(tscore)
        temp_list=local_cells[this_frame-1][tracked_cell_index]['tracking_score']
        new_list=[temp_list[x] for x in temp_indexs]
        local_cells[this_frame-1][tracked_cell_index]['tracking_score']=new_list
        
    return local_cells


def get_polys_coverage(poly1,poly2):

    a,c=get_max_x_angle(poly1)
    poly1=affinity.rotate(poly1,a,c)
    
    a,c=get_max_x_angle(poly2)
    poly2=affinity.rotate(poly2,a,c)

    (bminx, bminy, bmaxx, bmaxy)=poly1.bounds
    poly1=affinity.translate(poly1,-bminx,-bminy)
    (bminx, bminy, bmaxx, bmaxy)=poly2.bounds
    poly2=affinity.translate(poly2,-bminx,-bminy)

    poly1=poly1.convex_hull
    poly2=poly2.convex_hull
    
    pI=poly1.intersection(poly2)
 
    coverage_1=pI.area/poly1.area
    coverage_2=pI.area/poly2.area
    
    #coverage=max(coverage_1,coverage_2)
    coverage=min(coverage_1,coverage_2)
    
    return coverage

def fix_poly(p):
    buff=0
    while(p.is_valid==False or p.geom_type is not "Polygon"):
        p=p.buffer(buff)
        p=wkt.loads(wkt.dumps(p, rounding_precision=2)).simplify(0)
        #print(buff)
        buff+=.001
    
    return p


def track_cell_weight(this_cell, to_track_cells, weight_vec):
  
    neighbors=[]
    max_neighbors=10
    
    tracked_cell_index=-1
    
    this_cell_center=this_cell['center']
    
    # This cell calculations
    this_cell_axis=this_cell['axis']
    this_cell_box=this_cell_axis.bounds
    this_cell_angle=np.arctan((this_cell_box[3]-this_cell_box[1])/(this_cell_box[2]-this_cell_box[0]))    
    
    this_cell_axis=this_cell['axis']
    this_cell_poly=this_cell['roiPoly']
    this_cell_GFP=this_cell['GFP']
    this_cell_DsRed=this_cell['DsRed']
  
    #this_cell_poly=this_cell_poly.convex_hull    ####### <---- this works but turns elongated cells in bananas
    this_cell_poly=fix_poly(this_cell_poly)
    
    neighbors=get_N_neighbors(this_cell, to_track_cells,max_neighbors)
    
    # Neighbors calculations
    
    t_weight=float("-inf")
    chosen_w=[]
    for this_ng in neighbors:
        
        to_track_cell=to_track_cells[this_ng['index']]
        to_track_poly=to_track_cell['roiPoly']
        
        #to_track_poly=to_track_poly.convex_hull  ####### <---- this was
        to_track_poly=fix_poly(to_track_poly)
        
        to_track_axis=to_track_cell['axis']
        track_cell_box=to_track_axis.bounds
        track_cell_angle=np.arctan((track_cell_box[3]-track_cell_box[1])/(track_cell_box[2]-track_cell_box[0]))
        angle_diff=this_cell_angle-track_cell_angle
        angle_diff=abs(angle_diff)
        this_ng['angleDiff']=angle_diff
        
        
        polyInter=this_cell_poly.intersection(to_track_poly)
        axisFraction=max(polyInter.area/this_cell_poly.area,polyInter.area/to_track_poly.area)
        #axisFraction=min(polyInter.area/this_cell_poly.area,polyInter.area/to_track_poly.area)
        #axisFraction=polyInter.area#/this_cell_poly.area
        this_ng['axisFraction']=axisFraction
        
        coverage=get_polys_coverage(this_cell_poly,to_track_poly)
        
        this_ng['coverage']=coverage
        

        to_track_fl=to_track_cell['GFP']
        flDiff=this_cell_GFP-to_track_fl
        GFPDiff=abs(flDiff)
        this_ng['GFP']=GFPDiff
        to_track_fl=to_track_cell['DsRed']
        flDiff=this_cell_DsRed-to_track_fl
        DsRedDiff=abs(flDiff)
        this_ng['DsRed']=DsRedDiff
        
        this_vec=[this_ng['distance'],angle_diff,coverage,axisFraction,GFPDiff,DsRedDiff]

        this_ng['weightV']=this_vec

        this_weight=np.dot(this_vec,weight_vec)
        this_ng['weight']=this_weight
        
       
        if(this_weight>t_weight):
            t_weight=this_weight
            tracked_cell_index=this_ng['index']
            chosen_ng=this_ng['weightV']
            chosen_ng2=this_ng['weight']
        
    return [tracked_cell_index, (chosen_ng,chosen_ng2)]






######################################################

def clear_this_cell_track(this_cell,local_cells):
    this_cell_id=this_cell['cellID']
    this_cell_id_list=this_cell_id.split('.')
    

    this_cell_frame=int(this_cell_id_list[0])
    this_cell_index=int(this_cell_id_list[1])
    print("Clearing track of...",this_cell_frame,this_cell_index)
    
    for frame in reversed(range(this_cell_frame)):
        for tracked_index,this_tracked_cell in enumerate(local_cells[frame]):
            this_tracked_cell_trackedList=this_tracked_cell['trackedBy_next_frame']
            if(this_cell['cellID'] in this_tracked_cell_trackedList):
                print(this_tracked_cell['cellID'],end=' ')
               # print(" ",local_cells[frame][tracked_index]['trackID'],this_cell['trackID'])
                local_cells[frame][tracked_index]['trackedBy_next_frame'].remove(this_cell['cellID'])
                
                for this_track_id in this_cell['trackID']:
                    if(local_cells[frame][tracked_index]['trackID']):
                       # print(" ",local_cells[frame][tracked_index]['trackID'],this_track_id)
                        local_cells[frame][tracked_index]['trackID'].remove(this_track_id)
                
                
                this_cell=this_tracked_cell
    
    
    return local_cells



def clear_tracks(local_cells):
    for i in range(len(local_cells)):
        for j in range(len(local_cells[i])):
            local_cells[i][j]['trackID']=[]
            local_cells[i][j]['trackedBy_next_frame']=[]
            local_cells[i][j]['trackedBy_previous_frame']=[]
            local_cells[i][j]['motherID']=[]
    
    return local_cells


def get_N_neighbors(this_cell, to_track_cells,N):
    neighbors=[]
    distances=[]
    threshold_distance=0.
    this_cell_center=this_cell['center']
    
    #Compute threshold distance to the N closest cells
    for to_track_cell in to_track_cells:
        to_track_center=to_track_cell['center']
        this_d=this_cell_center.distance(to_track_center)
        distances.append(this_d)
    
    distances.sort()
    #print("*",len(distances))
    if(len(distances)<=N):
        threshold_distance=distances[len(distances)-1]
    else:
        threshold_distance=distances[N-1]
    
    #Find the N closets cells
    for to_track_index,to_track_cell in enumerate(to_track_cells):
        to_track_center=to_track_cell['center']
        this_d=this_cell_center.distance(to_track_center)
        
        to_track_id=to_track_cell['cellID'] 
        
        if(this_d<=threshold_distance and len(neighbors)<N):
            this_ng=new_neighbor(to_track_index,to_track_id,this_d)
            neighbors.append(this_ng)
    
    
    return neighbors
    
    







############### Lineages FUNCTIONS ##################
def make_cellLineages_all(tracked_cells,frame_experiment_start):
    cell_lineages=[]
    tracked_ids_list=[]
    max_frames=len(tracked_cells)-1  

    from_frame=max_frames
    to_frame=0
    last_lineageID=0
    for this_frame in range(from_frame,to_frame,-1):
        
        toMake_cells=[]
        
        for select_index, this_cell in enumerate(tracked_cells[this_frame]):
            
            if len(this_cell['trackID'])==0:
                print("***",this_cell['trackID'] ,"+", this_cell['trackID'][0],"+", this_cell['cellID'])
            if this_cell['trackID'][0]==this_cell['cellID']:
               # print("**",this_cell['trackID'] ,"+", this_cell['trackID'][0],"+", this_cell['cellID'])
                toMake_cells.append(this_cell)
            else:
                0
        print('\nMaking lineages of %s cells from frame %s out of %s'%(len(toMake_cells),this_frame,len(tracked_cells[this_frame])))
        cell_lineages,tracked_ids_list,last_lineageID=make_cellLineages(cell_lineages, tracked_cells, toMake_cells, this_frame,tracked_ids_list,frame_experiment_start,last_lineageID) 
    print(len(tracked_ids_list))
    print('\nLineage complete!')
    
    return cell_lineages


def make_cellLineages(cell_lineages, tracked_cells, toMake_cells, last_frame, tracked_ids_list,frame_experiment_start,last_lineageID):
    
    for this_cell in toMake_cells:
        this_track_id=this_cell['trackID']
        this_cellId=this_cell['cellID']
        if(this_track_id):
            print('\n',this_track_id, end='\t')
                   
            this_roiID=this_cell['roiID']
            this_roiPoly=this_cell['roiPoly']
            this_GFP=this_cell['GFP']
            this_DsRed=this_cell['DsRed']
            this_axis=this_cell['axis']
            this_center=this_cell['center']
            this_mother=this_cell['motherID']
            this_score=this_cell['tracking_score']
            this_lineage=new_cellLineage("lameID")
            status=0                 
            
            this_lineage=add_trackInfo_cellLineage(this_lineage,this_cellId, this_roiID, last_frame+frame_experiment_start, this_roiPoly, this_axis, this_center, this_GFP, this_DsRed, status,this_score, this_mother )
            
            last_mother=this_cellId
            aux_lineageID=0
            #lineageID=len(cell_lineages)
            lineageID=last_lineageID
            frame_flag=False
            
            for frame in reversed(range(last_frame)):
               # print(frame,end='.')
                if(frame_flag):
                    break
                print('.',end='')
                this_tracked_cells=tracked_cells[frame]
                
                for this_tracked_cell in this_tracked_cells:
                    
                    if len(this_tracked_cell['trackedBy_next_frame'])>0:
                        this_tracked_cellID_list=this_tracked_cell['trackedBy_next_frame']
                    else:
                        this_tracked_cellID_list=[]
                        
                    if((len(this_tracked_cell['trackedBy_next_frame'])>1) and (this_cellId in tracked_ids_list)):
                        print('*',end='')
                        for aux_lineage in cell_lineages:
                            aux_lienage_cellIds=aux_lineage['cellIDs']
                            if(this_cellId in aux_lienage_cellIds):
                                aux_lineageID=aux_lineage['lineageID']
                        frame_flag=True
                        break
                    
                    elif(this_cellId in this_tracked_cellID_list):
                    
                        index=this_tracked_cell['trackedBy_next_frame'].index(this_cellId)
                        this_score=this_tracked_cell['tracking_score'][index]
                        tracked_ids_list.append(this_cellId)
                        tracked_ids_list=list(set(tracked_ids_list))
                        
                        this_cellId=this_tracked_cell['cellID']
                        this_roiID=this_tracked_cell['roiID']
                        this_roiPoly=this_tracked_cell['roiPoly']
                        this_axis=this_tracked_cell['axis']
                        this_center=this_tracked_cell['center']
                        this_GFP=this_tracked_cell['GFP']
                        this_DsRed=this_tracked_cell['DsRed']
                        this_mother=this_tracked_cell['motherID']
                        
                        
                        
                        status=0                 
                        this_lineage=add_trackInfo_cellLineage(this_lineage,this_cellId, this_roiID, frame+frame_experiment_start, this_roiPoly, this_axis, this_center, this_GFP, this_DsRed,status,this_score,this_mother )
                        
                        last_mother=this_mother
                        last_mother=this_cellId
                
            if(aux_lineageID):
                print(lineageID,'-',aux_lineageID,end='\r')
                this_lineage['lineageID']=aux_lineageID
            else:
                last_lineageID=last_lineageID+1
                this_lineage['lineageID']=last_lineageID
            this_lineage['motherID']=last_mother

    
            #     ### Fix reverse tracking.....    
            this_lineage['cellIDs']=(this_lineage['cellIDs'][::-1])
            this_lineage['roiID']=(this_lineage['roiID'][::-1])
            this_lineage['roiFrames']=(this_lineage['roiFrames'][::-1])
            this_lineage['roiPolys']=(this_lineage['roiPolys'][::-1])
            this_lineage['axis']=(this_lineage['axis'][::-1])
            this_lineage['center']=(this_lineage['center'][::-1])
            this_lineage['GFP']=(this_lineage['GFP'][::-1])
            this_lineage['DsRed']=(this_lineage['DsRed'][::-1])
            this_lineage['state']=(this_lineage['state'][::-1])
            this_lineage['tracking_score']=(this_lineage['tracking_score'][::-1])
            
            for t_cell_id in this_lineage['cellIDs']:
                t_cell_id_list=t_cell_id.split('.')
                to_div_frame=int(t_cell_id_list[0])
                to_div_index=int(t_cell_id_list[1])
                
                to_div_frame=to_div_frame-frame_experiment_start
                to_div_cell_tracked_trackedByNext=tracked_cells[to_div_frame][to_div_index]['trackedBy_next_frame']
                if(len(to_div_cell_tracked_trackedByNext)>1):
                    #print('.',end='')
                    this_lineage['divisions'].append(1)                    
                else:
                    this_lineage['divisions'].append(0)
              
            this_lineage['trackID']=str("%s-%s"%(this_lineage['cellIDs'][0],this_lineage['cellIDs'][-1]))
            
            
            cell_lineages.append(this_lineage)
            
    return cell_lineages,tracked_ids_list,last_lineageID





############### MISC FUNCTIONS ##################


def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    # Prints the time difference yielded by generator instance TicToc
    tempTimeInterval = next(TicToc)
    if tempBool:
        print( "Elapsed time: %f seconds." %tempTimeInterval )

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)


    
    
    
print("loaded!")


