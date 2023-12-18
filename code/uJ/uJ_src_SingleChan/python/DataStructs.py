
import random
from shapely.geometry import Point
from shapely.geometry import MultiPoint
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon
from shapely.geometry import box
from shapely import affinity
from shapely.geometry import LineString



print("DataStructs...",end='')


def new_cell(this_cellID):
    this_cell={}
    this_cell['cellID']=this_cellID
    this_cell['cellColor']="#%06x" % random.randint(0, 0xFFFFFF)
    this_cell['roiID']=str()
    this_cell['trackID']=[]  #Just the tracking order change?
    this_cell['roiPoly']=Polygon()
    this_cell['center']=Point()
    this_cell['axis']=LineString()
    this_cell['GFP']=float()
    this_cell['DsRed']=float()
    this_cell['state']=0
    this_cell['trackedBy_previous_frame']=[]
    this_cell['trackedBy_next_frame']=[]
    this_cell['tracking_score']=[]
    this_cell['motherID']=this_cellID
    
    return this_cell

#def add_info_cell(this_cell, this_roiID, this_roiPoly, this_axis, this_center,this_GFP, this_DsRed, this_nGFP, this_nRed):
def add_info_cell(this_cell, this_roiID, this_roiPoly, this_axis, this_center,this_GFP, this_DsRed):
    
    this_cell['roiID']=this_roiID
    this_cell['roiPoly']=this_roiPoly
    this_cell['axis']=this_axis
    this_cell['center']=this_center
    this_cell['GFP']=this_GFP
    this_cell['DsRed']=this_DsRed
    
    return this_cell

def new_neighbor(index,ngID,distance):

    this_neighbor={}
    this_neighbor['ngID']=ngID
    this_neighbor['index']=index
    this_neighbor['distance']=distance
    this_neighbor['angleDiff']=float()
    this_neighbor['coverage']=float()
    this_neighbor['axisFraction']=float()
    this_neighbor['GFP']=float()
    this_neighbor['DsRed']=float()
    this_neighbor['weight']=float()
    this_neighbor['weightV']=[]
    this_neighbor['flDiff']=float()
    
    return this_neighbor 
    

############### LINEAGE FUNCTIONS ##################

def new_cellLineage(this_cellID):
    this_cell={}
    this_cell['trackID']=this_cellID
    this_cell['cellColor']="#%06x" % random.randint(0, 0xFFFFFF)
    this_cell['cellIDs']=[]
    this_cell['roiID']=[]
    this_cell['roiFrames']=[]    
    this_cell['roiPolys']=[]    
    this_cell['axis']=[]
    this_cell['center']=[]
    this_cell['GFP']=[]
    this_cell['DsRed']=[]
    this_cell['state']=[]
    this_cell['tracking_score']=[]
    this_cell['divisions']=[]
    this_cell['motherID']=str()
    this_cell['lineageID']=int()
    
    return this_cell


def add_trackInfo_cellLineage(this_cell,this_cellID, this_roiID, this_frame, this_roiPoly, this_axis, this_center, this_GFP, this_DsRed, status, score, this_motherID ):
    
    this_cell['cellIDs'].append(this_cellID)
    this_cell['roiID'].append(this_roiID)    
    this_cell['roiFrames'].append(this_frame)
    this_cell['roiPolys'].append(this_roiPoly)
    this_cell['axis'].append(this_axis)
    this_cell['center'].append(this_center)    
    this_cell['GFP'].append(this_GFP)
    this_cell['DsRed'].append(this_DsRed)
    this_cell['state'].append(status)
    this_cell['tracking_score'].append(score)
    this_cell['motherID']=this_motherID
    
    return this_cell




print("loaded!")