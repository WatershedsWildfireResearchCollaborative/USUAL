
import numpy as np
from scipy import spatial
import time
import arcpy
import os
#from def_USUAL import*


#-------------User Options----------------------------------------
arcpy.env.overwriteOutput=True
usegui= True

#%%-------------User Inputs----------------------------------------
# if not using the gui set your parameters below
if usegui==False:
    #input polygon
    infcpol=r'D:\test\VB3\VB.shp'
    #input centerline
    infccnt=r'D:\test\VB3\centerline.shp'
    # point spacing to generate transects
    tranden=10
    #polygon that overlaps input polygon where points will not be generated
    exludeaoi=r'D:\test\VB3\VB_RM.shp'#r'D:\test\area2exlclude.shp' #optional set to [] if not using
    # removes all points with x distance of intersection of centerline and polygon
    rmprad=25#optional set to [] if not using
    #output directory
    out_dir=r"D:\test\VB3"
    #basename for outputs
    out_fid='VB'

if usegui==True:
    #input polygon
    infcpol=arcpy.GetParameterAsText(0)
    
    #input centerline
    infccnt=arcpy.GetParameterAsText(1)

    # point spacing to generate transects
    tranden=arcpy.GetParameter(2)
    
    #output directory
    out_dir=arcpy.GetParameterAsText(3)
    
    #basename for outputs
    out_fid=arcpy.GetParameterAsText(4)
    
    #polygon that overlaps input polygon 
    #where points will not be used to build transects
    try:
        exludeaoi=arcpy.GetParameterAsText(5)#optional set to [] if not using
    except:
        exludeaoi=[]
    
    # search radius to exlude points near polygon and centline interesections
    try:
        rmprad=arcpy.GetParameter(6)#optional set to [] if not using
    except:
        rmprad=[]
        

    disp_outputs=arcpy.GetParameterAsText(7)

#%%------------------------Functions-------------------------------------------
# toggle between arcgis interface and command prompt printing
def print2(instr,usegui): 
    if usegui==True:
        arcpy.AddMessage(instr)
    else:
        print(instr)

# Check if a directory exists if not make it    
def chk_mk_dir(input_dir):
    try:
        os.mkdir(input_dir)
    except:
        pass
#report errors depending on commandline or gui version
def reporterror(errmsg,usegui):
    if usegui==True:
        arcpy.AddError(errmsg)
        raise arcpy.ExecuteError
    else:
        raise ValueError(errmsg)

#this corrects for arcpy.getparameterastext
def string2boolean(x):
        if x == 'true':
             x = True
        else:
             x = False
        return x

# check that spatial references are the same
def checkSpatialRefs(indata,usegui):
    SRall=[]
    for dataset in indata:
        SRtemp=arcpy.Describe(dataset).spatialReference
        SRall.append(SRtemp.Name)
    TFref=len(np.unique(SRall))# if >1 then the references aren't the same
    
    if TFref==1:# if only 1 value all spatial references match
        srTF=True
    else:
        srTF=False
        #print2(SRall)
        if not srTF:# srTF is false not all input data is projected the same
            srdisp=list(zip(SRall,indata))
            SRerr='Inputs do not have the same projection! \nPlease project all the '\
                'data into the same coordinate system.'
            print2(srdisp,usegui)
            reporterror(SRerr,usegui)
            
# function to find closest points
def do_kdtree(combined_x_y_arrays,points):
    mytree = spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes, dist
#%%--------------------Directory managment------------------------------------
# Check if output directory exists-- if not make the directory
chk_mk_dir(out_dir)
print2('Final outputs will be written to: \n'+out_dir,usegui)


# Make a temporary directory to store all intermediate outputs
tempdir=os.path.join(out_dir,'temp')
chk_mk_dir(tempdir)
out_temp_dir= os.path.join(out_dir,'temp','Tansects') # make directory name
chk_mk_dir(out_temp_dir)
print2('Temp outputs will be written to: \n'+out_temp_dir,usegui)

#.......Set up directory with base file id for outputs.......
out_tmp_id=out_temp_dir+"\\"+out_fid
out_id=out_dir+"\\"+out_fid



if not usegui:
    disp_outputs=False

if usegui:
    disp_outputs=string2boolean(disp_outputs)

#%%--------------------begin analysis-----------------------------------------
#%%-------------- Clean up inputs for analysis-------------
# first check the spatial refererance of all in puts are the same
if exludeaoi:
    SRcheckdata=[infcpol,infccnt,exludeaoi]
else:
    SRcheckdata=[infcpol,infccnt]
checkSpatialRefs(SRcheckdata,usegui)

print2('Beginning preprocessing data',usegui)
# preprocess data to get edge points split by polygon side
pfc=[infcpol,infccnt]
outpfc=out_tmp_id+'splitpoly.shp'
# seperate left and right polygon
arcpy.management.FeatureToPolygon(pfc, outpfc, None, "NO_ATTRIBUTES", None)

#apply unique line identifier to clean up data later
arcpy.AddField_management(outpfc, "LID", "LONG")
# write the origional ids to origid field
arcpy.CalculateField_management(outpfc, "LID", "!FID!","PYTHON","")

nodeletefield=["FID","LID","Shape"] # fields not to delete
# clean up all the unnecessary attributes
for f in arcpy.ListFields(outpfc):
    if f.name not in nodeletefield:
        try:
            arcpy.DeleteField_management(outpfc,f.name)
        except:
            print(f.name)



#convert to polygon edges to lines
polyperm=out_tmp_id+'LRpolyedge.shp'
arcpy.management.FeatureToLine(outpfc, polyperm, None, "ATTRIBUTES")

#delete the centerline from the polygon to line
plineedge=out_tmp_id+'edge.shp'
arcpy.analysis.Erase(polyperm, infccnt, plineedge, None)

# dissolve by UID to ensure not breaks along an edge
edgediss=out_tmp_id+'edgediss.shp'
arcpy.management.Dissolve(plineedge,edgediss, "LID", None,
                          "MULTI_PART", "DISSOLVE_LINES")


# genererate points along each edge
den_dist=str(tranden)#mapunits can forace a unit via +" Meters"
edgepoints=out_tmp_id+'points.shp'
arcpy.management.GeneratePointsAlongLines(edgediss, edgepoints, "DISTANCE", 
                                          den_dist, None, None)
#add a point unique id 
#apply unique identifier to clean up data later
arcpy.AddField_management(edgepoints, "UID", "LONG")
# write the origional ids to origid field
arcpy.CalculateField_management(edgepoints, "UID", "!FID!","PYTHON","")


#%%------------------ Optional filtering of points-----------------------------
# remove any points that fall under input polygon(s)
if exludeaoi:
    pclip=out_tmp_id+'points_aoiclip.shp'
    arcpy.analysis.Erase(edgepoints, exludeaoi, pclip, None)
    edgepoints=pclip
    
if rmprad:
    # points to initiate search radius--interection of center and polyine
    rmpl=out_tmp_id+'rmpointloc.shp'
    arcpy.analysis.Intersect(pfc, rmpl, "ALL", None, "POINT")
    
    #build buffer around the points
    rmrad=str(rmprad)#mapunits can forace a unit via +" Meters"
    rmaoi=out_tmp_id+"rmaoi.shp"
    arcpy.analysis.Buffer(rmpl, rmaoi, rmrad, "FULL", "ROUND", "NONE", None, "PLANAR")
    
    # delete the points inside the buffer
    rmbclip=out_tmp_id+'points_buffclip.shp'
    arcpy.analysis.Erase(edgepoints, rmaoi, rmbclip, None)
    edgepoints=rmbclip
    
#%%----------------------------Match up points---------------------------------
print2("Finding point pairs",usegui)
# fields to get from shape file
fields=["SHAPE@X","SHAPE@Y","UID","LID"]
# pull data from shape files and write them to arrays
arr = arcpy.da.FeatureClassToNumPyArray( edgepoints, fields, skip_nulls=True)
xp=arr["SHAPE@X"]
yp=arr["SHAPE@Y"]
uid=arr["UID"]
lid=arr["LID"]

# get unique line identifiers to loop over
stepper=np.unique(lid)

#xy = np.dstack([xp.ravel(),yp.ravel()])[0]
# build initial stack of fid xy data [2D array]
xyorig = np.dstack([xp,yp,lid,uid])[0]
xy=xyorig.copy()# make a copy to edit
c=0
out=[0,0,0]
for i in stepper:
    # xy coorinates
    #idx_q=np.logical_and(xy[:,2]==i,xy[:,4]==-9999)
    #idx_i=np.logical_and(xy[:,2]!=i,xy[:,4]==-9999)
    
    xyq=xy[xy[:,2]==i,0:2]#points along the line
    xyi=xy[xy[:,2]!=i,0:2]#other points
    #xyq=xy[idx_q,0:2]#points along the line
    #xyi=xy[idx_i,0:2]#other points
    xyq_uid=xy[xy[:,2]==i,3]
    xyi_uid=xy[xy[:,2]!=i,3]
    #xyq_uid=xy[idx_q,3]
    #    xyi_uid=xy[idx_i,3]
    [idx,dist] = do_kdtree(xyi,xyq)
    cxy=xyi[idx]
    
    # this is gross find elegant sol'n in future
    for j in np.arange(len(cxy)):
        cxytemp=np.hstack((cxy[j,:],c))
        xytemp=np.hstack((xyq[j,:],c))
        out=np.vstack((out,cxytemp,xytemp))
        c+=1
out=np.delete(out,0,0)



oid=np.arange(np.size(out,0))
X=out[:,0]
Y=out[:,1]
cid=out[:,2]
SR = arcpy.Describe(infcpol).spatialReference


array = np.array([(oid[0], (X[0], Y[0]),cid[0])],
                    np.dtype([('idfield',np.int32),('XY', '<f8', 2),('cid',np.int32)]))

for i in range(1,len(X)):
        array = np.concatenate((array,np.array([(oid[i], (X[i], Y[i]),cid[i])],
                            np.dtype([('idfield',np.int32),('XY', '<f8', 2),('cid',np.int32)]))))


#%%-----------------------------Generate final outputs-------------------------
print2('Building lines and final outputs',usegui)
outfc=out_id+"_transect_nodes.shp"
arcpy.da.NumPyArrayToFeatureClass(array, outfc, ['XY'], SR)

arcpy.env.outputCoordinateSystem = SR

tranfc=out_id+"_transects.shp"
arcpy.management.PointsToLine(outfc,tranfc , "cid", None, "NO_CLOSE")

arcpy.CalculateGeometryAttributes_management(tranfc, [["Length_m", "LENGTH"]],"METERS")

print2('Analysis Complete!',usegui)
print2('Final outputs are located at: \n'+out_dir,usegui)

#%%-----------------------------Display final outputs-------------------------

if disp_outputs and usegui:
    aprx = arcpy.mp.ArcGISProject('current')
    cmap = aprx.listMaps()[0]  #data to be added to first map listed
    cmap.addDataFromPath(tranfc)
    cmap.addDataFromPath(outfc)