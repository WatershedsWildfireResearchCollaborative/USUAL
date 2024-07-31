
# Code to discritize and size attributes to stream network

#execfile(r'D:\Box Sync\Wildfires\Watershed_preprocessing\Code\final\stream_discritization_attributes.py')

#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np
#from def_USUAL import*


#-------------User Inputs------------------------
# if not using the gui set your parameters below
usegui=True
arcpy.env.overwriteOutput=True

if usegui==False:
    indir=r'C:\Users\scott\Desktop\ESK_DATA\usualV2'
    infid="tan"
    # input stream network
    strm=os.path.join(indir,infid+"_stream_1km.shp")
    # input flow accumulation raster
    fac=os.path.join(indir,infid+"_fac.tif")
    # input filled dem
    demf=os.path.join(indir,infid+"_demf.tif")
    # length for each line segment
    s_len=500
    slope_thresh=10**-3 # minimum slope
    # output directory
    out_dir=indir
    #base name for output file naming
    out_fid="tan"
    
#--------------End of User Inputs---------------------
if usegui==True:
    # input stream network
    strm=arcpy.GetParameterAsText(0)
    # input flow accumulation raster
    fac=arcpy.GetParameterAsText(1)
    # input filled dem
    demf=arcpy.GetParameterAsText(2)
    # length for each line segment
    s_len=arcpy.GetParameter(3)
    slope_thresh=arcpy.GetParameter(4) # minimum slope
    # output directory
    out_dir=arcpy.GetParameterAsText(5)
    #base name for output file naming
    out_fid=arcpy.GetParameterAsText(6)
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
#-------------Directory Management-------------------------------

''' old code from autonaming stuff
print2(fac)

# seperate directory and file name from input file
out_dir,fid=os.path.split(fac)
size=len(fid) # get number of charcters
out_fid=fid[:size-7]#get base name
print2(out_fid)
fid=fid[:size-4]#get full name without extension
print2(fid)
'''

# Make a temporary directory to store all intermediate outputs
temp_dir=os.path.join(out_dir,'temp')
chk_mk_dir(temp_dir)
out_temp_dir= os.path.join(out_dir,'temp','RiverDisc') # make directory name
chk_mk_dir(out_temp_dir)
print2('All intermediate files will be output to: \n'+out_temp_dir,usegui)
print2('All final outputs will be written to: \n'+out_dir,usegui)

# Set up directory with base file id for outputs
out_tmp_id=out_temp_dir+"\\"+out_fid
out_id=out_dir+"\\"+out_fid

if usegui==False:
    disp_outputs=False
    
if usegui==True:
    disp_outputs=string2boolean(disp_outputs)
    
#-------------Begin Analysis--------------------------------------
#%%-------------- Clean up inputs for analysis-------------
# first check the spatial refererance of all in puts are the same
SRcheckdata=[strm,fac,demf]
checkSpatialRefs(SRcheckdata,usegui)

#.........Get Raster Cellsize............
# cell size x direction
dx_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEX")
dx=float(dx_temp.getOutput(0))

# cell size y direction
dy_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEY")
dy=float(dy_temp.getOutput(0))

# cell diagonal length
cell_diag=float(math.ceil(math.sqrt(dx**2 + dy**2))) 

# get area of a cell
cell_area=dx*dy

#toloerance to to find points to split the line
tol=cell_diag

#..........................splitting up the lines..............................
print2('Discritizing the river network and computing nodes',usegui)
# dissolve the line to only have spits at the confluences
strdiss=out_tmp_id+"_streamdiss.shp"
arcpy.Dissolve_management(strm,strdiss ,"","", "MULTI_PART", "UNSPLIT_LINES")
arcpy.GeneratePointsAlongLines_management(strdiss, out_tmp_id+"_XYbreakraw.shp", "DISTANCE", str(s_len)+" Meters", "", "")


# Get the endpoints of dissolved network
arcpy.management.FeatureVerticesToPoints(strdiss, out_tmp_id+"dissendpoints.shp", "BOTH_ENDS")

# build buffer around the endpoints
arcpy.analysis.Buffer(out_tmp_id+"dissendpoints.shp", out_tmp_id+"_endBuffer.shp", str(s_len/2)+" Meters", "FULL", "ROUND", "NONE", None, "PLANAR")

# erase the generate points near the confluence
arcpy.analysis.Erase(out_tmp_id+"_XYbreakraw.shp", out_tmp_id+"_endBuffer.shp", out_tmp_id+"_XYbreak.shp", None)

# split the polyline at the points and initial split
arcpy.SplitLineAtPoint_management(strdiss, out_tmp_id+"_XYbreak.shp", out_tmp_id+"splittmp.shp", str(tol))

# get segment lengths
arcpy.management.CalculateGeometryAttributes(out_tmp_id+"splittmp.shp", "LENGTH LENGTH", "METERS")

# make a layer and select all length longer than the minimum length
arcpy.MakeFeatureLayer_management(out_tmp_id+"splittmp.shp", out_fid+"_split_lyr")
arcpy.SelectLayerByAttribute_management(out_fid+"_split_lyr", "NEW_SELECTION", '"Length">'+str(s_len+0.1))
# find the midpoint of those line segments
arcpy.management.FeatureVerticesToPoints(out_fid+"_split_lyr", out_tmp_id+"midpoints.shp", "MIDPOINT")
# Clear the selection 
arcpy.SelectLayerByAttribute_management(out_fid+"_split_lyr", "CLEAR_SELECTION")

#merge midpoints and xybreak together and split liens
splitfid=out_tmp_id+"split.shp"
arcpy.management.Merge([out_tmp_id+"midpoints.shp",out_tmp_id+"_XYbreak.shp"],out_tmp_id+"_XYbreakall.shp" )
arcpy.SplitLineAtPoint_management(out_tmp_id+"_streamdiss.shp", out_tmp_id+"_XYbreakall.shp", splitfid, str(dx/2))


arcpy.management.CalculateGeometryAttributes(splitfid, "LENGTH LENGTH", "METERS")

# check if any lines remaining are longer than user defined length
# arr_check = arcpy.da.FeatureClassToNumPyArray( splitfid, ['FID','LENGTH'], skip_nulls=True)
# fid=arr_check['FID']# line id
arr_check = arcpy.da.FeatureClassToNumPyArray( splitfid, ['LENGTH'], skip_nulls=True)
l=arr_check['LENGTH']
l2=l[l>s_len+0.1]


# split those lines again
if len(l2)>0:
    # make a layer and select all length longer than the minimum length
    arcpy.MakeFeatureLayer_management(splitfid, out_fid+"_split_lyr2")
    arcpy.SelectLayerByAttribute_management(out_fid+"_split_lyr2", "NEW_SELECTION", '"Length">'+str(s_len+0.1))
    # find the midpoint of those line segments
    arcpy.management.FeatureVerticesToPoints(out_fid+"_split_lyr2", out_tmp_id+"midpoints.shp", "MIDPOINT")
    # Clear the selection 
    arcpy.SelectLayerByAttribute_management(out_fid+"_split_lyr2", "CLEAR_SELECTION")

    #merge midpoints and xybreak together and split liens
    splitfid=out_tmp_id+"split2.shp"
    arcpy.management.Merge([out_tmp_id+"midpoints.shp",out_tmp_id+"_XYbreakall.shp"],out_tmp_id+"_XYbreakall2.shp" )
    arcpy.SplitLineAtPoint_management(out_tmp_id+"_streamdiss.shp", out_tmp_id+"_XYbreakall2.shp", splitfid, str(dx/2))

#Add GridID to attribute table
arcpy.AddField_management(splitfid, "GridID", "Double")
arcpy.CalculateField_management(splitfid, "GridID", "!FID!+1","PYTHON","")
#print('if order is significant maybe need to Linear Referencing Tools > Create Routes and reorder the data')
#[repeat this to fix any issues of to long]



# Add lengths to the table and delete old data
arcpy.management.CalculateGeometryAttributes(splitfid, "LENGTH LENGTH", "METERS")
arcpy.AddField_management(splitfid, "Length_m", "Double")
arcpy.CalculateField_management(splitfid, "Length_m", "!LENGTH!","PYTHON","")
arcpy.DeleteField_management(splitfid, ["LENGTH","arcid","grid_code","from_node","to_node"])

# Get the starting points of the line segments
arcpy.FeatureVerticesToPoints_management(splitfid,out_tmp_id+"start_points.shp","START")
# get all the end pointss of the line segment
arcpy.FeatureVerticesToPoints_management(splitfid,out_tmp_id+"end_points_all.shp","END")

# get the max GridID value of start points
cursor = arcpy.da.SearchCursor(out_tmp_id+"start_points.shp", "GridID")
FirstRecord = True
for row in cursor:
    if FirstRecord:
        FirstRecord = False
        #MinValue = int(row[0])
        MaxValue = int(row[0])
    else:
        #MinValue = min(int(row[0]),MinValue)
        MaxValue = max(int(row[0]),MaxValue)

# Create a point at the end of the stream
# delete all endpoints that overlap with start points
arcpy.Erase_analysis(out_tmp_id+"end_points_all.shp", out_tmp_id+"start_points.shp", out_tmp_id+"end_point.shp", "2 Meters")



# set grid id to +1 largest grid id of 
arcpy.AddField_management(out_tmp_id+"end_point.shp", "GridID", "Double")
arcpy.CalculateField_management(out_tmp_id+"end_point.shp", "GridID", MaxValue+1,"PYTHON","")

# merge start and end points and clean up the data
arcpy.Merge_management(out_tmp_id+"end_point.shp;"+out_tmp_id+"start_points.shp",out_tmp_id+"_nodes.shp")
arcpy.DeleteField_management(out_tmp_id+"_nodes.shp", ["ORIG_FID"])#,"arcid","grid_code","from_node","to_node"])



### PUT THIS IN AN ERRORC CHECK
# snap pour points 
arcpy.gp.SnapPourPoint_sa(out_tmp_id+"_nodes.shp", fac, out_tmp_id+"_nodes_pp.tif", dx, "GridID")
# convert raster to points
arcpy.RasterToPoint_conversion(out_tmp_id+"_nodes_pp.tif", out_tmp_id+"_nodes_pp.shp", raster_field="Value")

#add flow accumulation and elevation data to the points
arcpy.gp.ExtractMultiValuesToPoints_sa( out_tmp_id+"_nodes_pp.shp", [[fac,"fac1"],[demf,"elev_m1"]], "NONE")

# check if there are any pour points with a drainage area of 0
# there are some odd instances where snapping in nested basins are causing the points to snap incorrectly
arrFACchk = arcpy.da.FeatureClassToNumPyArray( out_tmp_id+"_nodes_pp.shp", ['fac1'], skip_nulls=True)
mfac=min(arrFACchk['fac1'])

# get minimum flow accumulation under the dissolved line network
mfactbl=out_tmp_id+"_fac_tbl.dbf"
ZonalStatisticsAsTable(strdiss,"FID",fac,mfactbl, "DATA", "MINIMUM")
tbl=arcpy.da.TableToNumPyArray(mfactbl,"MIN")
mfacthresh=min(tbl["MIN"])
# mfacthresh=50000
# this is a fix
if mfac<mfacthresh: 
    print2("pour point snapped incorrectly, increasing search radius on problem point(s)",usegui)
    # make a selectable layer
    pplyr="ppnodes"
    arcpy.management.MakeFeatureLayer(out_tmp_id+"_nodes_pp.shp",pplyr)
    
    arcpy.management.SelectLayerByAttribute(pplyr,"NEW_SELECTION", "fac1 < "+str(mfacthresh),None)

    # snap pour points 
    arcpy.gp.SnapPourPoint_sa(pplyr, fac, out_tmp_id+"_nodes_pp2_sub.tif", cell_diag*2, "grid_code")
    
    # convert raster to points
    arcpy.RasterToPoint_conversion(out_tmp_id+"_nodes_pp2_sub.tif", out_tmp_id+"_nodes_pp2_sub.shp", raster_field="Value")

    #add flow accumulation and elevation data to the points
    arcpy.gp.ExtractMultiValuesToPoints_sa( out_tmp_id+"_nodes_pp2_sub.shp", [[fac,"fac1"],[demf,"elev_m1"]], "NONE")
    
    # delete old pour point and merge in the new one(s)
    arcpy.management.DeleteFeatures(pplyr)
    
    arcpy.management.Merge([pplyr,out_tmp_id+"_nodes_pp2_sub.shp"],out_tmp_id+"_nodes_pp2.shp")
    
    ppnodes=out_tmp_id+"_nodes_pp2.shp"
else:# proceed as usual
    ppnodes=out_tmp_id+"_nodes_pp.shp"
    
# Join the two point shapefiles *** Maybe change this to spatial join
arcpy.JoinField_management(out_tmp_id+"_nodes.shp","GridID", ppnodes, "grid_code", "FID;pointid;grid_code;fac1;elev_m1")

# export to final shape file of nodes
arcpy.FeatureClassToFeatureClass_conversion(out_tmp_id+"_nodes.shp",out_temp_dir,out_fid+"_nodes_att.shp")
# clean up the attribute table 
arcpy.DeleteField_management(out_tmp_id+"_nodes_att.shp", ["Length_m","pointid","grid_code"])

#......................computing attributes....................................
print2("Computing attributes",usegui)

# in areas where points are close together the snapping can loose some points here is my current work around to fix this
#arcpy.MakeFeatureLayer_management(out_tmp_id+"_nodes_att.shp", out_fid+"_node_lyr")
arcpy.gp.ExtractMultiValuesToPoints_sa(out_tmp_id+"_nodes_att.shp", [[fac,"fac2"],[demf,"elev_m2"]], "NONE")
arcpy.AddField_management(out_tmp_id+"_nodes_att.shp","fac", "Double")
arcpy.AddField_management(out_tmp_id+"_nodes_att.shp","elev_m", "Double")

calcDA="def calcDA( fac1 , fac2 , fac ):\n    if fac1==0:\n       return(fac2)\n    else:\n        return(fac1)"
arcpy.CalculateField_management(out_tmp_id+"_nodes_att.shp","fac","calcDA( !fac1! , !fac2! , !fac! )", "PYTHON_9.3",calcDA)

calcE="def calcelev( elev_m1 , elev_m2 , elev_m ):\n    if elev_m1==0:\n       return(elev_m2)\n    else:\n        return(elev_m1)"
arcpy.CalculateField_management(out_tmp_id+"_nodes_att.shp","elev_m","calcDA( !elev_m1! , !elev_m2! , !elev_m! )", "PYTHON_9.3",calcE)

arcpy.DeleteField_management(out_tmp_id+"_nodes_att.shp", ["elev_m1","elev_m2","fac1","fac2"])


arcpy.AddField_management(out_tmp_id+"_nodes_att.shp", "usarea_m2", "Double")
arcpy.CalculateField_management(out_tmp_id+"_nodes_att.shp", "usarea_m2","!fac!*"+str(cell_area),"PYTHON","")

# join all the data
arcpy.SpatialJoin_analysis(splitfid, out_tmp_id+"_nodes_att.shp", out_tmp_id+"_split_join.shp", "JOIN_ONE_TO_MANY", "KEEP_ALL",'', "INTERSECT", "", "")

# copyu the file over to a final output location
arcpy.FeatureClassToFeatureClass_conversion(out_tmp_id+"_nodes_att.shp",out_dir,out_fid+"_nodes_network.shp")


# Extract attributes to arrays
arr = arcpy.da.FeatureClassToNumPyArray( out_tmp_id+"_split_join.shp", ['GridID_1','GridID','Length_m','elev_m','usarea_m2'], skip_nulls=True)
gid=arr['GridID']# line id
l=arr['Length_m']
elev=arr['elev_m']
nodeid=arr['GridID_1']#['Join_FID']
usarea=arr['usarea_m2']


# Caluculate slopes and to node--- this explicitly assume drainage area increases downstream rather than Jon's approach
step=np.unique(gid)# get number of network segments
step=step.astype(int)

# allocate memory to write outputs to
l_i=np.zeros(len(step))
dz_i=np.zeros(len(step))
us_elev=np.zeros(len(step))
ds_elev=np.zeros(len(step))
to_link=np.zeros(len(step))
usarea_km2=np.zeros(len(step))
maxda=np.nanmax(usarea)
for i in step:
    tf=gid==i
    #print(i)
    nid_temp=nodeid[tf]
    elev_temp=elev[tf]
    l_temp=l[tf]
    nodeid_temp=nodeid[tf]
    usarea_temp=usarea[tf]
    #print2(i)
    # adding unique to help with resevoir data where values are the same
    us_elev[i-1]=np.unique(elev_temp[elev_temp==max(elev_temp)])
    ds_elev[i-1]=np.unique(elev_temp[elev_temp==min(elev_temp)])
    
    # ## updated version----does not work right
    # to_link[i-1]=np.min(nodeid_temp[usarea_temp==max(usarea_temp)])
    # # the downstream link should have to_link == gid
    # if  to_link[i-1]==i:
    #     # print(i)
    #     to_link[i-1]=0
    
    ## old version
    if np.max(usarea_temp)==maxda:
        # print(meow)
        to_link[i-1]=0
    else:
        to_link[i-1]=np.min(nodeid_temp[usarea_temp==max(usarea_temp)])
    usarea_km2[i-1]=np.unique(min(usarea_temp)*1e-6)
    dz_i[i-1]=us_elev[i-1]-ds_elev[i-1]
    l_i[i-1]=l_temp[0]

# check for to link =0 issues there should only be 1 with value 0
tchk=len(np.where(to_link==0)[0])

if tchk==0: # this error shouldn't arise; I've never created it but I went ahead and built it in
    terrmsg1=['No to link =0 was found. Please report this error.'
              ' \n A temporary potential fix for this is adjusting ' 
              'the location of your pour point further up or downstream.']
    reporterror(terrmsg1,usegui)
    
if tchk>1: # this one can and does occur need to correct the ones that are not right
    tolnk0idx=np.where(to_link==0)[0]# get all indices where tolink==0
    tusa=usarea_km2[tolnk0idx]# get the upstream areas
    tgid=step[tolnk0idx] # get the grid ids
    ctlnkidx=np.where(tusa==np.max(tusa))[0] # find which one has the max
    cgid=tgid[ctlnkidx]# get the correct tolink ==0 grid id
    fgid=np.delete(tgid,ctlnkidx)# grid ids to fix
    
    if len(fgid)>0:
        for k in fgid:
            to_link[step==k]=cgid
            
    if len(fgid)==0:
        warnmsg=['Multiple to link = 0 detected, the current fix these.'
                 '\n You will need to manualy fix this']
        print2(warnmsg,usegui)
    
    


# print(meow)
slope=dz_i/l_i

# remove slopes that are likely unrealistically small
slope[slope<slope_thresh]=slope_thresh

# build 2D array to pull data from
outdata = np.vstack((to_link,usarea_km2,us_elev,ds_elev,slope)).T

#............................Data Ouputing.....................................
print2('Writing all attributes to final output', usegui)
arcpy.FeatureClassToFeatureClass_conversion(splitfid, out_dir, out_fid+"_network.shp")

rivnetwork=out_id+"_network.shp"
# Add fields of interest to the shape file
flds=['ToLink','usarea_km2','uselev_m','dselev_m','Slope']
for fieldname in flds:
    arcpy.AddField_management(rivnetwork, fieldname, "Double")
    
jj=0
for fieldname in flds:
    with arcpy.da.UpdateCursor(rivnetwork,fieldname) as cursor: #Add data
        ii = 0
        for row in cursor:
            row[0]=outdata[ii,jj]
            cursor.updateRow(row)
            ii+=1
    jj+=1

print2("Finished! Final outputs are located: \n"+out_dir,usegui)

    
if disp_outputs and usegui:
    aprx = arcpy.mp.ArcGISProject('current')
    cmap = aprx.listMaps()[0]  #data to be added to first map listed
    cmap.addDataFromPath(rivnetwork)