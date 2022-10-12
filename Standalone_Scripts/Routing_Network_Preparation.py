
#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
#from def_USUAL import*


#-------------User Inputs------------------------
# if not using the gui set your parameters below
usegui=True
arcpy.env.overwriteOutput=True

if usegui==False:
    # input stream network
    in_strm=r"D:\test\a000\a000_network.shp"
    #output directory
    out_dir=r"D:\test"
    
    snapnodes=True 
    # input nodes of interest in a list
    inputnodes=[r"D:\test\a000\a000subbasins_pp.shp",
                r"D:\test\a000\a000_chillslope_nodes.shp",
                r"D:\test\a000\a000fhillslopes_pp.shp"]
    maxsnapdist=5000

        
    flagfeatures=True    
    in_aoi=r"D:\test\a000\a000excludedarea.shp"
    out_field="Res"

    
if usegui==True:
    # input stream network
    in_strm=arcpy.GetParameterAsText(0)

    #output directory
    out_dir=arcpy.GetParameterAsText(1)
    
    snapnodes=arcpy.GetParameterAsText(2)
    # input nodes of interest in a list
    inputnodes=arcpy.GetParameterAsText(3)
    
    maxsnapdist=arcpy.GetParameterAsText(4)

 
    flagfeatures=arcpy.GetParameterAsText(5)    
    in_aoi=arcpy.GetParameterAsText(6)
    out_field=arcpy.GetParameterAsText(7)

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

    
#%%---------------Correct booleans for arcgis boolean string------------------
if usegui:
    snapnodes=string2boolean(snapnodes)
    flagfeatures=string2boolean(flagfeatures)
    if inputnodes: # convert the arc interface sting to list
        inputnodes=inputnodes.split(";") 
    
#%% Check that all neccessary inputs are provided
if snapnodes: 
    if not inputnodes:
        errmsg="No input nodes provided"
        reporterror(errmsg,usegui)   
    if not out_dir:
        errmsg="No output directory provided"
        reporterror(errmsg,usegui)
    if not inputnodes:
        errmsg="No input nodes provided"
        reporterror(errmsg,usegui)   


if flagfeatures:
    if not in_aoi:
        errmsg="No input polygon provided provided"
        reporterror(errmsg,usegui)



#%%-------------Directory Management-------------------------------
# Check if output directory exists-- if not make the directory
chk_mk_dir(out_dir)
print2('Final outputs will be written to: \n'+out_dir,usegui)


# Make a temporary directory to store all intermediate outputs
tempdir=os.path.join(out_dir,'temp')
chk_mk_dir(tempdir)
out_temp_dir= os.path.join(out_dir,'temp','DataCleaning') # make directory name
chk_mk_dir(out_temp_dir)
print2('Temp outputs will be written to: \n'+out_temp_dir,usegui)
    
#%% -------------Flag the Resevoirs----------------------------------
if flagfeatures:
    # Add a field and flag all resevoir data with one
    arcpy.AddField_management(in_aoi,out_field, "Double")
    arcpy.CalculateField_management(in_aoi, out_field, "1", "PYTHON")
    
    out_feature_class=out_temp_dir+"networkres.shp"
    inputdat=[in_strm,in_aoi]
    arcpy.Intersect_analysis(inputdat, out_feature_class, "ALL", "", "INPUT")
    
    arcpy.JoinField_management(in_strm, "GridID", out_feature_class, 
                               "GridID", out_field)

#%%-------------------Node Snapping and Attribution--------------------------
if snapnodes:
    f2disp=[]
    for infc in inputnodes:
        #get just file name
        [d,fid]=os.path.split(infc)
        #remove extension
        fid=fid.split('.')
        fid=fid[0]
        print2(fid,usegui)
        # make copy to edit
        newfid=fid+"_snap.shp"
        print2(infc,usegui)
        print2(out_temp_dir,usegui)
        print2(newfid,usegui)
        arcpy.FeatureClassToFeatureClass_conversion(infc, out_temp_dir,newfid)
        newfid=out_temp_dir+"\\"+newfid
        #snap the points to the edge of the stream network
        arcpy.Snap_edit(newfid,[[in_strm,"EDGE",str(maxsnapdist)+" Meters"]])
        outfid=out_dir+'\\'+fid+"_snap_att.shp"
        #Add stream attributes to nodes
        arcpy.SpatialJoin_analysis(newfid, in_strm, outfid,"JOIN_ONE_TO_ONE", 
                                   "KEEP_ALL",'',"CLOSEST", str(10), "")
        f2disp.append(outfid) # build file for plotting
       