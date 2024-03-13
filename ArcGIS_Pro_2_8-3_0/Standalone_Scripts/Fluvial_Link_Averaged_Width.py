
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
    #input polyline
    infc_pline=r'D:\test\VB3\VB_network.shp'
    #field to be appended to input polyline with width data
    outfield="VB_width"
    # input area of interest polygon (the one used to generate the transects)
    infc_aoi=r'D:\test\VB3\VB.shp'
    #input transects
    infc_tran=r'D:\test\VB3\VBtransects.shp'
    #ouput directory
    out_dir=r"D:\test\VB3"
    #output basename for file naming
    out_fid="VB3"
    # point densification spacing (meters) along transects
    dx=10
    # output cell size for DEM
    dxr=5


#%%-------------------------------GUI based inputs-----------------------------
if usegui==True:
    #input polygon
    infc_pline=arcpy.GetParameterAsText(0)
    #field to be appended to input polyline with width data
    outfield=arcpy.GetParameterAsText(1)
    # input area of interest polygon (the one used to generate the transects)
    infc_aoi=arcpy.GetParameterAsText(2)
    #input transects
    infc_tran=arcpy.GetParameterAsText(3)
    # point densification spacing (meters) along transects
    dx=arcpy.GetParameter(4)    
    #ouput directory
    out_dir=arcpy.GetParameterAsText(5)
    #output basename for file naming
    out_fid=arcpy.GetParameterAsText(6)

    # output cell size for DEM
    dxr=arcpy.GetParameter(7)
    disp_outputs=arcpy.GetParameterAsText(8)
    
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
#%%--------------------Directory managment------------------------------------
# Check if output directory exists-- if not make the directory
chk_mk_dir(out_dir)
print2('Final outputs will be written to: \n'+out_dir,usegui)
print2("Final link averaged data will be appended to "+infc_pline,usegui)
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

#%%------------------------- Begin analysis------------------------------------
# first check the spatial refererance of all in puts are the same
SRcheckdata=[infc_pline,infc_aoi,infc_tran]
checkSpatialRefs(SRcheckdata,usegui)


# check for a Length_m field created in previous code.... if not create it
if not arcpy.ListFields(infc_tran, "Length_m"):
    arcpy.CalculateGeometryAttributes_management(infc_tran,
                                             [["Length_m", "LENGTH"]],"METERS")

#densify transects for point generation
tranden=out_fid+"_densetransect.shp"
#make duplicate line to not modify origional data
arcpy.conversion.FeatureClassToFeatureClass(infc_tran, out_temp_dir, tranden)
tranden=out_temp_dir+"\\"+tranden#add path to file name


# densify line to have sufficent number of point
arcpy.edit.Densify(tranden, "DISTANCE", str(dx))#+" Meters")

#convert line to points
tranpoints=out_tmp_id+"transectpoints.shp"
arcpy.management.FeatureVerticesToPoints(tranden, tranpoints, "ALL")

# get the data projection
SR = arcpy.Describe(tranpoints).spatialReference

#generate tin where "elevation" is the running width
indata=tranpoints+" Length_m Mass_Points <None>;"+infc_aoi+" <None> Hard_Clip <None>"
outTIN=out_tmp_id+"widthTIN"
arcpy.ddd.CreateTin(outTIN, SR, indata, "CONSTRAINED_DELAUNAY")

#convert the tin to a raster
widthRaster=out_id+"_widthraster.tif"
arcpy.ddd.TinRaster(outTIN, widthRaster, "FLOAT", "LINEAR", "CELLSIZE", 1, dxr)

# use zonal statistics to compute average along each link
outtbl=out_tmp_id+"_width"
arcpy.sa.ZonalStatisticsAsTable(infc_pline, 'FID', widthRaster, outtbl, '', 'MEAN')
arcpy.JoinField_management(infc_pline, 'FID', outtbl, 'FID','MEAN')
arcpy.AddField_management(infc_pline, outfield, "DOUBLE")
arcpy.management.CalculateField(infc_pline,outfield , '!MEAN!', 'PYTHON_9.3')
arcpy.management.DeleteField(infc_pline,'MEAN')

print2("Finished! Final outputs are located: \n"+out_dir,usegui)
print2("Final link averaged data has been appended to: \n"+infc_pline,usegui)

    
if disp_outputs and usegui:
    aprx = arcpy.mp.ArcGISProject('current')
    cmap = aprx.listMaps()[0]  #data to be added to first map listed
    cmap.addDataFromPath(widthRaster)

