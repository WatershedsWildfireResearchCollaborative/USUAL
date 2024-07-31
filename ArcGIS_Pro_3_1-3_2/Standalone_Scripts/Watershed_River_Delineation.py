#-------------------General notes about code --------------

# Watershed and river delineation script --> Step 1 in the workflow
# This code delineates watershed(s) and the rivers in them based on a drainage
# area threshold

# To run in python window in arc
#execfile(r'path to file with filename')
# example: execfile(r'D:\Code\final\Watershed_Delineation.py')

#--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np
import glob

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

#-------------User Options----------------------------------------
arcpy.env.overwriteOutput=True
# toogle to use as a standalone script or in a gui
usegui=True#True#
#-------------User Inputs----------------------------------------
# if not using the gui set your parameters below
if usegui==False:
    indir=r'D:\Box Sync\Wildfires\PLI_watersheds\Input_Data'
    # Input DEM (raster)
    in_dem=os.path.join(indir,'DEMS','PLIDEM.tif')
    
    #sq. km minimum drainage area for river channel network delineation
    channel_network_threshold=5 
    
    # Input Pour Points (watershed outlet points for each basin; points shapefile)
    in_pourpoint =os.path.join(indir,'aoi_pp','PLI_PP.shp')
    # name of field with pour point identifiers
    ppfieldin="UID"
    # pour point id INPUT AS STRING; 
    pp_id='73,75'#if using multiple seperate wtih comma or semicolon ex '0,3,97'
    # tolerance to shift pour point
    pp_tol=100 
    
    # Input area of interest shapefile ( polygon shapefile)
    in_aoi= os.path.join(indir,'aoi_pp','PLI_aoi.shp')
    # name of the field housing aoi of interest identifiers
    aoifield=[]#"UID"
    # aoi id input as a string 
    aoi_id=[]#'34,25' #if using multiple seperate with a comma or semicolon
    
    # Set output directories and filename
    # Directory to house output data
    out_dir=r"D:\test\nester"
    # Base file name (name or id of basin you are working with)
    out_fid='pprtest'
    
    isnested=True
    nestfid='mount,little'#names oor outputs of nested basins 
    #nestfid=[]# set this to empty and it will autoname as n1, n2, n3...


# End of Parameters--below is gui code 
if usegui==True:
    
    # Input DEM (raster)
    in_dem=arcpy.GetParameterAsText(0)
    
    #sq. km minimum drainage area for river channel network delineation
    channel_network_threshold=arcpy.GetParameter(1)
    
    # Input Pour Points (watershed outlet points for each basin; points shapefile)
    in_pourpoint = arcpy.GetParameterAsText(2)
    # name of field with pour point identifiers
    ppfieldin=arcpy.GetParameterAsText(3)
    # pour point id INPUT AS STRING; 
    pp_id=arcpy.GetParameterAsText(4)#if using multiple seperate wtih comma ex '0,3,97'
    # tolerance to shift pour point
    pp_tol=arcpy.GetParameter(5)
    
    # Input area of interest shapefile ( polygon shapefile) [optional]
    try:
        in_aoi= arcpy.GetParameterAsText(6)
    except:
        in_aoi=[]
    # name of the field housing aoi of interest identifiers
    try:
        aoifield=arcpy.GetParameterAsText(7)
    except:
        aoifield=[]
    # aoi id input as a string 
    try:
        aoi_id=arcpy.GetParameterAsText(8) #if using multiple seperate with a comma
    except:
        aoi_id=[]
    
    # Set output directories and filename
    # Directory to house output data
    out_dir=arcpy.GetParameterAsText(9) 
    # Base file name (name or id of basin you are working with)
    out_fid=arcpy.GetParameterAsText(10) 
    
    # nested parameters
    isnested=arcpy.GetParameterAsText(11) # toggle nested analysis on and off
    try:#names of outputs of nested basins 
        nestfid=arcpy.GetParameterAsText(12)
    except:# set this to empty and it will autoname as n1, n2, n3...
        nestfid=[]
    # option to add outputs to map
    disp_outputs=arcpy.GetParameterAsText(13)

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
            
#Calculate area of a cell
def findcellulararea_sqkm(in_ras):
    # function extract max value from a raster converts, it to a polygon and 
    #calculates area in sq km and returns a value
    # this approach makes it flexible for all input units
    
    # load a raster into an variable
    ras=arcpy.Raster(in_ras)
    # get maximum value of  raster
    maxrasval=ras.maximum
    #print2(maxrasval,usegui)
    del(ras)#delete unneeded data for computers with less ram
    # null out all other data than the maximum
    singcellras = arcpy.sa.SetNull(in_ras, 1, '"VALUE"<'+str(maxrasval))
    poly= r"in_memory\poly"
    arcpy.RasterToPolygon_conversion(singcellras,poly, 
         "NO_SIMPLIFY", "Value","MULTIPLE_OUTER_PART")
    del(singcellras)
    # Calculate area to subdivide basins by
    arcpy.AddField_management(poly,"area_km","Double")
    expression1 = "{0}".format("!SHAPE.area@SQUAREKILOMETERS!")        
    arcpy.CalculateField_management(poly, "area_km", expression1, "PYTHON","")
    arr = arcpy.da.FeatureClassToNumPyArray(poly, ['area_km'], skip_nulls=True)
    return arr['area_km'][0]# return the area of the cell

#%%-------------- Clean up inputs for analysis-------------
# first check the spatial refererance of all in puts are the same
if in_aoi:
    SRcheckdata=[in_dem,in_pourpoint,in_aoi]
else:
    SRcheckdata=[in_dem,in_pourpoint]
checkSpatialRefs(SRcheckdata,usegui)

# adjust for how arcgis inputs booleans
if usegui==True:
    isnested=string2boolean(isnested)
    disp_outputs=string2boolean(disp_outputs)

if usegui==False:
    disp_outputs=False

# create an error if aoi is input but other fields are empty
if in_aoi:
    if not aoifield:
        err1='You need to specify a field if using an AOI polygon'
        reporterror(err1,usegui)
    if not aoi_id:
        err1='You need to specify a value in the ID field if using an AOI polygon'
        reporterror(err1,usegui)

# adjust the inputs basin and pour point ids into tuples
if aoi_id:
    aoi_id=aoi_id.replace(';',',')# change semicolon to string if it exists
    aoi_id_temp=[int(x) for x in aoi_id.split(",")]# split string up
    aoi_id=tuple(aoi_id_temp)
    del(aoi_id_temp)

pp_id=pp_id.replace(';',',')# change semicolon to string if it exists
pp_id_temp=[int(x) for x in pp_id.split(",")]
pp_id=tuple(pp_id_temp)
del(pp_id_temp)

#if len(pp_id)>1 and isnested==False:
#    reporterror('The code only uses one watershed if not doing a nested analysis!')

# if using a nested name structure make sure the number of pour points
# is the same as the number of names listed
if isnested:
    if nestfid:
       nest_fid=nestfid.split(",")   
       if len(nest_fid)!=len(pp_id):
           nesterr="If using a naming structure for the nested basins you must"\
               "specify the same number of names as number of pour points"
           reporterror(nesterr,usegui)
           
# check if multiple pour points specified for nested base analysis
if isnested and len(pp_id)<2:
    nesterr="Nested analysis requires at least 2 pour points"
    reporterror(nesterr,usegui)

if isnested==True:
    print2('Doing a nested basin analysis',usegui)    

#%%-------------Directory Managment--------------------------------
# Check if output directory exists-- if not make the directory
chk_mk_dir(out_dir)
print2('Final outputs will be written to: \n'+out_dir,usegui)


# Make a temporary directory to store all intermediate outputs
out_temp_dir= os.path.join(out_dir,'temp') # make directory name
chk_mk_dir(out_temp_dir)

out_temp_dir= os.path.join(out_dir,'temp','wtrshd') # make directory name
chk_mk_dir(out_temp_dir)
print2('Temp outputs will be written to: \n'+out_temp_dir,usegui)

#.......Set up directory with base file id for outputs.......
out_tmp_id=out_temp_dir+"\\"+out_fid
out_id=out_dir+"\\"+out_fid
#%% ------------------------- get input properties ---------------------------
#.........Get Raster Cellsize............
dx_temp=arcpy.GetRasterProperties_management(in_dem, "CELLSIZEX")# cell size x direction
dx=float(dx_temp.getOutput(0))
dy_temp=arcpy.GetRasterProperties_management(in_dem, "CELLSIZEY")# cell size y direction
dy=float(dy_temp.getOutput(0))
#cell_diag=float(math.ceil(math.sqrt(dx**2 + dy**2))) # cell diagonal length


#%% ----------Preparing Data For Watershed delineation-------------------------

print2("Preparing input data",usegui)
demf_tmp=out_tmp_id+"demf_tmp.tif"#temp filled dem 

if in_aoi:
    #........Extract areas of interest.........
    aoi_tmp_fid="_aoi_tmp.shp"#temp aoi shapefile
    aoi_tmp=out_tmp_id+aoi_tmp_fid
    if isnested==True:
        #.........Extract and Disolve nested basins................
        arcpy.FeatureClassToFeatureClass_conversion(in_aoi, out_temp_dir,
                                                    out_fid+"_aoinb_tmp.shp",
                                                    aoifield+" IN"+str(aoi_id))
        # Next step is needed for nested basins
        arcpy.Dissolve_management(out_tmp_id+"_aoinb_tmp.shp", aoi_tmp,"FID",
                                  "", "MULTI_PART", "DISSOLVE_LINES")
    else:
        # get buffered watershed of interest (temp)
        arcpy.FeatureClassToFeatureClass_conversion(in_aoi, out_temp_dir, out_fid+aoi_tmp_fid,
                                                    aoifield+"="+str(aoi_id[0]))
        
    #..........get dem under the shape file.......
    
    tmpDEM=arcpy.sa.ExtractByMask(in_dem, aoi_tmp)
    DEMtmpout=out_tmp_id+"_DEM_aoi_tmp.tif"
    tmpDEM.save(DEMtmpout)
    del(tmpDEM)
    
    #.........Fill in the sinks...........
    
    fillras=arcpy.sa.Fill(DEMtmpout)
    fillras.save(demf_tmp)
    del(fillras)

else:
    #.........Fill in the sinks...........
    fillras=arcpy.sa.Fill(in_dem)
    fillras.save(demf_tmp)
    del(fillras)
    
#.......Compute flow direction (D8 algorithm)........
fdrtmp=out_tmp_id+"_fdr_tmp.tif"
arcpy.gp.FlowDirection_sa(demf_tmp, fdrtmp, "NORMAL")

#.......Compute buffered flow accumulation............
outFac = arcpy.sa.FlowAccumulation(fdrtmp, "", "Float")
facbuff=out_tmp_id+"_fac_tmp.tif"
outFac.save(facbuff)

#..........Adjust and prepare pourpoints........................
# make a temp file and add the origional field values to them for nested extraction
# essentially this corrects if someone uses fid for their input
arcpy.FeatureClassToFeatureClass_conversion(in_pourpoint, out_temp_dir, 
    out_fid+"_ppin.shp")
pporin=out_tmp_id+"_ppin.shp"
ppfield="_orig_ppid"
arcpy.AddField_management(pporin, ppfield, "Double")
arcpy.CalculateField_management(pporin, ppfield, "!"+ppfieldin+"!","PYTHON","")

ppaoi= out_id+"_pp.shp"
if len(pp_id)==1:
    if pp_tol>0:
        print2('Snapping Pour Points',usegui)
        #.........get pour points of interest........
        arcpy.FeatureClassToFeatureClass_conversion(pporin, out_temp_dir, 
            out_fid+"_pporig.shp", ppfield+"="+str(pp_id[0]) )
            
        #........Snap Pour Point to nearest cell.....
        arcpy.gp.SnapPourPoint_sa(out_tmp_id+"_pporig.shp", facbuff, 
            out_tmp_id+"_ppsnap.tif", pp_tol, ppfield)
            
        #.......Convert snap pour to point to point....
        arcpy.RasterToPoint_conversion(out_tmp_id+"_ppsnap.tif", ppaoi,"Value")
    else:
        arcpy.FeatureClassToFeatureClass_conversion(pporin, out_dir, 
            out_fid+"_pp.shp", ppfield+"="+str(pp_id[0]) )

else:  

        
    if pp_tol>0:
        print2('Snapping Pour Points',usegui)
        #.........get pour points of interest........
        arcpy.FeatureClassToFeatureClass_conversion(pporin, out_temp_dir, 
            out_fid+"_pporig.shp", ppfield+" IN "+str(pp_id))
            
        #........Snap Pour Point to nearest cell.....
        arcpy.gp.SnapPourPoint_sa(out_tmp_id+"_pporig.shp", facbuff, 
            out_tmp_id+"_ppsnap.tif", pp_tol, ppfield)
            
        #.......Convert snap pour to point to point....
        arcpy.RasterToPoint_conversion(out_tmp_id+"_ppsnap.tif", ppaoi,"Value")
        arcpy.AddField_management(ppaoi, ppfield, "Double")
        arcpy.CalculateField_management(ppaoi, ppfield, "!grid_code!","PYTHON","")
        
    else:
        arcpy.FeatureClassToFeatureClass_conversion(pporin, out_dir, out_fid+"_pp.shp", 
                                                    ppfield+" IN "+str(pp_id))


#%%---------------------------Watershed Delineation----------------------------
print2("Begining the Watershed Delineation Processes",usegui)
#........Extract the watershed raster..........
wtrshd= out_id+"_watershed.shp"# output watershed
if len(pp_id)==1:
    arcpy.gp.Watershed_sa(fdrtmp,ppaoi, out_tmp_id+"_w_rst.tif", "FID")
    #Convert the watershed raster to polygon
    arcpy.RasterToPolygon_conversion(out_tmp_id+"_w_rst.tif", wtrshd,
                                     "NO_SIMPLIFY", "Value","MULTIPLE_OUTER_PART")
else:
    wtr_rstall=""
    for i in np.arange(len(pp_id)):# using mulitple pour points
        arcpy.FeatureClassToFeatureClass_conversion(ppaoi, out_temp_dir, 
                                                    out_fid+"_pp"+str(i)+".shp", 
                                                    ppfield+" = "+str(pp_id[i]))
        ppi=out_tmp_id+"_pp"+str(i)+".shp"
        tmpwtrshdrst=out_tmp_id+"_w_rst"+str(i)+".tif"
        arcpy.gp.Watershed_sa(fdrtmp,ppi,tmpwtrshdrst, 'FID')#ppfield)
        wtr_rst= out_tmp_id+"_wtrsh"+str(i)+".shp"
        arcpy.RasterToPolygon_conversion(tmpwtrshdrst,wtr_rst, 
             "NO_SIMPLIFY", "Value","MULTIPLE_OUTER_PART")
        #arcpy.AddField_management(wtr_rst, "ppid", "Double")
        #arcpy.CalculateField_management(wtr_rst,"ppid", str(i),"PYTHON","")
        wtr_rstall=wtr_rstall+';'+wtr_rst   
    wtr_rstall=wtr_rstall[1:]
    wtrshdall=out_tmp_id+"_wtrshdall.shp"
    arcpy.analysis.Union(wtr_rstall, wtrshdall, "ONLY_FID", None, "GAPS")
    arcpy.management.Dissolve(wtrshdall,wtrshd, "FID", None, "MULTI_PART", "DISSOLVE_LINES")

#.........Repeat workflow to extract final dem and flow accumulation.........
# get dem under the watershed (final)
dem=out_id+"dem.tif"
arcpy.gp.ExtractByMask_sa(in_dem, out_id+"_watershed.shp", dem)

#.........Fill in the sinks in the new dem...........
# fillras=Fill(dem)
demf=out_id+"_demf.tif"
arcpy.gp.ExtractByMask_sa(demf_tmp, out_id+"_watershed.shp", demf)

# fillras.save(demf)
# del(fillras)
#arcpy.gp.Fill_sa(out_id+"dem.tif", out_id+"demf.tif")

#.......Compute flow direction (D8 algorithm)........
fdr=out_id+"_fdr.tif"
arcpy.gp.ExtractByMask_sa(fdrtmp, out_id+"_watershed.shp", fdr)
# arcpy.gp.FlowDirection_sa(demf,fdr, "NORMAL")

#.......Compute flow accumulation............
# outFlowAccumulation = arcpy.sa.FlowAccumulation(fdr, "", "Float")
fac=out_id+"_fac.tif"
arcpy.gp.ExtractByMask_sa(facbuff, out_id+"_watershed.shp", fac)

# outFlowAccumulation.save(fac)

#.........Convert square km to number of pixels.........
#Calculate area of a cell.
cellarea_km2=findcellulararea_sqkm(fac)

#convert km2 to number of cells
chan_net_thresh=channel_network_threshold/cellarea_km2

# Channel make names based on thresholds
if channel_network_threshold<1:
    chan_net_id='0'+str(int(round(channel_network_threshold*10)))+'km'
else:
    chan_net_id=str(int(round(channel_network_threshold)))+'km'


try:
    #.......create stream network raster.......
    outSetNull = arcpy.sa.SetNull(fac, 1, "VALUE<"+str(chan_net_thresh) )
    rivras=out_tmp_id+"_stream_"+chan_net_id+".tif"
    outSetNull.save(rivras)
except:
    reporterror("No data to set to null. \n The solution is to adjust your pour point location \n"
        "Check that the pour point is located in a cell with a large value in the fac file. \n "
        "Another option is to increase the Pour Point Tolerance",usegui)
    
#.......make polylines of the stream network.......
riv=out_id+"_stream_"+chan_net_id+".shp"
arcpy.gp.StreamToFeature_sa(rivras,fdr, riv, "NO_SIMPLIFY")

#%%------------------------Split out nested Basins----------------------
if isnested==True: 
    # make some layers to select interestection with
    # this helps with getting the naming right in the nested basins
    lyridwtr='wtshdlyr'
    lyridpp='pplyr'
    
    arcpy.MakeFeatureLayer_management(wtrshdall,lyridwtr)
    arcpy.MakeFeatureLayer_management(ppaoi,lyridpp)

    for i in np.arange(len(pp_id)):
        if nestfid: 
            nestdir=os.path.join(out_dir,nest_fid[i])
            chk_mk_dir(nestdir)
            nestout=nestdir+"\\"+nest_fid[i]
        else:
            nestdir=os.path.join(out_dir,'n'+str(i))
            chk_mk_dir(nestdir)
            nestout=nestdir+"\\"+'n'+str(i)        

        # select those portions with 2 cells
        arcpy.SelectLayerByAttribute_management(lyridpp, "NEW_SELECTION", ppfield+'='+str(pp_id[i]))
        arcpy.management.SelectLayerByLocation(lyridwtr, "INTERSECT", lyridpp, None, "NEW_SELECTION", "NOT_INVERT")
        
        
        # glob this for shp and raster and extract by mask and clip
        rasfiles=glob.glob(out_dir+"\*.tif")
        for file in rasfiles:
            indir,fileid=os.path.split(file)         
            outfile=nestout+fileid[len(out_fid):]
            # Execute ExtractByMask
            outExtractByMask = arcpy.sa.ExtractByMask(file, lyridwtr)
            # Save the output 
            outExtractByMask.save(outfile)

        shpfiles=glob.glob(out_dir+"\*.shp")
        for file in shpfiles:
            indir,fileid=os.path.split(file)
            outfile=nestout+fileid[len(out_fid):]
            # Execute Clip
            arcpy.Clip_analysis(file, lyridwtr, outfile)
        # clear the selected watershed
        arcpy.SelectLayerByAttribute_management(lyridwtr, "CLEAR_SELECTION")
        arcpy.SelectLayerByAttribute_management(lyridpp, "CLEAR_SELECTION")

#%%------------------Add ouputs to current arcpro map------------------------
if disp_outputs and usegui:
    outfilelist=["_demf.tif","_fac.tif","_fdr.tif",
             '_watershed.shp',"_stream_"+chan_net_id+".shp",'_pp.shp']
    if not isnested:
        for file in outfilelist:
            aprx = arcpy.mp.ArcGISProject('current')
            cmap = aprx.listMaps()[0]  #data to be added to first map listed
            #print(out_id+file)
            cmap.addDataFromPath(out_id+file)
    if isnested:
        for i in np.arange(len(pp_id)):
            if nestfid: 
                nestdir=os.path.join(out_dir,nest_fid[i])
                nestout=nestdir+"\\"+nest_fid[i]
            else:
                nestdir=os.path.join(out_dir,'n'+str(i))
                nestout=nestdir+"\\"+'n'+str(i)  
                
            for file in outfilelist:
                aprx = arcpy.mp.ArcGISProject('current')
                cmap = aprx.listMaps()[0]  #data to be added to first map listed
                #print(out_id+file)
                dispfile=nestout+file
                print2(dispfile,usegui)
                cmap.addDataFromPath(dispfile)
