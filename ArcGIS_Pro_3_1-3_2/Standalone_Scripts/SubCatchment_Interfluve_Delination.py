#Subcatchment and Intercfluve delineation --> Step 2 in the workflow
#%%--------------Import Packages----------------------------------
import os
import arcpy # needs the spatial analyist toolbox
from arcpy  import env
from arcpy.sa import *
import math
import numpy as np
from collections import Counter

#from def_USUAL import*
# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")


#%%------------------User Options----------------------------------------
arcpy.env.overwriteOutput=True # allow overwriting
usegui= True # toggle between gui or script
#...................User Inputs..........................................
# if not using the gui set your parameters below
if usegui==False:
    indir=r"D:\PLI\Output_Data\a060N\a070"
    bname="a070"
    
    # Filled DEM raster
    demf=os.path.join(indir,bname+"demf.tif")#r"D:\test\Output_Data\Logandemf.tif"
    #Flow Direction raster
    fdr=os.path.join(indir,bname+"fdr.tif")#r"D:\test\a000\a000fdr.tif"
    # Flow accumualtion raster
    fac=os.path.join(indir,bname+"fac.tif")#r"D:\test\a000\a000fac.tif"
    # raster of channel cells
    #chan_ras=os.path.join(indir,bname+"_stream_5km.tif")#r"D:\test\a000\a000_stream_5km.tif"
    str_network=os.path.join(indir,bname+"_stream_5km.shp")#r"D:\test\a000\a000_stream_5km.shp"
    # watershed shapefile
    wtr_shd=os.path.join(indir,bname+"_watershed.shp")#r"D:\test\a000\a000_watershed.shp"

    # areas to exclude from delineation shapefile (e.g. reservoir polygons)
    excl_shp=r"D:\Box Sync\Wildfires\PLI_watersheds\Input_Data\waterbodies\Reservoirs.shp"
    #excl_shp=[]# use this if you do not want to include an area to exclude
    buffcells=0
    prevent_holes=False
    delinatesub=True
    delinatehill=True
    delinatehillbroad=True
    mincell4coarse=5
    # area threshold (sq. km) for sub-basin delineation
    subbasin_threshold=0.1
    
    # output directory
    out_dir=r'D:\PLI\Output_Data\test_data\finaltest'
    #base name for output file naming
    out_fid="T2"

if usegui==True:
    # End of Parameters--below is gui code 
    demf=arcpy.GetParameterAsText(0)
    fdr=arcpy.GetParameterAsText(1)
    fac=arcpy.GetParameterAsText(2)
    str_network=arcpy.GetParameterAsText(3)
    wtr_shd=arcpy.GetParameterAsText(4)
    subbasin_threshold=arcpy.GetParameter(5)
    out_dir=arcpy.GetParameterAsText(6)
    out_fid=arcpy.GetParameterAsText(7)
    try:
        excl_shp=arcpy.GetParameterAsText(8)
    except:
        excl_shp=[]
    try:
        buffcells=arcpy.GetParameter(9)
    except:
        buffcells=0
    mincell4coarse=arcpy.GetParameter(10)
    prevent_holes=arcpy.GetParameterAsText(11)
    delinatesub=arcpy.GetParameterAsText(12)
    delinatehill=arcpy.GetParameterAsText(13)
    delinatehillbroad=arcpy.GetParameterAsText(14)
    # option to add outputs to map
    disp_outputs=arcpy.GetParameterAsText(15)
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

#%%-------------Directory Managment--------------------------------
# Check if output directory exists-- if not make the directory
chk_mk_dir(out_dir)
print2('Final outputs will be written to: \n'+out_dir,usegui)


# Make a temporary directory to store all intermediate outputs
tempdir=os.path.join(out_dir,'temp')
chk_mk_dir(tempdir)
out_temp_dir= os.path.join(out_dir,'temp','Subbasin') # make directory name
chk_mk_dir(out_temp_dir)
print2('Temp outputs will be written to: \n'+out_temp_dir,usegui)

#.......Set up directory with base file id for outputs.......
out_tmp_id=out_temp_dir+"\\"+out_fid
out_id=out_dir+"\\"+out_fid


if usegui==False: # shuts off loading outputs if not using arcgis gui
    disp_outputs=False

if usegui==True: # convert 'true' or 'false' to True or False 
    delinatesub=string2boolean(delinatesub)
    delinatehill=string2boolean(delinatehill)
    delinatehillbroad=string2boolean(delinatehillbroad)
    prevent_holes=string2boolean(prevent_holes)
    # option to add outputs to map
    disp_outputs=string2boolean(disp_outputs)

#%%-------------Begin Analysis-------------------------------------
# first check the spatial refererance of all in puts are the same
if excl_shp:
    SRcheckdata=[demf,fdr,fac,str_network,wtr_shd,excl_shp]
else:
    SRcheckdata=[demf,fdr,fac,str_network,wtr_shd]
srTF=checkSpatialRefs(SRcheckdata,usegui)

#SR=arcpy.Describe(demf).spatialReference
#units=SR.linearUnitName

#.........Get Raster Cellsize............
# cell size x direction
dx_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEX")
dx=float(dx_temp.getOutput(0))

# cell size y direction
dy_temp=arcpy.GetRasterProperties_management(fac, "CELLSIZEY")
dy=float(dy_temp.getOutput(0))

# cell diagonal length
cell_diag=float(math.ceil(math.sqrt(dx**2 + dy**2))) 

# make sure everything snaps to flow accumulation raster cell
arcpy.env.extent=fdr
arcpy.env.snapRaster=fdr



#%% ----------------Extract Pour points adjacent to rivers---------------------
print2('Extracting Pour Points',usegui)


# build raster mask of river channel
flagfield="oneflagtmp"
arcpy.AddField_management(str_network, flagfield, "DOUBLE")
arcpy.management.CalculateField(str_network, flagfield, "1", "PYTHON3")
chan_ras=out_tmp_id+"_chanras.tif"
arcpy.conversion.PolylineToRaster(str_network, flagfield,chan_ras, 
                                  "MAXIMUM_LENGTH", "NONE", dx, "BUILD")
arcpy.management.DeleteField(str_network, flagfield, "DELETE_FIELDS")


fdist=out_tmp_id+"_fdist.tif"

if excl_shp:
    # do this if a reservoir or other area is being exlcuded
    excl_aoi=out_id+'_excludedarea.shp'
    # Clip polygon to watershed extent
    arcpy.analysis.Clip(excl_shp,wtr_shd,excl_aoi)
    
     
    # buffereing the data to help with odd locations
    if buffcells!=0:
        exbuff=out_tmp_id+'_excludedbuff.shp'
        arcpy.analysis.Buffer(excl_aoi, exbuff, str(dx*buffcells), "FULL", "ROUND", "NONE", None, "PLANAR")
    else:
        exbuff=excl_aoi

    
    try:# if a field of res doesnt exist add it
        arcpy.AddField_management(exbuff,"RES", "LONG")       
    except:
        print2('Field RES already exists it will be replaced with constant values of 1',usegui)       
    arcpy.CalculateField_management(exbuff,"RES", "1", "PYTHON_9.3")
    
    # get the polygon edges to aline with raster cells
    # make raster of no go area
    excras= out_tmp_id+'_excl.tif'
    arcpy.conversion.PolygonToRaster(exbuff,'RES',excras, '', '', str(dx))
    
    if prevent_holes:
        # convert back to polygon now polygon edges should align with rasters
        exbuff1=out_tmp_id+'_excl_buffcell.shp'
        arcpy.conversion.RasterToPolygon(excras, exbuff1, "NO_SIMPLIFY",
                                         "Value", "SINGLE_OUTER_PART", None)     
        # run union to fill holes
        exbuff2=out_tmp_id+'excl_buffnoholes.shp'
        arcpy.analysis.Union([exbuff1,exbuff1], exbuff2, "ALL", None, "NO_GAPS")
        # now dissolve the filled holes into the larger polygon
        excell=out_id+ "_exclude_cell.shp"
        arcpy.management.Dissolve(exbuff2, excell, None, None, 
                                  "SINGLE_PART", "DISSOLVE_LINES", '')
    else:
        # convert back to polygon now polygon edges should align with rasters
        excell=out_id+"_exclude_cell.shp"
        arcpy.conversion.RasterToPolygon(excras, excell, "NO_SIMPLIFY", 
                                         "Value", "SINGLE_OUTER_PART", None)
    
    
    try:# if a field of res doesnt exist add 
        arcpy.AddField_management(excell,"RES", "LONG")
    except:
        print2('Field RES already exists it will be replaced with constant values of 1',usegui)
        
    arcpy.CalculateField_management(excell,"RES", "1", "PYTHON_9.3")
    
    #build a 1 and 0 raster for the watershed (0) and  exclude aoi and river(1)
    excludeAOI=out_tmp_id+"_excl_wtrshd.shp"
    arcpy.Union_analysis([excell,wtr_shd],excludeAOI)
    AOIexcludetif=out_tmp_id+'_AOIexclude.tif'
    arcpy.conversion.PolygonToRaster(excludeAOI, 'RES', AOIexcludetif, '', '', str(dx))
    arcpy.management.MosaicToNewRaster([AOIexcludetif,chan_ras],out_temp_dir,
          out_fid+"_excladj.tif","","",dx,1,"MAXIMUM")
    
    #.......Compute distances to excludedaoi & stream network.......
    exadj=out_tmp_id+"_excladj.tif"
    arcpy.gp.FlowDistance_sa(exadj, demf, fdist, fdr, "HORIZONTAL", "D8", "MINIMUM")

if not excl_shp:
    #.......Compute distances to stream network.......
    arcpy.gp.FlowDistance_sa(chan_ras, demf,fdist ,fdr, "HORIZONTAL", "D8", "MINIMUM")

#.......find pixels in and flowing into to stream network.......
fd_adj=out_tmp_id+"_fd_adj.tif"
arcpy.gp.Con_sa(fdist, "1",fd_adj , "0", "VALUE<="+str(cell_diag) )

#.......get river channel pixels.......
rivchan=out_tmp_id+"_rivchan.tif"
arcpy.gp.Con_sa(fdist, "1", rivchan, "0", "VALUE<="+str(0) )

#.......Get pour points for outlets.......
# get cellts that flow into the river
R=Raster(fd_adj)-Raster(rivchan)
outlets=out_tmp_id+"_outlets.tif"
R.save(outlets)
del R

# null out zero data
ppr=out_tmp_id+"_ppr.tif"
arcpy.gp.SetNull_sa(outlets, "1",ppr , '"Value"=0')

#convert in to point data
abpp=out_tmp_id+"_allbasin_pp.shp" # pour points for all subbasins
arcpy.RasterToPoint_conversion(ppr, abpp, raster_field="Value")

#.......Extract the all sub-catchments and interfluves.......

if excl_shp:
    fdrbasin=out_tmp_id+'_fdrnull.tif'
    arcpy.gp.SetNull_sa(exadj, fdr,fdrbasin, '"Value"=1')    
else:
    fdrbasin=fdr
    
allbasintif=out_tmp_id+"_allbasins.tif"
arcpy.gp.Watershed_sa(fdrbasin, abpp,allbasintif , "pointid")

#.......Convert all sub-catchments and interfluves to vector.......
allbasins1=out_tmp_id+"_all_basinscellular.shp"
arcpy.RasterToPolygon_conversion(allbasintif,allbasins1,"NO_SIMPLIFY", "Value")

allbasins=out_tmp_id+"_all_basins.shp"
arcpy.Dissolve_management(allbasins1,allbasins,"gridcode","",
                          "MULTI_PART","DISSOLVE_LINES")


#%%-------------Seperate sub-catcments and interfluves ---------------------
# Calculate area to subdivide basins by
arcpy.AddField_management(allbasins,"area","Double")
expression1 = "{0}".format("!SHAPE.area@SQUAREKILOMETERS!")        
arcpy.CalculateField_management(allbasins, "area", expression1, "PYTHON", )

#..................Delineate Sub-catchments.....................................
if delinatesub:
    print2("Delineating sub-catchments",usegui)
    # extract sub-catchments based on area threshold
    arcpy.conversion.FeatureClassToFeatureClass(allbasins, out_dir, 
        out_fid+"_subcatchments.shp","area >="+str(subbasin_threshold))
    
    # extract the respective pour points
    subbasinppout=out_id+"_subcatchments_pp.shp"
    subbasinout=out_id+"_subcatchments.shp"
    arcpy.analysis.Clip(abpp,subbasinout, subbasinppout, None)

    #.......Clean up attribute tables..........   
    # Add a new subbasin and origonal id field
    arcpy.AddField_management(subbasinppout, "origid", "LONG")
    arcpy.AddField_management(subbasinppout, "sub_ID", "LONG")
    arcpy.AddField_management(subbasinout, "origid", "LONG")
    
    # write the origional ids to origid field
    arcpy.CalculateField_management(subbasinppout, "origid", "!pointid!","PYTHON","")
    arcpy.CalculateField_management(subbasinout, "origid", "!gridcode!","PYTHON","")
    
    # write new fids to them
    arcpy.CalculateField_management(subbasinppout, "sub_ID", "!FID!","PYTHON","")
    arcpy.management.JoinField(subbasinout, "origid", subbasinppout, "origid", "sub_ID")
    
    #Delete fields that we do not need
    arcpy.DeleteField_management(subbasinppout, "pointid;grid_code;origid")
    arcpy.DeleteField_management(subbasinout, "gridcode;origid")   

#..................Delineate fine resolution hillslopes........................
if delinatehill:
    print2("Delineating fine scale interfluves",usegui)
    
    # extract hillslope aka below threshold
    arcpy.conversion.FeatureClassToFeatureClass(allbasins, out_dir, 
        out_fid+"_f_interfluves.shp","area <"+str(subbasin_threshold))
    
    hillfout=out_id+"_f_interfluves.shp"
    hillfppout=out_id+"_f_interfluves_pp.shp"
    
    arcpy.analysis.Clip(abpp,hillfout, hillfppout, None)
    # think about this section later not sure if we need all of this here
    #.......Clean up attribute tables..........   
    # Add a new subbasin and origonal id field
    arcpy.AddField_management(hillfppout, "origid", "LONG")
    arcpy.AddField_management(hillfppout, "influv_ID", "LONG")
    arcpy.AddField_management(hillfout, "origid", "LONG")
    
    # write the origional ids to origid field
    arcpy.CalculateField_management(hillfppout, "origid", "!pointid!","PYTHON","")
    arcpy.CalculateField_management(hillfout, "origid", "!gridcode!","PYTHON","")
    
    # write new fids to them
    arcpy.CalculateField_management(hillfppout, "influv_ID", "!FID!","PYTHON","")
    arcpy.management.JoinField(hillfout, "origid", hillfppout, "origid", "influv_ID")
    
    #Delete fields that we do not need
    arcpy.DeleteField_management(hillfppout, "pointid;grid_code;origid")
    arcpy.DeleteField_management(hillfout, "gridcode;origid")   
#..................Delineate coarse resolution hillslopes........................
if delinatehillbroad:
    print2("Delineating coarse interfluves",usegui)
    fluvestempid='_fineintefluvestemp.shp'
    fluvestemp=out_tmp_id+fluvestempid
    
    fluvespptempid='_fineintefluvesPPtemp.shp'
    fluvespptemp=out_tmp_id+fluvespptempid

    # we need the fine interfluves copy over output to temp folder or create them
    # incorportate this above
    if delinatehill==True:
        arcpy.conversion.FeatureClassToFeatureClass(hillfout, out_temp_dir, 
                    out_fid+fluvestempid)

        arcpy.conversion.FeatureClassToFeatureClass(hillfppout, out_temp_dir, 
                    out_fid+fluvespptempid)
        
    if delinatehill==False: # turns out we still need this to process the coarse resolution
        # extract hillslope aka below threshold
        arcpy.conversion.FeatureClassToFeatureClass(allbasins, out_temp_dir, 
            out_fid+fluvestempid,"area <"+str(subbasin_threshold))
        
        hillfout=fluvestemp
        hillfppout=fluvespptemp
        
        arcpy.analysis.Clip(abpp,hillfout, hillfppout, None)
        
        #.......Clean up attribute tables..........   
        # Add a new subbasin and origonal id field
        arcpy.AddField_management(hillfppout, "origid", "LONG")
        arcpy.AddField_management(hillfppout, "influv_ID", "LONG")
        arcpy.AddField_management(hillfout, "origid", "LONG")
        
        # write the origional ids to origid field
        arcpy.CalculateField_management(hillfppout, "origid", "!pointid!","PYTHON","")
        arcpy.CalculateField_management(hillfout, "origid", "!gridcode!","PYTHON","")
        
        # write new fids to them
        arcpy.CalculateField_management(hillfppout, "influv_ID", "!FID!","PYTHON","")
        arcpy.management.JoinField(hillfout, "origid", hillfppout, "origid", "influv_ID")
        
        #Delete fields that we do not need
        arcpy.DeleteField_management(hillfppout, "pointid;grid_code")
        arcpy.DeleteField_management(hillfout, "gridcode")   
           
    # convert the river to a polygon
    rivnul=out_tmp_id+"_strm_null.tif"
    arcpy.gp.SetNull_sa(chan_ras, chan_ras, rivnul, '"Value"=0')        
    str_shpnul=out_tmp_id+"_strm_pol.shp"    
    arcpy.RasterToPolygon_conversion(rivnul,str_shpnul, 
                             "NO_SIMPLIFY", "Value","MULTIPLE_OUTER_PART")
    arcpy.DeleteField_management(str_shpnul, "gridcode")
    
    # build a series of buffers to isolate interfluve groups
    # buffer the river to overlap with interfluves
    # large buffer around the stream (2 cell buffer)
    Lbuff_strm=out_tmp_id+"_Lstrmbuff.shp"
    arcpy.analysis.Buffer(str_shpnul, Lbuff_strm,str(2*dx), "FULL", "ROUND", "NONE", None, "PLANAR")
    
    # buffer the exclusion area, erase overlap with large river buffer
    # and merge it with the large river buffer
    if excl_shp:
        #buffer the aoi to exlucde
        aoi2exbuff=out_tmp_id+"_LbuffexAOI.shp"
        arcpy.analysis.Buffer(excell, aoi2exbuff, str(2*dx), "FULL", "ROUND", "NONE", None, "PLANAR")        
        #merge the two buffers together
        mbuffdata=Lbuff_strm+";"+aoi2exbuff
        mbuff1=out_tmp_id+"_buffmerge.shp"
        arcpy.management.Merge(mbuffdata, mbuff1,  None, "NO_SOURCE_INFO")
        mbuff=out_tmp_id+"_buffmergediss.shp"
        arcpy.management.Dissolve(mbuff1,mbuff,None, None, "MULTI_PART", "DISSOLVE_LINES")
        
    if not excl_shp:
       mbuff=Lbuff_strm
      
    # now create a small river buffer to erase the center of the channel
    # this helps with splitting up the areas later (1/10th cell buffer)
    Sbuff_strm=out_tmp_id+"_Sbuffriver.shp"
    arcpy.analysis.Buffer(str_shpnul, Sbuff_strm,str(dx/10), "FULL", "ROUND", "NONE", None, "PLANAR")
    Ebuff=out_tmp_id+"_bufchaninErase.shp"
    arcpy.analysis.Erase(mbuff, Sbuff_strm, Ebuff, None)
    

    # now clip the buffer to the inferfluve extent
    clipbuf=out_tmp_id+"_buffaoidiss.shp"
    arcpy.analysis.Clip(Ebuff, fluvestemp, clipbuf, None)
    
    #split the isolated features into their own part
    singbuff=out_tmp_id+"_buffaoi.shp"
    arcpy.management.MultipartToSinglepart(clipbuf, singbuff)
    
    # Add a unique identifier
    arcpy.AddField_management(singbuff, "Zoneid", "DOUBLE")
    # write the origional ids to origid field
    arcpy.CalculateField_management(singbuff, "Zoneid", "!FID!","PYTHON","")
    arcpy.DeleteField_management(singbuff, "gridcode;BUFF_DIST;RES;ORIG_FID")
    
    
    # above polygon represent continous interfluves
    # now need to split them based on flowing into river segments and 
    #aoi to exlude... This keeps interfluves that bump up to drainage devides 
    # from being lumped together
 
    #First erase river channels in the area 2 exclude if provided
    if excl_shp:
        tribchan=out_tmp_id+"_sr_linexc.shp"   
        arcpy.Erase_analysis(str_network, excell,tribchan,"")
    else:
        tribchan=str_network
    
    # dissolve line to single part to ensure all lines are split at the confluence
    channet=out_tmp_id+"_sr_id.shp"   
    arcpy.management.Dissolve(tribchan, channet, None, None, 
                              "SINGLE_PART", "DISSOLVE_LINES", '')

    # build a line version of the area to exclude
    if excl_shp:
        arcpy.AddField_management(excell, "excID", "LONG")     
        # write a aoi2exlude id field to each polygon
        arcpy.CalculateField_management(excell, "excID", "!FID!","PYTHON","")
        
        # convert the polygon to a line feature
        exaoilinetemp=out_tmp_id+"_exaoiline.shp"   
        arcpy.management.PolygonToLine(excell, exaoilinetemp, "IDENTIFY_NEIGHBORS")
        # need to spatial join the polygon to the line as polygon to line is not preserving attributs
        
        exaoilinejoin=out_tmp_id+"_exaoilinejoin.shp"
        arcpy.analysis.SpatialJoin(exaoilinetemp, excell, exaoilinejoin,
                                   "JOIN_ONE_TO_ONE", "KEEP_ALL", None, 
                                   "INTERSECT", str(dx/2), None)
        
        # NOTE IF THERE ARE HOLES IN THE AOI2EXCLUDE POLYGON USE UNION NO GAPS
        # THEN DISSOLVE AND THIS WILL FIX THIS
        exaoilinediss=out_tmp_id+"_exaoilinediss.shp"
        arcpy.management.Dissolve(exaoilinejoin,exaoilinediss, "excID",
                                  None, "MULTI_PART", "DISSOLVE_LINES", '')

        mLdata=channet+";"+exaoilinediss
        Rivdata2join=out_tmp_id+"_lines2join.shp"
        arcpy.management.Merge(mLdata, Rivdata2join,  None, "NO_SOURCE_INFO")

    else:
        Rivdata2join=channet

    # add unique identifier for each line
    arcpy.AddField_management(Rivdata2join, "Lineid", "LONG")
    arcpy.management.CalculateField(Rivdata2join, "Lineid", "!FID!", 
                                    "PYTHON3", '', "Double", "NO_ENFORCE_DOMAINS")
    
         
     # join the zones of the aoi polygon to the pour points
    ppzone=out_tmp_id+"_ppzoned1.shp"
    arcpy.analysis.SpatialJoin(fluvespptemp, singbuff, ppzone,"JOIN_ONE_TO_ONE", 
                               "KEEP_ALL", None, "CLOSEST", None, '')
    # join the lines of the river/aoi2exlclude polylines to the pour points
    ppzone2=out_tmp_id+"_ppzoned2.shp"
    arcpy.analysis.SpatialJoin(ppzone, Rivdata2join, ppzone2,"JOIN_ONE_TO_ONE", 
                                "KEEP_ALL", None, "CLOSEST", None, '')     
    f2delete="TARGET_FID;JOIN_COU_1;TARGET_F_1;EXCID;Join_Count"
    arcpy.DeleteField_management(ppzone2,f2delete)
  
    # join the pour points with zone id to the interfluve at fine res. 
    fluvzone=out_tmp_id+"_fluvzoned.shp"
    arcpy.analysis.SpatialJoin(fluvestemp, ppzone2, fluvzone,"JOIN_ONE_TO_ONE", 
                                "KEEP_ALL", None, "CLOSEST", None, '')
    
    # dissolve by line ID aka river and aoi2exclude
    Lintfluve_diss=out_tmp_id+"_fluvedissL.shp"
    arcpy.management.Dissolve(fluvzone,Lintfluve_diss, 
                               "Lineid", None, "MULTI_PART", "DISSOLVE_LINES")
     
    # dissolve by zone id groups of interfluves along a river reach
    Zintfluve_diss=out_tmp_id+"_fluvedissZ.shp"
    arcpy.management.Dissolve(fluvzone,Zintfluve_diss, 
                               "Zoneid", None, "MULTI_PART", "DISSOLVE_LINES")
    
    # bring the two dissolves together to build a prelimary c_interlfuce dataset
    cfluveU=out_tmp_id+"_c_union_interfluves.shp"
    arcpy.analysis.Union([Lintfluve_diss,Zintfluve_diss], cfluveU, "ALL", None, "GAPS")
    # this can create some odd coarse interfluves; run one more round of buffer
    # and join
    # buffer prelimenary coarse interfluves to remove square edges
    fluvbuff=out_tmp_id+"_c_ubuff_interfluves.shp"
    arcpy.analysis.Buffer(cfluveU, fluvbuff, str(dx/4), "FULL", "ROUND", 
                          "NONE", None, "PLANAR")
    fluvsing=out_tmp_id+"_c_using_interfluves.shp"
    # now all corner overlap run multipart to singlepart
    arcpy.management.MultipartToSinglepart(fluvbuff,fluvsing)
    # add a unique values
 
    # now clip the buffer to the buffered inferfluve extent
    arcpy.AddField_management(fluvsing, "diss_ID", "Double")
    arcpy.management.CalculateField(fluvsing, "diss_ID", "!FID!", 
                                    "PYTHON3", '', "Double", "NO_ENFORCE_DOMAINS")
    ppzone3=out_tmp_id+"_ppzoned3.shp"
    arcpy.analysis.SpatialJoin(ppzone2, fluvsing, ppzone3,"JOIN_ONE_TO_ONE", 
                                "KEEP_ALL", None, "CLOSEST", None, '')  
    # join the pour points with dissid id to the interfluve at fine res. 
    fluvs=out_tmp_id+"_fluvejoingroups.shp"
    arcpy.analysis.SpatialJoin(fluvestemp, ppzone3, fluvs,"JOIN_ONE_TO_ONE", 
                                "KEEP_ALL", None, "CLOSEST", None, '') 
    
    # dissolve by line ID aka river and aoi2exclude
    cfluve=out_id+"_c_interfluves.shp"
    arcpy.management.Dissolve(fluvs, cfluve, "diss_ID", None, "MULTI_PART",
                              "DISSOLVE_LINES")
    
    # now let's clean up attributes
    arcpy.AddField_management(cfluve, "influv_ID", "Double")
    # converting LineID to influve_ID
    arcpy.management.CalculateField(cfluve, "influv_ID", "!FID!", 
                                    "PYTHON3", '', "Double", "NO_ENFORCE_DOMAINS")
    arcpy.DeleteField_management(cfluve, "diss_ID")
    # the new workflow should fix the missing interfluve issue, however if it 
    # does show back up insert codeblock from previous version here
    # or increase the buffer radius another cell
    
    # get area in map units
    arcpy.management.CalculateGeometryAttributes(cfluve, "area AREA", '', '', 
                                                 None, "SAME_AS_INPUT")
    if mincell4coarse>0:
        arcpy.MakeFeatureLayer_management(cfluve, 'lyrI') #interfluve
        s_exp="area <="+str(mincell4coarse*dx*dy)
        arcpy.management.SelectLayerByAttribute("lyrI", "NEW_SELECTION", s_exp)
        # if features selected delete them
        if int(arcpy.GetCount_management("lyrI")[0]) > 0:
            arcpy.DeleteFeatures_management("lyrI")
            
#%%------------------Add ouputs to current arcpro map------------------------
print2("Finished! Final outputs are located: \n"+out_dir,usegui)

if disp_outputs and usegui:
    outfilelist=[]
    
    if delinatesub==True:
        sbfiles=['_subcatchments.shp','_subcatchments_pp.shp']
        for filename in sbfiles:
            outfilelist.append(filename)
    
    if delinatehill==True:
        fhillfiles=['_f_interfluves.shp','_f_interfluves_pp.shp']
        for filename in fhillfiles:
            outfilelist.append(filename)
    
    if delinatehillbroad==True:
        chillfiles=['_c_interfluves.shp']
        #outfilelist.append(chillfiles)
        for filename in chillfiles:
            outfilelist.append(filename)
    
    
    for file in outfilelist:
        aprx = arcpy.mp.ArcGISProject('current')
        cmap = aprx.listMaps()[0]  #data to be added to first map listed
        #print(out_id+file)
        cmap.addDataFromPath(out_id+file)

