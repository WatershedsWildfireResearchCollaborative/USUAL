# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 22:06:59 2023

@author: Scott
"""

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
    print2(maxrasval,usegui)
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