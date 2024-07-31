# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:58:19 2022

@author: Scott
"""
import os
import arcpy # needs the spatial analyist toolbox
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#from def_USUAL import*


#--------------Booleans & License Checkout----------------------------
arcpy.env.overwriteOutput=True
arcpy.CheckOutExtension("Spatial")
usegui=True

#--------------User Inputs----------------------------------
if usegui==False:
    infc = r"D:\test\rivers.shp"
    in_field = "usarea_km2"
    a=0.5
    b=0.5
    fittype="powerlaw" # options: powerlaw, exponential, linear
    out_field="Riv_width"

if usegui==True:
    infc= arcpy.GetParameterAsText(0)
    in_field=arcpy.GetParameterAsText(1)
    fittype= arcpy.GetParameterAsText(2)
    a = arcpy.GetParameter(3)
    b = arcpy.GetParameter(4)
    out_field=arcpy.GetParameterAsText(5)
    
    
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
def power_lawfn(x, a, b):
    return a*np.power(x, b)#y=ax**b

def exponetialfn(x, a, b):
    return a*np.exp(b*x)	

def linearfn(x, a, b):
    return a*x+b



#%% Get data of interest


# extract data of interest
arr = arcpy.da.FeatureClassToNumPyArray( infc, in_field, skip_nulls=True)
x=arr[in_field]

if fittype== "powerlaw":
    y=power_lawfn(x, a, b)
    print2("Computing data as power law fit",usegui)

elif fittype == "exponential":
    y=exponetialfn(x, a, b)
    print2("Computing data as expontential fit",usegui)

elif fittype=="linear":
    y=linearfn(x, a, b)
    print2("Generating data as linear fit",usegui)
else:
    errmsg=("You did not specify a valid fit function.\n"
        "Please choose one of the following:\n"
        "powerlaw \nexponential \nlinear")
    reporterror(errmsg,usegui)

arcpy.AddField_management(infc, out_field, "Double")


with arcpy.da.UpdateCursor(infc,out_field) as cursor: #Add data
    i=0    
    for row in cursor:
        row[0]=y[i]
        cursor.updateRow(row)
        i+=1
print2('Analysis completed!',usegui)
