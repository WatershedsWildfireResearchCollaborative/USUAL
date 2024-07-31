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
#from def_USUAL import*

#--------------Booleans & License Checkout----------------------------
arcpy.env.overwriteOutput=True
arcpy.CheckOutExtension("Spatial")
usegui=True

#--------------User Inputs----------------------------------
if usegui==False:
    infc = r"D:\test\rivers.shp"    
    #input fit type 
    x_field = "usarea_km2"
    y_field="Width_m"
    fittype="powerlaw" # current choices--> powerlaw,exponential,linear

    fidfigsave=r"D:\test\fig.png"
    # any of these extensions will work
    # eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff
    developregression=True
    
if usegui==True:
    infc= arcpy.GetParameterAsText(0)
    x_field = arcpy.GetParameterAsText(1)
    y_field = arcpy.GetParameterAsText(2)
    fittype=arcpy.GetParameterAsText(3)
    try:
        fidfigsave=arcpy.GetParameterAsText(4) 
    except: 
        fidfigsave=[]

    developregression=True

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
# fitting functions
def power_lawfn(x, a, b):
    return a*np.power(x, b)#y=ax**b

def exponetialfn(x, a, b):
    return a*np.exp(b*x)	

def linearfn(x, a, b):
    return a*x+b

#%%-----------------------Directory Management---------------------------------

# get the output directories and make sure they exist
#[d1,f1]=os.path.split(outfc)
[d2,f2]=os.path.split(fidfigsave)

# check if directories exist if not make them
#chk_mk_dir(d1)
chk_mk_dir(d2)

#%% specify the data fitting type

if fittype== "powerlaw":
    infittype=power_lawfn
    print2("Computing fit as power law function",usegui)
    
elif fittype == "exponential":
    infittype=exponetialfn
    print2("Computing fit as expontential law function",usegui)

elif fittype=="linear":
    infittype=linearfn
    print2("Computing fit as linear law function",usegui)
else:
    errmsg=("You did not specify a valid fit function.\n"
        "Please choose one of the following:\n"
        "powerlaw \nexponential \nlinear")
    reporterror(errmsg,usegui)

#%% bein data fitting
#---------------------Extract Data-------------------------
fields=[x_field,y_field]
arr = arcpy.da.FeatureClassToNumPyArray( infc, fields, skip_nulls=True)
xdat=arr[x_field]
ydat=arr[y_field]

# Fit the data
pars, cov = curve_fit(infittype, xdat, ydat, p0=[0, 0], bounds=(-np.inf, np.inf))

# R^2 value
residuals = ydat- infittype(xdat, *pars)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ydat-np.mean(ydat))**2)
r_squared = 1 - (ss_res / ss_tot)



#%%-----------------------Data Plotting---------------------------------------
# plot data and output fits to table
plt.plot(xdat, ydat, 'bo', label='data')


plt.plot(np.sort(xdat), infittype(np.sort(xdat), *pars), 'k--', 
         label='fit: a=%5.3f, b=%5.3f' % tuple(pars))

#plt.annotate("r-squared = {:.3f}".format(r_squared), (0, 1))

plt.xlabel(x_field)
plt.ylabel(y_field)
plt.legend()
if fittype== "powerlaw":
    fittitle="y=%5.3f*x$^{%5.3f}$ " % tuple(pars)

elif fittype == "exponential":
    fittitle="%5.3f e$^{%5.3f*x}$" % tuple(pars)

elif fittype=="linear":
    fittitle="y=%5.3fx+%5.3f" % tuple(pars)
  
plt_title=fittitle + " ; r-squared = {:.3f}".format(r_squared)
plt.title(plt_title)
    
if fidfigsave:# save plot to an output location
    plt.savefig(fidfigsave)
    print2('Plot output to: \n' +fidfigsave,usegui)

plt.show()
    
    
#%%%------------Output coeffients and r2 to the input shapefile----------------
# add output fields
# future considerations dynamic naming set in while loop and solve for name then add 
try:
    arcpy.AddField_management(infc, "coef_a", "DOUBLE")
    arcpy.AddField_management(infc, "coef_b", "DOUBLE")
    arcpy.AddField_management(infc, "r_squared", "DOUBLE")
    
    # write coefficients to the output shapefile
    arcpy.CalculateField_management(infc, "coef_a", str(pars[0]),"PYTHON","")
    arcpy.CalculateField_management(infc, "coef_b", str(pars[1]),"PYTHON","")
    arcpy.CalculateField_management(infc, "r_squared", str(r_squared),"PYTHON","")
except:
    # the code should overwrite but if not this code block will take over
    errmsg=('Input shapefile already contains coef_a, coef_b, and/or r_squared.'
           'Delete or rename these fields and rerun the tool')
    reporterror(errmsg,usegui)
print2('Analysis Complete!',usegui)