#!/usr/bin/env python

import assessed_cloud_fbks.cal_CloudRadKernel as CRK
import assessed_cloud_fbks.compute_ECS as CE
import assessed_cloud_fbks.organize_jsons as OJ
import assessed_cloud_fbks.cld_fbks_ecs_assessment_v3 as dataviz
import os

# User Input:
#================================================================================================
model = 'GFDL-CM4'	
institution = 'NOAA-GFDL'
variant = 'r1i1p1f1'
grid_label = 'gr1'
version = 'v20180701'
path = '/p/css03/esgf_publish/CMIP6'

# Flag to compute ECS
# True: compute ECS using abrupt-4xCO2 run
# False: do not compute, instead rely on ECS value present in the json file (if it exists)
get_ecs = True
#================================================================================================

if get_ecs:
    exps = ['amip','amip-p4K','piControl','abrupt-4xCO2']
else:
    exps = ['amip','amip-p4K']

# generate xmls pointing to the cmorized netcdf files 
os.system('mkdir ../xmls/')
filenames={}
for exp in exps:
    filenames[exp]={}
    if exp=='amip-p4K':
        activity = 'CFMIP'
    else:
        activity = 'CMIP'
    if 'amip' in exp:
        fields = ['tas','rsdscs','rsuscs','wap','clisccp'] # necessary for cloud feedback calcs
    else:
        fields = ['tas', 'rlut', 'rsut', 'rsdt'] # needed for ECS calc
    for field in fields:
        if field=='clisccp':
            table='CFmon'
        else:
            table='Amon'
        searchstring = path+'/'+activity+'/'+institution+'/'+model+'/'+exp+'/'+variant+'/'+table+'/'+field+'/'+grid_label+'/'+version+'/*.nc'
        xmlname = '../xmls/'+exp+'.'+model+'.'+variant+'.'+field+'.'+version+'.xml'
        os.system('cdscan -x '+xmlname+' '+searchstring)
        filenames[exp][field] = xmlname

# calculate all feedback components and Klein et al (2013) error metrics:
fbk_dict,obsc_fbk_dict,err_dict = CRK.CloudRadKernel(filenames) 

# add this model's results to the pre-existing json file containing other models' results: 
updated_fbk_dict,updated_obsc_fbk_dict = OJ.organize_fbk_jsons(fbk_dict,obsc_fbk_dict,model,variant)
updated_err_dict = OJ.organize_err_jsons(err_dict,model,variant)

ecs = None
if get_ecs:
    # calculate ECS and add it to the pre-existing json file containing other models' results:
    ecs = CE.compute_ECS(filenames) 
updated_ecs_dict = OJ.organize_ecs_jsons(ecs,model,variant)

# plot this model alongside other models and expert assessment:
os.system('mkdir ../figures/')
result = dataviz.make_all_figs(updated_fbk_dict,updated_obsc_fbk_dict,updated_err_dict,updated_ecs_dict,model)

print('Done!')
