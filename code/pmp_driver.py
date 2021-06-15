import cal_CloudRadKernel as CRK
#import organize_jsons as OJ
import cld_fbks_ecs_assessment_v3 as dataviz
import os

import sys
import json

# User Input:
#================================================================================================
model = 'GFDL-CM4'	
institution = 'NOAA-GFDL'
realization = 'r1i1p1f1'
grid_label = 'gr1'
version = 'v20180701'
path = '/p/css03/esgf_publish/CMIP6'
#exp_list = ['amip','amip-p4K']
exp_list = ['amip']
#================================================================================================


# generate xmls pointing to the cmorized netcdf files 
os.system('mkdir ../xmls/')
filenames={}
for exp in exp_list:
    filenames[exp]={}
    if exp=='amip':
        activity = 'CMIP'
    else:
        activity = 'CFMIP'
    for field in ['tas','rsdscs','rsuscs','wap','clisccp']:
        """
        if field=='clisccp':
            table='CFmon'
        else:
            table='Amon'
        searchstring = path+'/'+activity+'/'+institution+'/'+model+'/'+exp+'/'+variant+'/'+table+'/'+field+'/'+grid_label+'/'+version+'/*.nc'
        """
        xmlname = '../xmls/'+exp+'.'+model+'.'+realization+'.'+field+'.'+version+'.xml'
        #os.system('cdscan -x '+xmlname+' '+searchstring)
        filenames[exp][field] = xmlname

print('filenames', json.dumps(filenames, indent=4))
with open('filenames.json', 'w') as f:
    json.dump(filenames, f, sort_keys=True, indent=4)

# calculate all feedback components and Klein et al (2013) error metrics:
fbk_dict, err_dict = CRK.CloudRadKernel(filenames) 

print('calc done')

with open('fbk_dict.json', 'w') as f:
    json.dump(fbk_dict, f, sort_keys=True, indent=4)

with open('err_dict.json', 'w') as f:
    json.dump(err_dict, f, sort_keys=True, indent=4)

# add this model to the pre-existing json file containing other model results:
updated_err_dict = OJ.organize_err_jsons(err_dict, model, realization) 
updated_fbk_dict = OJ.organize_fbk_jsons(fbk_dict, model, realization)

with open('updated_err_dict.json', 'w') as f:
    json.dump(updated_err_dict, f, sort_keys=True, indent=4)

with open('updated_fbk_dict.json', 'w') as f:
    json.dump(updated_fbk_dict, f, sort_keys=True, indent=4)

# plot this model alongside other models and expert assessment:
os.system('mkdir ../figures/')
result = dataviz.make_all_figs(updated_fbk_dict, updated_err_dict, model)

print('Done!')
