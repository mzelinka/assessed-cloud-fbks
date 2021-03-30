import cal_CloudRadKernel as CRK
import organize_jsons as OJ
import cld_fbks_ecs_assessment_v3 as dataviz
import os

# USER INPUT:
#================================================================================================
model = 'YOUR_MODEL'
variant = 'r1i1p1f1'
filenames={}
filenames['amip']={}
filenames['amip-p4K']={}
path = '/p/user_pub/xclim/CMIP6/CMIP/amip/atmos/mon/'
path4K = '/p/user_pub/xclim/CMIP6/CFMIP/amip-p4K/atmos/mon/'

filenames['amip']['tas'] = path+'tas/CMIP6.CMIP.amip.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.tas.atmos.glb-z1-gr1.v20180701.0000000.0.xml'
filenames['amip']['rsdscs'] = path+'rsdscs/CMIP6.CMIP.amip.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.rsdscs.atmos.glb-2d-gr1.v20180701.0000000.0.xml'
filenames['amip']['rsuscs'] = path+'rsuscs/CMIP6.CMIP.amip.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.rsuscs.atmos.glb-2d-gr1.v20180701.0000000.0.xml'
filenames['amip']['wap'] = path+'wap/CMIP6.CMIP.amip.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.wap.atmos.glb-p19-gr1.v20180701.0000000.0.xml'
filenames['amip']['clisccp'] = path+'clisccp/CMIP6.CMIP.amip.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.clisccp.atmos.glb-p7c-gr1.v20180701.0000000.0.xml'

filenames['amip-p4K']['tas'] = path4K+'tas/CMIP6.CFMIP.amip-p4K.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.tas.atmos.glb-z1-gr1.v20180701.0000000.0.xml'
filenames['amip-p4K']['rsdscs'] = path4K+'rsdscs/CMIP6.CFMIP.amip-p4K.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.rsdscs.atmos.glb-2d-gr1.v20180701.0000000.0.xml'
filenames['amip-p4K']['rsuscs'] = path4K+'rsuscs/CMIP6.CFMIP.amip-p4K.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.rsuscs.atmos.glb-2d-gr1.v20180701.0000000.0.xml'
filenames['amip-p4K']['wap'] = path4K+'wap/CMIP6.CFMIP.amip-p4K.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.wap.atmos.glb-p19-gr1.v20180701.0000000.0.xml'
filenames['amip-p4K']['clisccp'] = path4K+'clisccp/CMIP6.CFMIP.amip-p4K.NOAA-GFDL.GFDL-CM4.r1i1p1f1.mon.clisccp.atmos.glb-p7c-gr1.v20180701.0000000.0.xml'
#================================================================================================

    

# calculate all feedback components and Klein et al (2013) error metrics:
fbk_dict,err_dict = CRK.CloudRadKernel(filenames) 

# add this model to the pre-existing json file containing other model results:
updated_err_dict = OJ.organize_err_jsons(err_dict,model,variant) 
updated_fbk_dict = OJ.organize_fbk_jsons(fbk_dict,model,variant)

# plot this model alongside other models and expert assessment:
os.system('mkdir ../figures/')
result = dataviz.make_all_figs(updated_fbk_dict,updated_err_dict,model)

print('Done!')
