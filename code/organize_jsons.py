# Load in the dictionary containing feedbacks or error metrics
# Put its elements in a more natural order (e.g., Tropical marine low cloud amount = region/sfc/sec/name)
# Compute NET = LW+SW feedbacks
# Append this dictionary to the existing json containing feedbacks and error metrics

import numpy as np
import glob
import json
from datetime import date 

datadir = '../data/'

def organize_fbk_jsons(new_dict,mo,ripf):

    # Load in the existing file containing pre-computed CMIP6 feedbacks
    file = datadir+'cmip6_amip-p4K_cld_fbks.json'
    f = open(file,'r')
    old_dict = json.load(f)
    f.close()

    names = ['cld_tot','cld_amt','cld_alt','cld_tau','cld_err']  
    old_dict[mo]={} 
    old_dict[mo][ripf]={}
    for region in new_dict['ALL']['LWcld_tot'].keys():
        old_dict[mo][ripf][region]={}
        for sfc in ['all','ocn','lnd','ocn_asc','ocn_dsc']:
            old_dict[mo][ripf][region][sfc]={}
            for sec in ['ALL','HI680','LO680']:
                old_dict[mo][ripf][region][sfc][sec]={}
                for name in names:
                    # place in a natural order (e.g., Tropical marine low cloud amount = region/sfc/sec/name)
                    old_dict[mo][ripf][region][sfc][sec]['LW'+name] = new_dict[sec]['LW'+name][region][sfc]
                    old_dict[mo][ripf][region][sfc][sec]['SW'+name] = new_dict[sec]['SW'+name][region][sfc]
                    old_dict[mo][ripf][region][sfc][sec]['NET'+name] = new_dict[sec]['LW'+name][region][sfc] + new_dict[sec]['SW'+name][region][sfc]
                # end name loop
            # end section loop:
        # end sfc type loop
    # end region loop
    old_dict['metadata'] = {}
    meta = {
    "date_modified" :   str(date.today()),
    "author"        :   "Mark D. Zelinka <zelinka1@llnl.gov>",
    }
    old_dict['metadata'] = meta

    return(old_dict) # now updated to include info from input dictionary 




def organize_err_jsons(new_dict,mo,ripf):

    # Load in the existing file containing pre-computed CMIP6 error metrics
    file = datadir+'cmip6_amip_cld_errs.json'
    f = open(file,'r')
    old_dict = json.load(f)
    f.close()

    names = ['E_TCA','E_ctpt','E_LW','E_SW','E_NET']
    old_dict[mo]={} 
    old_dict[mo][ripf]={}
    for region in new_dict['ALL'].keys():
        old_dict[mo][ripf][region]={}
        for sfc in ['all','ocn','lnd','ocn_asc','ocn_dsc']:
            old_dict[mo][ripf][region][sfc]={}
            for sec in ['ALL','HI680','LO680']:
                old_dict[mo][ripf][region][sfc][sec]={}
                for name in names:
                    # place in a natural order (e.g., Tropical marine low cloud amount = region/sfc/sec/name)
                    try:
                        old_dict[mo][ripf][region][sfc][sec][name] = new_dict[sec][region][sfc][name] 
                    except:
                        old_dict[mo][ripf][region][sfc][sec][name] = np.nan
                # end name loop
            # end section loop:
        # end sfc type loop
    # end region loop   
    old_dict['metadata'] = {}
    meta = {
    "date_modified" :   str(date.today()),
    "author"        :   "Mark D. Zelinka <zelinka1@llnl.gov>",
    }
    old_dict['metadata'] = meta

    return(old_dict) # now updated to include info from input dictionary 
