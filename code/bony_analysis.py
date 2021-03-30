import cdutil
import numpy as np
import MV2 as MV

###############################################################################################
def bony_sorting_part1(w500,binedges):

    A,B,C = w500.shape
    dx=np.diff(binedges)[0]
    # Compute composite:    
    OKwaps=nanarray((A,B,C,2+len(binedges))) # add 2 for the exceedances
    xx=0
    for x in binedges:
        xx+=1
        w500_bin=MV.masked_less(w500,x)
        OKwaps[...,xx]=MV.masked_greater_equal(w500_bin,x+dx)

    # do the first wap bin:
    OKwaps[...,0]=MV.masked_greater_equal(w500,binedges[0])
    # do the last wap bin:
    OKwaps[...,-1]=MV.masked_less(w500,binedges[-1]+dx)

    return OKwaps   # [month,lat,lon,wapbin]

###############################################################################################
def bony_sorting_part2(OKwaps,data,OKland,WTS,binedges):
    
    # this function maps data from [time,lat,lon] to [time,wapbin]
    A,B,C = data.shape 

    DATA = nanarray((A,3+len(binedges))) # add 2 for the exceedances, 1 for land
    CNTS = MV.zeros((A,3+len(binedges)))
    for xx in range(2+len(binedges)):    
        A1 = MV.masked_where(OKwaps[...,xx].mask,WTS)
        A2 = MV.masked_where(data.mask,A1)
        denom = np.ma.sum(np.ma.sum(A2,-1),-1)
        DATA[...,xx] = np.sum(np.sum(data*A2,-1),-1)/denom # bin-avg data is computed where both data and wap are defined
        CNTS[...,xx] = np.sum(np.sum(A1,-1),-1) # fractional area of this bin includes regions where data is undefined
            
    # now do the land-only average:
    xx+=1
    A1 = MV.masked_where(OKland.mask,WTS)
    A2 = MV.masked_where(data.mask,A1)
    denom = np.ma.sum(np.ma.sum(A2,-1),-1)
    DATA[...,xx] = np.sum(np.sum(data*A2,-1),-1)/denom # bin-avg data is computed where both data and wap are defined
    CNTS[...,xx] = np.sum(np.sum(A1,-1),-1) # fractional area of this bin includes regions where data is undefined

    # Ensure that the area matrix has zeros rather than masked points
    CNTS[CNTS.mask]=0
    
    if np.allclose(0.5,np.sum(CNTS,-1))==False:
        print('sum of fractional counts over all wapbins does not equal 0.5 (tropical fraction)')
        moot
        
    # DATA contains area-weighted averages within each bin
    # CNTS contains fractional areas represented by each bin
    # so summing (DATA*CNTS) over all regimes should recover the tropical contribution to the global mean
    v1 = np.sum(DATA*CNTS,1)
    v2a = 0.5*cdutil.averager(data, axis='xy', weights='weighted')
    v2b = np.ma.sum(np.ma.sum(WTS*data,1),1)
    
    if np.allclose(v1,v2a)==False or np.allclose(v1,v2b)==False:
        print('Cannot reconstruct tropical average via summing regimes')
            
    return DATA,CNTS #[time,wapbin]



###########################################################################
def nanarray(vector):
    """
    this generates a masked array with the size given by vector
    example: vector = (90,144,28)
    similar to this=NaN*ones(x,y,z) in matlab
    """

    this=MV.zeros(vector)
    this=MV.masked_where(this==0,this)

    return this

