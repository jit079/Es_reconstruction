import h5py
import numpy as np
import matplotlib.pyplot as plt
from luts import *
from netCDF4 import Dataset

# read Look-up Table for Tg
def read_Tg_lut():
    
    hf=h5py.File('../auxdata/LUT_Tg.h5')
    #print(list(hf.keys()))
    
    Tg=hf['Tg'][()]

    sz=hf['axis/sz'][()]
    u_o3=hf['axis/u_o3'][()]
    u_o2=hf['axis/u_o2'][()]
    u_h2o=hf['axis/u_h2o'][()]
    wl=hf['axis/wl'][()]

    hf.close()
    
    lut= LUT(Tg, axes=[sz,u_o3,u_o2,u_h2o,wl], names=['solar zenith (deg)','ozone (Dobson)','oxygen (mb)','water vapor (g/cm2)','wavelength (nm)'])
    
    return lut

def read_linear_coeff():
    # linear coefficients, 1st column - intercept, 2nd column - 412 nm, 3rd column - 489 nm
    # 4th column - 555 nm, 4th column - 705 nm
    coeff=np.genfromtxt('../auxdata/coeff.txt',skip_header=0)
    #print(coeff.shape)
    
    return coeff

# extraterrestrial solar irradiance, unit W/m2/nm
def read_Es_toa():
    
    data=np.genfromtxt('../auxdata/E0.txt',skip_header=0)
    wl=data[:,0]
    Es_toa=data[:,1]
    
    return wl,Es_toa

def nn_index(lat,lon,lat0,lon0):
    lat_index_max=np.searchsorted(lat,lat0)
    lon_index_max=np.searchsorted(lon,lon0)

    if lat_index_max>lat.size-1: lat_index_max=lat.size-1
    if lon_index_max>lon.size-1: lon_index_max=lon.size-1

    lat_index_min=lat_index_max-1
    lon_index_min=lon_index_max-1

    if lat_index_min<0: lat_index_min=0
    if lon_index_min<0: lon_index_min=0
    
    lat_index_closest=np.where(np.abs(lat[lat_index_max]-lat0)>np.abs(lat[lat_index_min]-lat0),lat_index_min,lat_index_max)
    lon_index_closest=np.where(np.abs(lon[lon_index_max]-lon0)>np.abs(lon[lon_index_min]-lon0),lon_index_min,lon_index_max)

    return lat_index_closest,lon_index_closest

# t: observation time corresponding to the multispectral measurement (GMT, h)
# lat0: latitude corresponding to the multispectral measurement (degree)
# lon0: longitude corresponding to the multispectral measurement (degree)
# merra2 data corresponding to the day of the multispectral measurement

def get_MERRA2(merra2_file,t,lat0,lon0):

    root=Dataset(merra2_file)
    
    # MERRA2 grid
    lat=root['lat'][:]
    lon=root['lon'][:]
    
    # find the indices correponding to the nearest neighbor
    x_idx,y_idx=nn_index(lat,lon,lat0,lon0)
               
    params=['TO3','TQV','PS']
    
    res=np.zeros(len(params),dtype='float')+np.nan
    for i in range(len(params)):
        val=root.variables[params[i]][:,:,:]
        val=val[:,x_idx,y_idx]

        t1=np.floor(t).astype('int')
        t2=t1+1
        
        if t2==24: t2=t1
        
        #print(t1,t2)

        val1=val[t1]
        val2=val[t2]
        
        #print(val1,val2)

        if params[i]=='TQV':
            val1=val1/10.0  # in the unit of g/cm2
            val2=val2/10.0
        
        if params[i]=='PS':
            val1=val1/100.0 # in the unit of mb 
            val2=val2/100.0


        slope=(val2-val1)/1.
        dist=t-t1
        
        res[i]=val1+slope*dist
   
        
    return res

    

# input_Ed: multispectral Ed at 412, 489, 555, and 705 nm (W/m2/nm)
# output_wl:  the hyperpectral wavelengths at which Ed will be reconstructed
# sz: solar zenith angle (degree)
# ozone: ozone amount (Dobson)
# pressure: sea surface pressure (mb)
# water: water vapor (g/cm2)

def reconstruct_Es(sz,ozone,pressure,water,input_Ed,output_wl):
    
    input_wl=np.array([412,489,555,705])
    
    Tg_lut=read_Tg_lut()
    
    coeff=read_linear_coeff()
    
    wl,Es_toa=read_Es_toa()

    # ozone range: [250,450], pressure range: [1000,1025], water vapor range: [0.1,7]
    # if inputs are outside the defined ranges, interpolation will be performed using the closest values.
    if ozone<250: ozone=250
    if ozone>450: ozone=450
    if pressure<1000: pressure=1000
    if pressure>1025: pressure=1025
    if water<0.1: water=0.1
    if water>7: water=7

    Tg=Tg_lut[Idx(sz+np.zeros(input_wl.size)),Idx(ozone+np.zeros(input_wl.size)),\
                  Idx(pressure+np.zeros(input_wl.size)),Idx(water+np.zeros(input_wl.size)),Idx(input_wl)]

    # normalize input Ed by Tg*Es_toa
    
    Ed0=input_Ed/Tg 
    
    Es_toa_4=np.zeros(4,dtype='float')+np.nan
    
    for i in range(input_wl.size):
        loc=np.where((wl>=input_wl[i]-5) & (wl<=input_wl[i]+5))[0]
        Es_toa_4[i]=np.mean(Es_toa[loc])
    
    Ed_4=Ed0/Es_toa_4

    
    x=np.insert(Ed_4,0,1) # add constant
    
    Ed_norm_full=np.matmul(coeff,x.T) # normalized Ed at 0.5 nm resolution (full resolution)
    
    Tg_full=Tg_lut[Idx(sz+np.zeros(wl.size)),Idx(ozone+np.zeros(wl.size)),\
                Idx(pressure+np.zeros(wl.size)),Idx(water+np.zeros(wl.size)),Idx(wl)]
    
    Ed_full=Ed_norm_full*Tg_full*Es_toa


    Ed_output=np.zeros(output_wl.size,dtype='float')+np.nan
    
    for i in range(len(output_wl)):
        loc=np.where((wl>=output_wl[i]-5) & (wl<=output_wl[i]+5))[0]
        Ed_output[i]=np.mean(Ed_full[loc])

    return Ed_output

