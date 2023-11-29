from reconstruct_Es import *
import numpy as np
import matplotlib.pyplot as plt

Ed=[0.2,0.25,0.23,0.2]

sz=30.

merra2_file='../merra2/MERRA2_401.tavg1_2d_slv_Nx.20200929.SUB.nc'
ozone,water,pressure=get_MERRA2(merra2_file,12,0,0)


print(ozone,water,pressure)

output_wl=np.arange(315,900,5)

Ed_reconstructed=reconstruct_Es(sz,ozone,pressure,water,Ed,output_wl)

np.savetxt('../output/reconstructed_Ed.txt',np.stack((output_wl,Ed_reconstructed)).T)

plt.figure()
plt.plot(output_wl,Ed_reconstructed,'k',label='reconstructed')
plt.scatter([412,489,555,705],Ed,30,'r',label='measured')
plt.xlabel('Wavelength, nm')
plt.ylabel('Ed, W/m2/nm')
plt.legend()
plt.tight_layout()
plt.savefig('../output/reconstructed_Ed.png')
plt.close()
