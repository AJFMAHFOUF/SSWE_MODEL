import matplotlib.pyplot as plt
import numpy as np
trunc='063'
step1='0000'
step2='0120'
expid1='t10'
expid2='t10'
plt.figure(figsize=(6,6))
x1,y1=np.loadtxt('../data_out/phi_time_series_t10.dat',unpack=True)
plt.plot(x1*900/3600,y1/9.81 + 4450)
plt.xlim(0,120)
plt.ylim(9700,10150)
plt.title('Geopotential 500 hPa')
plt.xlabel('hours')
plt.ylabel('h (m)')
#plt.legend()
plt.savefig('../plots/timeseries_'+trunc+'_'+expid1+'.png')
plt.show()
