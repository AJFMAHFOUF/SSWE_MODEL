import matplotlib.pyplot as plt
import numpy as np
trunc1='042'
trunc2='042'
step1='0720'
step2='0720'
expid1='tr1'
expid2='tr2'
x1,y1=np.loadtxt('../data_out/SSWE_model_spectrum_T'+trunc1+'_step_'+step1+'_expid_'+expid1+'.dat',unpack=True)
x2,y2=np.loadtxt('../data_out/SSWE_model_spectrum_T'+trunc2+'_step_'+step2+'_expid_'+expid2+'.dat',unpack=True)
plt.loglog(x1,y1,label='step='+step1+'- T'+expid1)
plt.loglog(x2,y2,label='step='+step2+'- T'+expid2)
plt.title('Kinetic energy spectrum')
#plt.xlim(1,1000)
#plt.ylim(10E-9,10E3)
plt.legend()
plt.savefig('../plots/KE_spectrum_'+step1+'_'+expid1+'.png')
plt.show()
