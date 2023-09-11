import matplotlib.pyplot as plt
import numpy as np
trunc='063'
step1='0000'
step2='0120'
expid1='t10'
expid2='t10'
plt.figure(figsize=(5,5))
x1,y1=np.loadtxt('../data_out/SSWE_model_spectrum_T'+trunc+'_step_'+step1+'_expid_'+expid1+'.dat',unpack=True)
x2,y2=np.loadtxt('../data_out/SSWE_model_spectrum_T'+trunc+'_step_'+step2+'_expid_'+expid2+'.dat',unpack=True)
plt.loglog(x1,y1,label='step='+step1+'-'+expid1)
plt.loglog(x2,y2,label='step='+step2+'-'+expid2)
plt.xlim(1,1000)
plt.ylim(1E-4,1E2)
plt.title('Kinetic energy spectrum')
plt.legend()
plt.savefig('../plots/KE_spectrum2_'+trunc+'_'+expid1+'.png')
plt.show()
