from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
nlats=96 ; nlons=193
nstep='0024'
trunc='063'
expid='t10'
var='phi'
#plt.figure(figsize=(12,6))
plt.figure(figsize=(6,6))
map = Basemap(projection='stere',width=18000000,height=18000000,lat_0=90.0,lon_0=-90.0,lat_ts=60.0,resolution='l')
#map = Basemap(projection='ortho',lat_0=60,lon_0=-90,resolution='l')           
map.drawcoastlines(linewidth=0.25)
map.drawcountries(linewidth=0.25)
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))            
lon1,lat1,vor,div,phi,u,v,psi,khi,qv=np.loadtxt('../data_out/SSWE_model_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.dat',unpack=True)
if var == 'psi':
	field = psi.reshape((nlats,nlons))
if var == 'phi':
	field = phi.reshape((nlats,nlons))
if var == 'khi':
	field = khi.reshape((nlats,nlons))
if var == 'vor':
	field = vor.reshape((nlats,nlons))
if var == 'div':
	field = div.reshape((nlats,nlons))
if var == 'u':
	field = u.reshape((nlats,nlons))
if var == 'v':
	field = v.reshape((nlats,nlons))
if var == 'qv':
	field = qv.reshape((nlats,nlons))
lon2 = lon1.reshape((nlats,nlons))
lat2 = lat1.reshape((nlats,nlons))
x,y = map(lon2,lat2)
if var == 'phi':
	#cs = map.contour(x,y,field/9.81,20,cmap='jet')
	#cs = map.contour(x,y,field/98.1,levels=np.arange(460,590,10),cmap='winter')
	cs = map.contour(x,y,(field)/9.81+4450,levels=np.arange(9200,10300,50),cmap='winter')#colors='black')
else:
	#cs = map.contourf(x,y,field,levels=np.arange(0.0,1.1,0.1),cmap='jet')
	cs = map.contourf(x,y,field,10,cmap='seismic')
plt.colorbar(shrink=0.5)
plt.title(var+' 500 hpa '+nstep+'h - T'+trunc+' expid='+expid)
plt.savefig('../plots/'+var+'_500_T'+trunc+'_step_'+nstep+'_expid_'+expid+'_NHa.pdf')
#plt.savefig('../plots/'+var+'_500_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.pdf')
plt.show()
