from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
nlats=32 ; nlons=64
nstep='000'
trunc='21'
expid='h02'
var='phi'
plt.figure(figsize=(12,6))
#map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
#            llcrnrlon=0,urcrnrlon=360,resolution='c')
#map = Basemap(projection='ortho',lat_0=40,lon_0=0,resolution='l')           
#map.drawcoastlines(linewidth=0.25)
#map.drawcountries(linewidth=0.25)
#map.drawmeridians(np.arange(0,360,30))
#map.drawparallels(np.arange(-90,90,30))            
lon1,lat1,phi=np.loadtxt('fort.255',unpack=True)
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
lon2 = lon1.reshape((nlats,nlons))
lat2 = lat1.reshape((nlats,nlons))
print (lon2)
#x,y = map(lon2,lat2)
if var == 'phi':
	#cs = map.contourf(x,y,field/9.81,levels=np.arange(4600,6200,100),cmap='seismic')
	plt.contourf(field,10,cmap='seismic')
else:
	cs = map.contourf(x,y,field,10,cmap='seismic')
plt.colorbar(shrink=0.5)
plt.title(var+' 500 hpa '+nstep+'h - T'+trunc+' expid='+expid)
plt.savefig('test.png')
plt.show()
