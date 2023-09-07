from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
nlats=64 ; nlons=129
nstep='0000'
trunc='042'
expid='tr2'
var='qv'
#
zh0=2000.0
zr0=np.pi/9.0
zlonc=0.5*np.pi
zlatc=np.pi/6.0
plt.figure(figsize=(12,6))
map = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
#map = Basemap(projection='ortho',lat_0=40,lon_0=-90,resolution='l')           
#map.drawcoastlines(linewidth=0.25)
#map.drawcountries(linewidth=0.25)
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))  
lon1,lat1,vor,div,phi,u,v,psi,khi,qv=np.loadtxt('../data_out/SSWE_model_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.dat',unpack=True)
# define orography
zr2=(lon1*np.pi/180. - zlonc)**2 + (lat1*np.pi/180. - zlatc)**2
zr02=zr0*zr0*np.ones(nlats*nlons)
zr2=np.minimum(zr2,zr02)
zr=np.sqrt(zr2)
orog=zh0*(1 - zr/zr0)
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
relief = orog.reshape((nlats,nlons))
lon2 = lon1.reshape((nlats,nlons))
lat2 = lat1.reshape((nlats,nlons))
x,y = map(lon2,lat2)
print ('max=',np.max(field))
if var == 'phi':
	cs = map.contour(x,y,relief,linestyles='solid',levels=[100,500,1000,1500,2000],colors=['black','black','black','black','black'])
	cs = map.contour(x,y,field/9.81,levels=np.arange(5050,5950,50),cmap='hsv')
else:
	cs = map.contour(x,y,relief,linestyles='solid',levels=[100,500,1000,1500,1800],colors=['black','black','black','black','black'])
	cs = map.contour(x,y,field,linestyles='solid',levels=[0.05,0.25,0.5,0.75,0.9],colors=['red','red','red','red','red'])
	cs = map.contourf(x,y,field,levels=np.arange(0.0,1.1,0.1),cmap='YlGnBu')
	#cs = map.contourf(x,y,field,10,cmap='YlGnBu')
plt.colorbar(shrink=0.8)
plt.title(var+' 500 hpa '+nstep+'h - T'+trunc+' expid='+expid)
plt.xlabel('Longitudes')
plt.ylabel('Latitudes')
plt.savefig('../plots/'+var+'_500_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.pdf')
plt.show()
