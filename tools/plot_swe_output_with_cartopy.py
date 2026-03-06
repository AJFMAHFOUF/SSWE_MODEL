import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
nlats=48 ; nlons=97
nstep='0072'
trunc='031'
expid='001'
var='phi'
plt.figure(figsize=(8,6))        
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

proj = ccrs.Orthographic(central_longitude=-90.0,central_latitude=40.0)
ax = plt.axes(projection=proj)
ax.coastlines('110m', linewidth=0.5)
#ax.coastlines()
#ax.add_feature(cfeature.LAND.with_scale('10m'), edgecolor='red', linewidth=0.2, facecolor='black', zorder=0)
#ax.add_feature(cfeature.LAKES, edgecolor='navy', linewidth=0.2, facecolor='none')
#ax.add_feature(cfeature.OCEAN)
#ax.add_feature(cfeature.NaturalEarthFeature(category='cultural',
#					name='admin_0_countries',
#					scale='110m',
#					edgecolor='black',
#					facecolor='none',
#					linewidth=0.2))
gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=1,color='black',alpha=0.5)
gl.ylocator = mticker.FixedLocator(np.arange(-90,90,30))
gl.xlocator = mticker.FixedLocator(np.arange(-180,180,30))

#x,y = map(lon2,lat2)

if var == 'phi':
	cbar_ticks = np.arange(4600,6100,100)
	#cs = map.contourf(x,y,field/9.81,levels=np.arange(4600,6200,100),cmap='jet')
	cs = ax.contourf(lon2,lat2,field/9.81,levels=15,transform=ccrs.PlateCarree(),cmap='jet')
else:
	cs = ax.contourf(lon2,lat2,field,levels=10,transform=ccrs.PlateCarree(),cmap='seismic')
#cbar = plt.colorbar(cs,orientation='horizontal',shrink=0.6,pad=0.073,extendrect=True)
#cbar.ax.tick_params(labelsize=8)
cbar = plt.colorbar(cs,orientation='vertical',shrink=0.5,pad=0.07)
cbar.ax.tick_params()
plt.title(var+' 500 hpa '+nstep+'h - T'+trunc+' expid='+expid, fontsize=17)
plt.savefig('../plots/'+var+'_500_T'+trunc+'_step_'+nstep+'_expid_'+expid+'.pdf')
plt.show()
