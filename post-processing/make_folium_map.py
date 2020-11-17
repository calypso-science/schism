import netCDF4
from tvtk.api import tvtk
from mayavi import mlab
import numpy
from scipy.interpolate import griddata,RegularGridInterpolator
from numpy import mgrid, empty, sin, pi, array, meshgrid, arange, prod
import rasterio
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from pyproj import Proj, transform
import fiona
from matplotlib.path import Path
from folium.plugins import TimestampedGeoJson
import mplleaflet
import folium
import branca
import dateutil
from datetime import timedelta
from matplotlib.dates import num2date,date2num
inProj = Proj(init='epsg:'+str(2193))
outProj = Proj(init='epsg:'+str(4326))


	# some path to your local netcdf file
#path = '/home/remy/Calypso/Projects/006_Southport/hydro/child/hg6_dt30_newD/outputs/schout_14.nc'
T0=dateutil.parser.parse('2019-09-02T12:00:00Z')
fileout='Fov.html'

path = '/home/remy/Calypso/Projects/006_Southport/hydro/mother/schism/run5/outputs/schout_14.nc'
#T0=dateutil.parser.parse('2020-01-03T09:00:00Z')#T0=datenum(2020,01,03,09,00,00)+0/24;ext=3;
#fileout='Neap_tide.html'

lim=[1241352-400,1245886+2500,4825135,4830428]
lim=[1110740,1304120,4704419,4874442]
#lim=[1236590,1254365,4821428,4836312]
quiver_res=300
quiver_res=50
quiver_scale=80
# try and import the dataset, prefer netcdf4, you might want to use pydap also here if you access a dap server.

ds = netCDF4.Dataset(path)
X=ds.variables['SCHISM_hgrid_node_x'][:]
Y=ds.variables['SCHISM_hgrid_node_y'][:]

X,Y=transform(inProj,outProj,X,Y)
lim[0],lim[2]=transform(inProj,outProj,lim[0],lim[2])
lim[1],lim[3]=transform(inProj,outProj,lim[1],lim[3])



off=0
gd = (X>=lim[0]-off) & (X<=lim[1]+off) & (Y>=lim[2]-off) & (Y<=lim[3]+off)
X=X[gd]
Y=Y[gd]
XY=numpy.array((X,Y)).T

d= ds.variables['depth'][:]
d=d[gd]
dtime = netCDF4.num2date(ds.variables['time'][:],ds.variables['time'].units)

dtime=[date2num(x._to_real_datetime()) for x in dtime]


idxTs=[]

for dt in numpy.arange(-6,6.5,1):
    to=date2num(T0+timedelta(hours=dt))
    idxTs.append(dtime.index(min(dtime, key=lambda x:abs(x-to))))
    



#D=D[gd]
Xreg, Yreg = numpy.meshgrid(numpy.linspace(lim[0],lim[1], quiver_res), numpy.linspace(lim[2], lim[3], quiver_res))




mapa = folium.Map(location=[Yreg.mean(), Xreg.mean()], tiles="Cartodb Positron",
                  zoom_start=10)

# shape = fiona.open('/home/remy/Calypso/Projects/006_Southport/hydro/child/animation/poly.shp')
verts=[]
# for poly in shape:
#     if poly['geometry']!=None :
#             nvert=len(poly['geometry']['coordinates'][0])
#             verts.append([(poly['geometry']['coordinates'][0][x][0],poly['geometry']['coordinates'][0][x][1]) for x in range(0,nvert)])

is_first=True

for i,idxT in enumerate(idxTs):
    print(i)
    if (i!=2) and (i!=9):
        continue
    
     # read the variables (these are just links to the objects for now)
    u = ds.variables['hvel'][idxT,:,-1,0] # velocity in u direction
    v = ds.variables['hvel'][idxT,:,-1,1] # v direction
    e = ds.variables['elev'][idxT,:] # v direction
    e=e[gd]
    u=u[gd]
    v=v[gd]
    bad=e+d<0.1

    U = griddata(XY[~bad,:], u[~bad],(Xreg,Yreg),method='linear')
    V = griddata(XY[~bad,:], v[~bad],(Xreg,Yreg),method='linear')
    E = griddata(XY[~bad,:], e[~bad],(Xreg,Yreg),method='linear')
    D = griddata(XY[~bad,:], d[~bad],(Xreg,Yreg),method='linear')

    fig = plt.figure(figsize=(30,18))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')

    mag=numpy.sqrt(U**2+V**2)*1.94

    U[numpy.isnan(U)]=0
    V[numpy.isnan(V)]=0
    mag[numpy.isnan(mag)]=0   


    #Xreg[numpy.isnan(Xreg)]=0
    #Yreg[numpy.isnan(Yreg)]=0

    Q = ax.quiver(Xreg, Yreg, U*1.94, V*1.94,mag,scale=quiver_scale,color=cm.viridis(64),clim=[0,4])


    gj = mplleaflet.fig_to_geojson(fig=fig)

    if i==2:
        name='Ebb tide'
    elif i==9:
        name='Flood tide'
    feature_group0 = folium.FeatureGroup(name=name,show=is_first)
    #'%.1fH' % ((dtime[idxTs[i]]-date2num(T0))*24)
    is_first=False

    spd=RegularGridInterpolator([Yreg[:,0],Xreg[0,:]],mag)

    for feature in gj['features']:
        if feature['geometry']['type'] == 'Point':
            lon, lat = feature['geometry']['coordinates']

            div = feature['properties']['html']
            flag=False
            for vert in verts:
                p=Path(vert)
                flag = p.contains_points([[lon,lat]])
                if flag:
                   break

            Spd=spd([lat,lon])[0]
            if not flag and Spd>0.05:
                #print(Spd)
                icon_anchor = (feature['properties']['anchor_x'],
                               feature['properties']['anchor_y'])

                icon = folium.features.DivIcon(html=div,
                                               icon_anchor=icon_anchor)
                

               
                #spd=
                marker = folium.Marker([lat, lon], icon=icon)
                
                marker.add_child(folium.Popup('%.1f' % Spd))
                



                feature_group0.add_child(marker)
        else:
            msg = "Unexpected geometry {}".format
            raise ValueError(msg(feature['geometry']))

    mapa.add_child(feature_group0)



tile = folium.TileLayer(
        tiles = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr = 'Esri',
        name = 'Esri Satellite',
        overlay = False,
        control = False
       ).add_to(mapa)


colormap = branca.colormap.linear.viridis.scale(0, 4)
colormap = colormap.to_step(index=numpy.arange(0,4.01,.01))
colormap.caption = 'Tidal speed [knt]'
colormap.add_to(mapa)


folium.LayerControl(collapsed=False).add_to(mapa)
mapa.save(fileout)


