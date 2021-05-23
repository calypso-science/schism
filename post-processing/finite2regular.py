import datetime
import sys
import logging
from netCDF4 import Dataset
import netCDF4
import xarray as xr
import os
import numpy as np
import copy
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from pyschism.mesh import Hgrid
from pyproj import Proj,transform
from descartes.patch import PolygonPatch
from shapely.ops import cascaded_union
from shapely.geometry import Point, Polygon
import glob
import re
try:
    import interpz
except:
    print('no vertical interpolation')
        # dset = xr.Dataset(dset_dict)
        # dset.attrs = dict(type="SCHISM standard output file")
        # dset.to_netcdf(path=join(self.nest.outdir, self.filename_dict[outtype]), format='NETCDF4')

schism_mapping = {
    'elev' : 'ssh',
    'WWM_1' : 'hs', # time-varying vertical coordinates
    'WWM_11': 'tp',
    'WWM_18':'dpm',
    }


# -------------------------------------------- BASE MODEL CLASSES -----------------------------
WGS84 = 4326 # Spatial Reference System
def set_up_logger(logdir='/home/dataout'):
    if not os.path.isdir(logdir):
        os.makedirs(logdir)

    logger=logging.getLogger('logger')
    logger.setLevel(10)
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    fh = logging.FileHandler(os.path.join(logdir, 'log_transform.txt'), mode='w', encoding='utf-8')
    fh.setLevel(10)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger



#-------------------------------------------- GRID METHODS -----------------------------
def transform_proj(x, y, in_espg=2193, out_espg=4326):
    ''' Performs cartographic transformations (converts from
        longitude,latitude to native map projection x,y coordinates and
        vice versa) using proj (http://trac.osgeo.org/proj/).
    '''
    in_proj = Proj(init='epsg:%s'%in_espg)
    out_proj = Proj(init='epsg:%s'%out_espg)
    x2, y2 = transform(in_proj, out_proj, x, y)
    return x2, y2

class MakeMeshMask():   
    '''Docstring'''
    def __init__(self,lim, resolution,filgrid,hgrid_proj, **kwargs):
        super(MakeMeshMask, self).__init__(**kwargs)

        self.mesh = self.load_mesh(filgrid,hgrid_proj)
        x, y = self.get_coords(filgrid,hgrid_proj)
 
        if  hgrid_proj != WGS84:
            self.lon, self.lat = transform_proj(x, y, in_espg=hgrid_proj, out_espg=WGS84)
        else:
            self.lon, self.lat = x, y 

        self.resolution = resolution
        self.lim= lim
        
    def load_mesh(self,filgrid,epsg):
        hgrid = Hgrid.open(filgrid,crs="EPSG:%i" % epsg)
        return hgrid

    def get_coords(self,filgrid,epsg):
        mesh = self.load_mesh(filgrid,epsg)


        x = np.array(mesh.x)
        y = np.array(mesh.y)
        return x, y

    def get_boundary_segments(self, type):
        '''Docstring'''

        coast_segment, obc_segment = list(), list()
        bnd_segment, island_segment = list(), list()
        #
        # [island_segment.append(x) for x in self.mesh.boundaries._interior['indexes']]
        # [coast_segment.append(x) for x in self.mesh.boundaries._land['indexes']]
        # [bnd_segment.append(x) for x in self.mesh.boundaries._land['indexes']]
        # [obc_segment.append(x) for x in self.mesh.boundaries._ocean['indexes']]
        # [bnd_segment.append(x) for x in self.mesh.boundaries._ocean['indexes']]

        for flag in self.mesh.boundaries:
            for bnd in self.mesh.boundaries[flag]:
                segment = self.mesh.boundaries[flag][bnd]['indexes']
                segment=[int(x)-1 for x in segment]

                if flag==0: # coastline
                    coast_segment.append(segment)
                    bnd_segment.append(segment)
                elif flag==1: #islands
                    island_segment.append(segment)
                else:# 'open boundary' in flag: # ocean
                    obc_segment.append(segment)
                    bnd_segment.append(segment)
        if type == 'mesh_edge':
            return np.array(bnd_segment)
        elif type == 'obc':
            return np.array(obc_segment)
        elif type == 'coastline':
            return np.array(coast_segment)
        elif type == 'island':
            return np.array(island_segment)


    def order_segments(self, segments):
        consume_nodes = segments.tolist()
        # import pdb;pdb.set_trace()
        # ordered_nodes=[consume_nodes[3]]+[consume_nodes[6]]+\
        # [consume_nodes[1]]+[consume_nodes[7]]+\
        # [consume_nodes[0]]+[consume_nodes[8]]+\
        # [consume_nodes[2]]+[consume_nodes[5]]+\
        # [consume_nodes[4][::-1]]+[consume_nodes[9][::-1]]
        ordered_nodes=[consume_nodes.pop(0)]
        I=0
        # import pdb;pdb.set_trace()
        while consume_nodes and I<1000:
            I+=1
            for s,seg in enumerate(consume_nodes):
                if ordered_nodes[-1][-1]==seg[0]:
                    nodes=consume_nodes.pop(s)
                    ordered_nodes.append(nodes)
                elif ordered_nodes[-1][-1]==seg[-1]:
                    nodes=consume_nodes.pop(s)[::-1]
                    ordered_nodes.append(nodes)
                elif ordered_nodes[0][0]==seg[0]:
                    nodes=consume_nodes.pop(s)[::-1]
                    ordered_nodes=[nodes]+ordered_nodes
                elif ordered_nodes[0][0]==seg[-1]:
                    nodes=consume_nodes.pop(s)
                    ordered_nodes=[nodes]+ordered_nodes                  


        return ordered_nodes

    def get_mesh_polygon(self):
        '''Docstring'''

        # loading mesh edge segments
        mesh_edge = self.get_boundary_segments(type='mesh_edge')
        mesh_edge = self.order_segments(mesh_edge)
        mesh_nodes = list()

        for seg in mesh_edge:
            for node in seg:
                mesh_nodes.append((self.lon[node],self.lat[node]))

        # load  island edge segments
        island_edge = self.get_boundary_segments(type='island')
        island_nodes = list()
        for seg in island_edge:
            island_nodes.append([(self.lon[node],self.lat[node]) for node in seg])

        return Polygon(mesh_nodes,island_nodes ) #

    def _inpolygon(self, polygon, xp, yp):
        return np.array([Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],dtype=np.bool)

    def _outpolygon(self, polygon, xp, yp):
        return ~np.array([Point(x, y).intersects(polygon) for x, y in zip(xp, yp)],dtype=np.bool)

    def make_grid_mask(self):
        '''Docstring'''

        polygon = self.get_mesh_polygon()
        if self.lim == None:
            xi = np.arange(self.lon.min(),self.lon.max(),self.resolution)
            yi = np.arange(self.lat.min(),self.lat.max(),self.resolution)
        else:
            xi = np.arange(self.lim[0],self.lim[1],self.resolution)
            yi = np.arange(self.lim[2],self.lim[3],self.resolution)

        xig, yig = np.meshgrid(xi,yi)

        mask = self._outpolygon(polygon, xig.ravel(), yig.ravel())
        ds = xr.Dataset({'mask': (['rgrid_size'], mask.astype(bool))},
                         coords={'rgrid_size':(['rgrid_size'], np.arange(xig.size)),
                                 'rgrid_lon': (['rgrid_lon'], xi),
                                 'rgrid_lat': (['rgrid_lat'], yi),
                                 'ugrid_lon': (['ugrid_lon'], self.lon),
                                 'ugrid_lat': (['ugrid_lat'], self.lat)})
        ds.attrs['Grid resolution [deg]'] = self.resolution

        outfilestr = '/static/{mode}/schism/{imp}/regular_grid-mask_{res}.nc'

        return ds
def vertical_interpolation(zcor,e,lev):
    E=np.zeros((e.shape[0],len(lev)))

    if lev[0]==0.:
        E[:,0]=e[:,-1]
        lev2=lev[1:]
    else:
        lev2=lev
    
    e2=e.data
    e2[e.data==e.fill_value]=0
    z2=zcor.data


    z2[z2>1e36]=-99999
    z2[np.isnan(z2)]=-99999
    surf=np.zeros((zcor.shape[0],1))
    surf[:,0]=z2[:,-1]
    z2=z2-surf
    z2[z2<-90999]=np.nan
    for la in range(z2.shape[1]-2,-1,-1):
        ind = np.where(np.isnan(z2[:,la]))[0]
        z2[ind,la] = z2[ind,la+1]


    EE=interpz.interpz1d(e2,z2,lev2,np=e2.shape[0],nzin=e2.shape[1],nzout=len(lev2),kz=1, null_value=-9.9e15)
    
    for n in range(0,len(lev2)):
        E[:,n+1]=EE[:,n]


    E[E==-99999]=np.nan
    return E

def create_dataset(times,unit,X,Y,Vars,depth,lev=0):
    #import pdb;pdb.set_trace()
    # create dummy dataframe
    dset={}
    dset['time']=xr.DataArray(data=times,
                            dims=['time'],
                            attrs = {"units": unit
                            }
                            )



    dset['dep']=xr.DataArray(
            data   = depth,   # enter data here
            dims   = ['lat','lon'],
            coords = {"lat": (["lat"], Y[:,0]),
                     "lon": (["lon"], X[0,:])},                        
            attrs  = {
                '_FillValue': 1e20,
                'units'     : 'm',
                'standard_name': 'sea_floor_depth_below_mean_sea_level',
                'long_name': 'Depth',
                }
            )


    for var in Vars:
        if var == 'elev':
            dset['ssh']=xr.DataArray(
                    data   = np.random.random((len(times),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time': times,
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},                        
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm',
                        'long_name': 'Sea surface height',
                        'standard_name': 'sea_surface_height_above_geoid',
                        }
                    )

        if var == 'WWM_1':
            dset['hs']=xr.DataArray(
                    data   = np.random.random((len(times),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time': times,
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},                        
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm',
                        'long_name': 'significant height of wind and swell wave',
                        'standard_name': 'sea_surface_wave_significant_height',
                        }
                    )

        if var == 'WWM_11':
            dset['tp']=xr.DataArray(
                    data   = np.random.random((len(times),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time': times,
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},                        
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 's',
                        'long_name': 'discrete peak period',
                        'standard_name': 'discrete_peak_period',
                        }
                    )

        if var == 'WWM_18':
            dset['dpm']=xr.DataArray(
                    data   = np.random.random((len(times),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time': times,
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},                        
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'deg',
                        'long_name': 'peak direction',
                        'standard_name': 'sea_surface_wave_peak_direction',
                        }
                    )



        if var == 'salt': 
            dset['salt']=xr.DataArray(
                    data   = np.random.random((len(times),len(lev),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lev','lat','lon'],
                    coords = {'time': times,
                             "lev": (['lev'],lev),
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'PSU',
                        'long_name': 'Water salinity',
                        'standard_name': 'sea_water_salinity',
                        }
                    )
        if var == 'temp': 
            dset['temp']=xr.DataArray(
                    data   = np.random.random((len(times),len(lev),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lev','lat','lon'],
                    coords = {'time': times,
                             "lev": (['lev'],lev),
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'degC',
                        'long_name': 'Water temperature',
                        'standard_name': 'sea_water_temperature',
                        }
                    )
        
        if var == 'vertical_velocity': 
            dset['vertical_velocity']=xr.DataArray(
                    data   = np.random.random((len(times),len(lev),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lev','lat','lon'],
                    coords = {'time': times,
                             "lev": (['lev'],lev),
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm.s^{-1}',
                        'long_name': 'vertical velocity',
                        'standard_name': 'vertical_velocity',
                        }
                    )
        if var == 'TKE': 
            dset['TKE']=xr.DataArray(
                    data   = np.random.random((len(times),len(lev),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lev','lat','lon'],
                    coords = {'time': times,
                             "lev": (['lev'],lev),
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : '',
                        'long_name': 'turbulent kinetic energy',
                        'standard_name': 'turbulent_kinetic_energy',
                        }
                    )

        if var == 'hvel':
            dset['u']=xr.DataArray(
                    data   = np.random.random((len(times),len(lev),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lev','lat','lon'],
                    coords = {'time': times,
                             "lev": (['lev'],lev),
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm.s^{-1}',
                        'long_name': 'Eastward current',
                        'standard_name': 'eastward_sea_water_velocity',
                        }
                    )


            dset['v']=xr.DataArray(
                    data   = np.random.random((len(times),len(lev),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lev','lat','lon'],
                    coords = {'time': times,
                             "lev": (['lev'],lev),
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm.s^{-1}',
                        'long_name': 'Northward current',
                        'standard_name': 'northward_sea_water_velocity',
                        }
                    )

        if var == 'dahv':
            dset['um']=xr.DataArray(
                    data   = np.random.random((len(times),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time': times,
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm.s^{-1}',
                        'long_name': 'Eastward depth-averaged current',
                        'standard_name': 'barotropic_eastward_sea_water_velocity',
                        }
                    )



            dset['vm']=xr.DataArray(
                    data   = np.random.random((len(times),X.shape[0],X.shape[1])),   # enter data here
                    dims   = ['time','lat','lon'],
                    coords = {'time': times,
                             "lat": (["lat"], Y[:,0]),
                             "lon": (["lon"], X[0,:])},  
                    attrs  = {
                        '_FillValue': 1e20,
                        'units'     : 'm.s^{-1}',
                        'long_name': 'Northward depth-averaged current',
                        'standard_name': 'barotropic_nortward_sea_water_velocity',
                        }
                    )


    ds = xr.Dataset(dset,attrs={'type':"SCHISM standard output file"})
    return ds

def extract_raw(file,params, lim,idx,lev):

    nc=netCDF4.Dataset(file)
    nt=len(nc.variables['time'])

    ds={}
    for n in range(0,nt):      
        for param in params:

            data=nc.variables[param][n]
            data = np.ma.masked_where(data==-9999,data)
            i23d=nc.variables[param].i23d
            ivs=nc.variables[param].ivs

            
            if ivs==2:
                if param=='dahv':
                    Param=['um','vm']
                if param=='hvel':
                    Param=['u','v']
            else:
                Param=[param]

            if n==0:
                ds[Param[0]]=np.ones((nt,len(idx.nonzero()[0])))
                if i23d>1:
                    ds[Param[0]]=np.ones((nt,len(lev),len(idx.nonzero()[0])))
                    if ivs==2:
                        ds[Param[1]]=np.ones((nt,len(lev),len(idx.nonzero()[0])))
                elif ivs==2:
                    ds[Param[0]]=np.ones((nt,len(idx.nonzero()[0])))
                    ds[Param[1]]=np.ones((nt,len(idx.nonzero()[0])))

            if ivs==1 and i23d==1: # elev
                ds[Param[0]][n,:]=data[idx]
            if ivs==2 and i23d==1: # dahv
                ds[Param[0]][n,:]=data[idx,0]
                ds[Param[1]][n,:]=data[idx,1]

            if i23d>1: # zcor
                zcor=nc.variables['zcor'][n,:,:]
                zcor=zcor[idx]
                if ivs==1: # temp
                    tmp=data[idx]
                    tmp=vertical_interpolation(zcor,tmp,lev)
                    ds[Param[0]][n,:,:]=tmp.T
                else:
                    for k in range(0,2):
                        tmp=data[idx,:,k]
                        tmp=vertical_interpolation(zcor,tmp,lev)
                        ds[Param[k]][n,:,:]=tmp.T

    return ds,nt

def read_initial_netcdf_file(file0,file1,epsg,lim,min_depth):

    inProj = Proj(init='epsg:'+str(epsg))
    outProj = Proj(init='epsg:'+str(WGS84))
    nc=netCDF4.Dataset(file0)
    nc1=netCDF4.Dataset(file1)

    unit=nc.variables['time'].units
    depth=nc.variables['depth'][:]
    X=nc.variables['SCHISM_hgrid_node_x'][:]
    Y=nc.variables['SCHISM_hgrid_node_y'][:]
    X,Y=transform(inProj,outProj,X[:],Y[:]) 

    bufX=(lim[1]-lim[0])*5/100
    bufY=(lim[3]-lim[2])*5/100

    gd=(X>=lim[0]-bufX) & (X<=lim[1]+bufX) & (Y>=lim[2]-bufY) & (Y <=lim[3]+bufY) & (depth>=min_depth )

    X=X[gd]
    Y=Y[gd]
    depth=depth[gd]


    t0 = nc.variables['time'][0] #netCDF4.num2date(nc.variables['time'][0],unit)
    t1 = nc1.variables['time'][-1] #netCDF4.num2date(nc1.variables['time'][-1],unit)



    
    dt=(nc.variables['time'][2]-nc.variables['time'][1])


    times = np.arange(t0,t1+dt,dt)# pd.date_range(start=t0.strftime(),freq='%fH'%dt,end=t1.strftime())


    nc.close()
    nc1.close()

    return times,X,Y,depth,gd,unit
def get_INDstart(dirout,prefix):
    all_file=glob.glob(os.path.join(dirout,prefix+'*'))
    all_file.sort(key=lambda f: int(re.sub('\D', '', f)))
    return int(all_file[0].split('_')[-1].replace('.nc',''))

def get_INDend(dirout,prefix):
    all_file=glob.glob(os.path.join(dirout,prefix+'*'))  
    all_file.sort(key=lambda f: int(re.sub('\D', '', f)))

    return int(all_file[-1].split('_')[-1].replace('.nc',''))



def process(fileout,hgrid,dirout,INDstart,INDend,params,res,levs,min_depth,lim,prefix,epsg,log_dir=False):
    if log_dir:
        logger=set_up_logger(log_dir)
    else:
        logger=False


    if type(levs)!=type([]):
        levs=[levs]


    if INDstart==0:
        INDstart=get_INDstart(dirout,prefix)

    if INDend==0:
        INDend=get_INDend(dirout,prefix)

    if logger:
        logger.info('\tReading the mask: %s' % (hgrid))

    msk=MakeMeshMask(lim,res,hgrid,epsg)
    Mask=msk.make_grid_mask()
    mask = Mask.mask.values
    if lim == None:
        lim=[Mask.ugrid_lon.values.min(), Mask.ugrid_lon.values.max(), \
             Mask.ugrid_lat.values.min(), Mask.ugrid_lat.values.max()]

    rloni, rlati = np.meshgrid(Mask.rgrid_lon.values, Mask.rgrid_lat.values)
    rgrid = (rloni, rlati)

    if logger:
        logger.info('\tReads Netcdf file from %i to %i' % (INDstart,INDend))

    Ts,X,Y,depth,gd,unit=read_initial_netcdf_file(os.path.join(dirout,prefix+str(INDstart)+'.nc'),\
                                          os.path.join(dirout,prefix+str(INDend)+'.nc'),\
                                          epsg,lim,min_depth)

    ugrid=(X,Y)
    interpo = LinearNDInterpolator(ugrid, depth)
    tmp=interpo(rgrid)
    masked_vari = np.ma.masked_array(tmp, mask=mask.reshape(tmp.shape))


    if logger:
        logger.info('\tInitialize' )

    ds=create_dataset(Ts,unit,rloni,rlati,params,masked_vari,lev=levs)
    ds['dep']=ds['dep'][:].fillna(1e20)

    

    for nfile in range(INDstart,INDend+1):
        fileIN=os.path.join(dirout,prefix+str(nfile)+'.nc')
        if logger:
            logger.info('\tReading file %s' % fileIN )      
        Z,nt=extract_raw(fileIN,params, lim, gd,lev=levs)

        for n in range(0,nt):
            nn=n+(nt*(nfile-INDstart))
            if logger:
                logger.info('\t\tDoing timestep %.i' % n )         

            for vv in Z:
                if logger:
                    logger.info('\t\t\tDoing variable %s' % vv ) 
                if vv in schism_mapping.keys():
                    v=schism_mapping[vv]
                else:
                    v=copy.deepcopy(vv)



                if len(Z[vv].shape)==3:
                    for i in range(0,Z[vv].shape[1]):
                        Z[vv][n,i,Z[vv][n,i,:]<=-990]=np.nan
                        interpo.values[:,0]=Z[vv][n,i,:]
                        masked_vari = np.ma.masked_array(interpo(rgrid), mask=mask.reshape(tmp.shape))
                        import pdb;pdb.set_trace()
                        # if i>2:
                        #     masked_vari[masked_vari<=-990.0]=1e20
                        ds[v][nn,i,:,:]=masked_vari
                else:
                        interpo.values[:,0]=Z[vv][n,:]
                        masked_vari = np.ma.masked_array(interpo(rgrid), mask=mask.reshape(tmp.shape))
                        ds[v][nn,:,:]=masked_vari
            

                ds[v]=ds[v][:].fillna(1e20)


    fileout=fileout.replace('.nc','')+netCDF4.num2date(Ts,unit)[0].strftime('%Y%m%d_%Hz.nc')
    if logger:
        logger.info('\tSaving file %s '%  fileout)
    ds.to_netcdf(fileout)

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(prog='cons2ts.py', usage='%(prog)s fileout consfile params tstart tend [options]')
    ## main arguments
    parser.add_argument('fileout', type=str,help='name of the output tidal file')
    parser.add_argument('hgrid', type=str,help='name of the output tidal file')
    parser.add_argument('dirout', type=str,help='Folder with all the SCHSIM files')
    parser.add_argument('INDstart', type=int,help='First file to take')
    parser.add_argument('INDend', type=int,help='Last file to take')
    parser.add_argument('params', type=str,nargs='+',help='Parameters (i.e e um vm')
    parser.add_argument('res', type=float,help='resolution')

    ## options  
    parser.add_argument('--lev', type=float,nargs='+',default=0,help='outputs level 0 -1 -2 ... ')
    parser.add_argument('--lim', type=float,help='Xlim,Ylim',nargs='+')
    parser.add_argument('--min_depth', type=float,default=1,help='min_depth')
    parser.add_argument('--prefix', type=str,default='schout_',help='prefix')
    parser.add_argument('--epsg', type=int,default=2193,help='epsg')
    parser.add_argument('--log_dir', type=str,default=None,help='logger')
    
    args = parser.parse_args()

    process(args.fileout,\
        args.hgrid,\
        args.dirout,\
        args.INDstart,args.INDend,\
        args.params,\
        args.res,\
        args.lev,\
        args.min_depth,
        args.lim,
        args.prefix,
        args.epsg,
        args.log_dir)
