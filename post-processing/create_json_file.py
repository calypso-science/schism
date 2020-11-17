import datetime
import sys
sys.path.append('/home/remy/Calypso/Software/schism-dev/pyschism') 
from netCDF4 import Dataset
from ttide.t_tide import t_tide
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf
from ttide import t_utils as tu
import xarray as xr
import os
import numpy as np
from matplotlib.dates import num2date,date2num
import copy
import pandas as pd
from scipy.interpolate import griddata
from pyschism.mesh import Hgrid
from pyproj import Proj,transform
from descartes.patch import PolygonPatch
from shapely.ops import cascaded_union
from shapely.geometry import Point, Polygon
import pandas as pd
import rfc3339
import geojson
        # dset = xr.Dataset(dset_dict)
        # dset.attrs = dict(type="SCHISM standard output file")
        # dset.to_netcdf(path=join(self.nest.outdir, self.filename_dict[outtype]), format='NETCDF4')
# -------------------------------------------- BASE MODEL CLASSES -----------------------------
WGS84 = 4326 # Spatial Reference System
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
    def __init__(self,lim, resolution,filgrid, **kwargs):
        super(MakeMeshMask, self).__init__(**kwargs)

        hgrid_proj = 2193
        self.mesh = self.load_mesh(filgrid)
        x, y = self.get_coords(filgrid)
 
        if  hgrid_proj != WGS84:
            self.lon, self.lat = transform_proj(x, y, in_espg=hgrid_proj, out_espg=WGS84)
        else:
            self.lon, self.lat = x, y 

        self.resolution = resolution
        self.lim= lim
        
    def load_mesh(self,filgrid):
        hgrid = Hgrid.open(filgrid,crs="EPSG:%i" % 2193)
        return hgrid

    def get_coords(self,filgrid):
        mesh = self.load_mesh(filgrid)

        x = np.array([ row[1][0][0] for row in mesh.nodes])
        y = np.array([ row[1][0][1] for row in mesh.nodes])
        return x, y

    def get_boundary_segments(self, type):
        '''Docstring'''

        coast_segment, obc_segment = list(), list()
        bnd_segment, island_segment = list(), list()
        #
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
        '''Docstring'''
        #segments=np.flipud(segments)
        ordered_nodes = list()
        consume_nodes = segments.tolist()
        #import pdb;pdb.set_trace()
        #ordered_nodes.append(consume_nodes[1])
        #ordered_nodes.append(consume_nodes[0])
        # ordered_nodes.append(consume_nodes[0])
        # ordered_nodes.append(consume_nodes[2])

        while consume_nodes:
            for s, segment in enumerate(segments):
                if not ordered_nodes:
                    print(1)
                    ordered_nodes.append(consume_nodes.pop(consume_nodes.index(segment)))
                    continue
                for seg in consume_nodes:
                    if seg[0] == ordered_nodes[-1][-1]:
                        ordered_nodes.append(consume_nodes.pop(consume_nodes.index(seg)))
                        print('bon')
                        break
                for seg in consume_nodes:
                    #import pdb;pdb.set_trace()
                    if seg[0] == ordered_nodes[-1][-1]+1:
                        ordered_nodes.append(consume_nodes.pop(consume_nodes.index(seg)))
                        break
                if len(consume_nodes)==1:
                    ordered_nodes.append(consume_nodes.pop(0))
                    break


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

        return Polygon(mesh_nodes, island_nodes) #

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

def extract_HC(consfile,Vars, lim, conlist=None,min_depth=1):
    """
    Extract harmonic constituents and interpolate onto points in lon,lat
    set "z" to specifiy depth for transport to velocity conversion
    set "constituents" in conlist
    Returns:
        u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist
    """
    ###
    # Read the filenames from the model file
    pathfile = os.path.split(consfile)
    path = pathfile[0]

    f = Dataset(consfile,'r')

    ###
    # Read the grid file
    X=f.variables['Easting'][:]
    Y=f.variables['Northing'][:]

    if 'Depth' in f.variables:
        Z=f.variables['Depth'][:]
    else:
        Z=f.variables['dep'][:]

    X, Y = transform_proj(X, Y, in_espg=2193, out_espg=WGS84)
    gd=(X>=lim[0]) & (X<=lim[1]) & (Y>=lim[2]) & (Y <=lim[3]) & (Z>=min_depth )

    X=X[gd]
    Y=Y[gd]
    Z=Z[gd]
    ###
    # Check that the constituents are in the file
    conList = []
    conIDX=[]
    if conlist is not None:
        for ncon in range(0,len(f.variables['cons'])):
            if ''.join([y.decode('UTF-8') for y in f.variables['cons'][ncon] if y!=b' ']) in conlist:
                x=''
                conList.append(''.join([x+n.decode('UTF-8') for n in f.variables['cons'][ncon]]))
                conIDX.append(ncon)
    else:
        for ncon in range(0,len(f.variables['cons'])):
            x=''
            conList.append(''.join([x+n.decode('UTF-8') for n in f.variables['cons'][ncon].data]))
            conIDX.append(ncon)

    const = t_getconsts(np.array([]))
    Const= [con.decode('UTF-8') for con in const[0]['name']] 

    consindex = [Const.index(con.ljust(4)) for con in conList]

    tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.

    if 'lev' in f.variables:
        levs=f.variables['lev'][:]
        Ilev=(levs==0).nonzero()[0]
        if len(Ilev)==0:
            print('Lev %f not found' %lev)
            Ilev=0
        else:
            Ilev=Ilev[0]
    else:
        Ilev=0
    

    var={}
    # interpolating to ROMS grid requires special care with phases!!
    #    this is subjected to parallelization via muitiprocessing.Process - do it!
    #    there is a sandbox in roms repo with a starting exercise
    for var0 in Vars:
        var[var0]=np.ones(shape=(len(tfreq),4,len(X)))*-1e-8
        N=-1
        for con,ncon in zip(const[0]['name'][consindex],conIDX):
            N=N+1
            con=con.decode('UTF-8')
            print("extracting %s for %s" %(var0, con))   
            amp =  f.variables[var0+"_amp"][ncon]
            pha = f.variables[var0+"_amp"][ncon]

            if len(amp.shape)>1:
                amp=amp[Ilev,gd]
                pha=pha[Ilev,gd]
            else:
                amp=amp[gd]
                pha=pha[gd]


            var[var0][N,0,:] = amp 
            var[var0][N,2,:] = pha#*180./np.pi



    return var,tfreq,consindex,X,Y,Z

def get_tide(ju,freq,tidecon0,t_time,lat0):
    tidecon=copy.deepcopy(tidecon0)
    nodes=tidecon.shape[2]

    if t_time.dtype.name.startswith('datetime64') or t_time.dtype is np.dtype("O"):
        t_time = tm.date2num(t_time)

    t_time = t_time.reshape(-1, 1)



    snr = (tidecon[:, 0,0] / tidecon[:, 1,0]) ** 2
    I = snr > 2

    
    tidecon = tidecon[I, :]
    ju = np.asarray(ju)[I]
    freq = freq[I]


    ap = np.multiply(tidecon[:, 0] / 2.0,np.exp(-1j * tidecon[:, 2] * np.pi / 180))
    am = np.conj(ap)


    jdmid = np.mean(t_time[0:((2 * int((max(t_time.shape) - 1) / 2)) + 1)])

    #import pdb;pdb.set_trace()
    v, u, f = t_vuf('nodal', [jdmid], ju, lat0)

    ap = ap * np.kron(np.ones((ap.shape[1],1)),f * np.exp(+1j * 2 * np.pi * (u + v))).T
    am = am * np.kron(np.ones((ap.shape[1],1)),f * np.exp(-1j * 2 * np.pi * (u + v))).T

    t_time = t_time - jdmid

    n, m = t_time.shape
    yout = np.zeros([ap.shape[1], 1], dtype='complex128')
    touter = np.outer(24 * 1j * 2 * np.pi * freq, t_time[0])

    touter=np.kron(np.ones((1,ap.shape[1])),touter)

    yout[:,0]=np.sum(np.multiply(np.exp(touter), ap), axis=0)+np.sum(np.multiply(np.exp(-touter), am), axis=0)

    tide=np.real(yout)


    return tide

def data2geojson(method,df,fileout):

    features = []
    if method=='velocity':

        insert_features = lambda X: features.append(
                geojson.Feature(geometry=geojson.Point((X["lon"],
                                                        X["lat"])),
                                 properties=dict(timestamp=rfc3339.format(X["time"],utc=True),
                                                 U=X["u"],V=X["v"])))
    else:
        insert_features = lambda X: features.append(
            geojson.Feature(geometry=geojson.Point((X["lon"],
                                                    X["lat"])),
                             properties=dict(timestamp=rfc3339.format(X["time"],utc=True),
                                             E=X["e"])))       
    df.apply(insert_features, axis=1)
    with open(fileout, 'w', encoding='utf8') as fp:
        geojson.dump(geojson.FeatureCollection(features), fp, sort_keys=True, ensure_ascii=False)

def process(fileout,gridfile,consfile,tstart,tend,dt,method,res,Cons,levs,min_depth,lim):

    TimeSeries=np.arange(date2num(tstart),date2num(tend)-(dt/24.),dt/(24.))


    msk=MakeMeshMask(lim,res,gridfile)
    Mask=msk.make_grid_mask()
    mask = Mask.mask.values
    if lim == None:
        lim=[Mask.ugrid_lon.values.min(), Mask.ugrid_lon.values.max(), \
             Mask.ugrid_lat.values.min(), Mask.ugrid_lat.values.max()]

    rloni, rlati = np.meshgrid(Mask.rgrid_lon.values, Mask.rgrid_lat.values)
    rgrid = (rloni, rlati)

    if method=='velocity':
        params=['um','vm']
    else:
        params=['e']


    var,tfreq,consindex,X,Y,Z=extract_HC(consfile,params, lim, conlist=Cons,min_depth=min_depth)
    ugrid = (X, Y)
    idx = pd.period_range(tstart,tend,freq='%is'%int(dt*3600))
    idx=idx.to_timestamp()
    

    df=pd.DataFrame()
    for n in range(0,len(TimeSeries)):
        df0=pd.DataFrame()
        for i,v in enumerate(var): 

            if n % 10 ==0:
                print('doing variable %s timestep %i / %i' % (v,n,len(TimeSeries)))

            tmp=get_tide(consindex,tfreq,var[v],np.array(TimeSeries[n]),np.mean(rlati[0,:]))
            tmp[tmp[:,0]<=-990,0]=np.nan
            tmp=griddata(ugrid, tmp[:,0], rgrid, method='linear')

            #masked_vari = np.ma.masked_array(tmp, mask=mask.reshape(tmp.shape))
            Z=tmp.ravel()#masked_vari.data[:].flatten()
            #import pdb;pdb.set_trace()
            X=rloni.ravel()
            Y=rlati.ravel()
            gd=(~mask) 
            X=X[gd]
            Y=Y[gd]
            Z=Z[gd]
            gd=(~np.isnan(Z))
            X=X[gd]
            Y=Y[gd]
            Z=Z[gd]            

            if i==0:
                df0['lon']=pd.Series(X)
                df0['lat']=pd.Series(Y)
                df0['time']=np.tile(idx[n],(len(X),1))

            df0[v.replace('m','')]=pd.Series(Z)
        df0.set_index('time',inplace=True)
        df=df.append(df0)
    df.reset_index(inplace=True,drop=False)
    data2geojson(method,df,fileout)




if __name__ == "__main__":
    
    # import argparse
    # parser = argparse.ArgumentParser(prog='cons2ts.py', usage='%(prog)s fileout consfile params tstart tend [options]')
    # ## main arguments
    # parser.add_argument('fileout', type=str,help='name of the output tidal file')
    # parser.add_argument('consfile', type=str,help='Folder with all the SCHSIM files')
    # parser.add_argument('params', type=str,nargs='+',help='Parameters (i.e e um vm')
    # parser.add_argument('tstart', type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),help='start time')
    # parser.add_argument('tend', type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),help='end time')
    # parser.add_argument('res', type=float,help='resolution')
    # parser.add_argument('gridfile',type=str,help='schism grid')
    # #'/home/remy/Calypso/Projects/003_HG_America/hgrid2.gr3'


    # ## options  
    # parser.add_argument('--dt', type=float,default=1,help='timestep in hour')
    # parser.add_argument('--cons', type=str,nargs='+',help='Constituents (i.e M2 S2')
    # parser.add_argument('--lev', type=float,nargs='+',default=0,help='outputs level 0 -1 -2 ... ')
    # parser.add_argument('--lim', type=float,help='Xlim,Ylim',nargs='+')
    # parser.add_argument('--min_depth', type=float,default=1,help='min_depth')

    # args = parser.parse_args()
    

    # process(args.fileout,\
    #     args.gridfile,\
    #     args.consfile,\
    #     args.tstart,args.tend,args.dt,\
    #     args.params,\
    #     args.res,\
    #     args.cons,\
    #     args.lev,\
    #     args.min_depth,
    #     args.lim)

    # process('E_bluff.geojson',
    # '/home/remy/Calypso/Projects/006_Southport/hydro/child/hg6_dt15_NEW/hgrid.gr3',
    # '/home/remy/Calypso/Projects/006_Southport/hydro/child/SP3D.nc',
    # '2020-01-01','2020-01-02',6,
    # 'elev',0.05,None,0,2,None)

    # process('UV_bluff.geojson',
    #     '/home/remy/Calypso/Projects/006_Southport/hydro/child/hg6_dt15_NEW/hgrid.gr3',
    #     '/home/remy/Calypso/Projects/006_Southport/hydro/child/SP3D.nc',
    #     '2020-01-01','2020-01-02',6,
    #     'velocity',0.01,None,0,2,None)

    # process('UV_Foveaux.geojson',
    # '/home/remy/Calypso/Projects/006_Southport/hydro/mother/schism/hgrid1.gr3',
    # '/home/remy/Calypso/Projects/006_Southport/hydro/mother/schism/run5/outputs/foveaux.nc',
    # '2020-01-01','2020-01-02',6,
    # 'velocity',0.05,None,0,2,None)


    process('E_Foveaux.geojson',
    '/home/remy/Calypso/Projects/006_Southport/hydro/mother/schism/hgrid1.gr3',
    '/home/remy/Calypso/Projects/006_Southport/hydro/mother/schism/run5/outputs/foveaux.nc',
    '2020-01-01','2020-01-02',6,
    'elev',0.1,None,0,2,None)