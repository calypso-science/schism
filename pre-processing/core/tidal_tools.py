import os
import numpy as np
from scipy import interpolate
from scipy.spatial import Delaunay
from netCDF4 import Dataset
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf
from matplotlib.tri import Triangulation
from ttide.t_predic import t_predic
from interp2D import mask_interp
#from vcmq import regrid2d,create_grid,MV2,set_grid,griddata,grid2xy
#from vacumm.misc.grid.regridding import fill2d

import copy
from scipy.interpolate import interp1d

# def missing_interp(lon,lat,arr):
#     grid0 = create_grid(lon[0,:], lat[:,0])
#     varri = MV2.asarray(arr)
#     varri=set_grid(varri, grid0)
#     tempf = fill2d(varri, method='carg')
#     lon=np.ma.masked_where(tempf==tempf.fill_value, lon)
#     lat=np.ma.masked_where(tempf==tempf.fill_value, lat)
#     arri = griddata(lon.compressed(),lat.compressed(), tempf.compressed(), (lon[0,:], lat[:,0]), method='nat', ext=True, sub=10)

#     return arri




def extrapolate_amp(lon,lat,lonr, latr, maskr, arr):
    arr=mask_interp(lon,lat,arr)
    arri = arr(np.vstack((lonr,latr)).T, nnear=6, p=2)

    bad=np.isnan(arri)
    if np.any(bad):
        import pdb;pdb.set_trace()

    return arri
def extrapolate_pha(lon,lat,lonr, latr, maskr, arr):

    cx, cy = np.cos(arr), np.sin(arr)
    cx=mask_interp(lon,lat,cx)
    cy=mask_interp(lon,lat,cy)
    cx = cx(np.vstack((lonr,latr)).T, nnear=6, p=2)
    cy = cy(np.vstack((lonr,latr)).T, nnear=6, p=2)
    bad=np.logical_or(np.isnan(cy),np.isnan(cx))
    if np.any(bad):
        import pdb;pdb.set_trace()

    arr = np.arctan2(cy, cx)
    return arr

def extract_HC(modfile,Vars, lon, lat, conlist=None, logger=None):
    """
    Extract harmonic constituents and interpolate onto points in lon,lat
    set "z" to specifiy depth for transport to velocity conversion
    set "constituents" in conlist
    Returns:
        u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist
    """

    ###
    lon = np.asarray(lon)
    lat = np.asarray(lat)

    # Make sure the longitude is between 0 and 360
    lon = np.mod(lon,360.0)

    ###
    # Read the filenames from the model file
    pathfile = os.path.split(modfile)
    path = pathfile[0]


    f = Dataset(modfile,'r')

    ###
    # Read the grid file
    X=np.mod(f.variables['lon'][:],360)
    Y=f.variables['lat'][:]
    if len(X.shape)==1:
        X, Y = np.meshgrid(X, Y)

    if 'dep' in f.variables:
        depth=f.variables['dep'][:]
    else:
        depth=np.zeros((X.shape))

    
    ###
    # Check that the constituents are in the file
    conList = []
    conIDX=[]
    if conlist is not None:
        for ncon in range(0,len(f.variables['cons'])):
            if f.variables['cons'][ncon][:] in conlist:
                x=''
                conList.append(''.join([x+n for n in f.variables['cons'][ncon]]))
                conIDX.append(ncon)
    else:
        for ncon in range(0,len(f.variables['cons'])):
            x=''
            try:
                conList.append(''.join([x+n.decode('UTF-8') for n in f.variables['cons'][ncon].data]))
            except:
                conList.append(''.join([x+n for n in f.variables['cons'][ncon]]))
            conIDX.append(ncon)


    const = t_getconsts(np.array([]))
    Const= [con.decode('UTF-8') for con in const[0]['name']] 
    consindex = [Const.index(con.ljust(4).upper()) for con in conList]

    tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.

    ###
    # Now go through and read the data for each

    
    maskr = []
    
    Vars=list(set([x.replace('_amp','').replace('_pha','') for x in Vars]))
    var={}
    # interpolating to ROMS grid requires special care with phases!!
    #    this is subjected to parallelization via muitiprocessing.Process - do it!
    #    there is a sandbox in roms repo with a starting exercise
    for var0 in Vars:
        var[var0]=np.ones(shape=(len(tfreq),4,len(lon)))*-1e-8
        N=-1

        for con,ncon in zip(const[0]['name'][consindex],conIDX):
            N=N+1
            con=con.decode('UTF-8')
            logger.info("Interpolating %s for %s" %(var0, con))        
            var[var0][N,0,:] = extrapolate_amp(X, Y, lon, lat, maskr, f.variables[var0+"_amp"][ncon,:,:]) 
            var[var0][N,2,:] = extrapolate_pha(X, Y, lon, lat, maskr, f.variables[var0+"_pha"][ncon,:,:])*180./np.pi



    return var,tfreq,consindex

def get_tide(ju,freq,tidecon0,t_time,lat0):
    tidecon=copy.deepcopy(tidecon0)
    nodes=tidecon.shape[2]

    
    #tidecon[:,0,:]=tidecon[:,0,:]+tidecon[:,0,:]*50/100

    #print('!!!!!!ADDING 50%!!!!!!!!!!!!!!')




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

    v, u, f = t_vuf('nodal', jdmid, ju, lat0)

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
