import datetime
import sys,os
from netCDF4 import Dataset
from ttide.t_tide import t_tide
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf
from ttide import t_utils as tu
import numpy as np
import copy
import matplotlib.dates as mpld
mpld.set_epoch('0000-12-31T00:00:00')
from matplotlib.dates import date2num,num2date
import matplotlib.pyplot as plt
def extract_HC(consfile):
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
    ds = Dataset(consfile)
    X=f.variables['lon_u'][:]
    Y=f.variables['lat_u'][:]

    rloni, rlati = np.meshgrid(X,Y)
    ###
    # Check that the constituents are in the file
    conList = []
    conIDX=[]
    for ncon in range(0,len(f.variables['con'])):
        x=''
        conList.append(''.join([x+n.decode('UTF-8') for n in f.variables['con'][ncon].data]))
        conIDX.append(ncon)

    const = t_getconsts(np.array([]))
    Const= [con.decode('UTF-8') for con in const[0]['name']] 

    consindex = [Const.index(con.ljust(4).upper()) for con in conList]

    tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.   

    var={}
    # interpolating to ROMS grid requires special care with phases!!
    #    this is subjected to parallelization via muitiprocessing.Process - do it!
    #    there is a sandbox in roms repo with a starting exercise

    Vars=['u','v']
    for var0 in Vars:
        var[var0]=np.ones(shape=(len(tfreq),4,rloni.shape[0]*rloni.shape[1]))*-1e-8
        N=-1
        for con,ncon in zip(const[0]['name'][consindex],conIDX):
            N=N+1
            con=con.decode('UTF-8')
            print("extracting %s for %s" %(var0, con))   
            Uim =  f.variables[var0+"Im"][ncon]
            Ure = f.variables[var0+"Re"][ncon]
            h= Ure+Uim*1j
            amp=np.abs(h)
            pha=np.arctan2(-np.imag(h),np.real(h))

            var[var0][N,0,:] = amp.flatten()
            var[var0][N,2,:] = pha.flatten()#*180./np.pi



    return var,tfreq,consindex,rloni.flatten(),rlati.flatten()

def get_tide(ju,freq,tidecon0,t_time,lat0):
    tidecon=copy.deepcopy(tidecon0)
    nodes=tidecon.shape[2]

    t_time=date2num(t_time.astype('datetime64'))

    t_time = t_time.reshape(-1, 1)



    #snr = (tidecon[:, 0,0] / tidecon[:, 1,0]) ** 2
    I = np.arange(0,len(ju)).astype(int)#snr > 2

    
    tidecon = tidecon[I, :]
    ju = np.asarray(ju)[I]
    freq = freq[I]


    ap = np.multiply(tidecon[:, 0] / 2.0,np.exp(-1j * tidecon[:, 2] * np.pi / 180))
    am = np.conj(ap)


    jdmid = t_time[0]#np.mean(t_time[0:((2 * int((max(t_time.shape) - 1) / 2)) + 1)])
    v, u, f = t_vuf('nodal', jdmid, ju, lat0)

    ap = ap * np.kron(np.ones((ap.shape[1],1)),f * np.exp(+1j * 2 * np.pi * (u + v))).T
    am = am * np.kron(np.ones((ap.shape[1],1)),f * np.exp(-1j * 2 * np.pi * (u + v))).T



    n, m = t_time.shape
    yout = np.zeros([ap.shape[1], 1], dtype='complex128')
    touter = np.outer(24 * 1j * 2 * np.pi * freq, t_time[0])

    touter=np.kron(np.ones((1,ap.shape[1])),touter)

    yout[:,0]=np.sum(np.multiply(np.exp(touter), ap), axis=0)+np.sum(np.multiply(np.exp(-touter), am), axis=0)

    tide=np.real(yout)


    return tide
def process(consfile,tstart):



    var,tfreq,consindex,X,Y=extract_HC(consfile)
    TS={}
    for v in var.keys():
        TS[v]=get_tide(consindex,tfreq,var[v],np.array(tstart),np.mean(Y))

    fig, ax = plt.subplots()
    bad=np.isnan( TS['u'][:,0]) | np.isnan(TS['v'][:,0])

    q = ax.quiver(X[~bad], Y[~bad], TS['u'][~bad,0], TS['v'][~bad,0])
    ax.quiverkey(q, X=0.3, Y=1.1, U=10,
                 label='Quiver key, length = 10', labelpos='E')

    plt.show()
if __name__ == "__main__":
    
    # import argparse
    # parser = argparse.ArgumentParser(prog='make_plot_from_cons.py', usage='%(prog)s consfile tstart')
    # parser.add_argument('consfile', type=str,help='Folder with all the SCHSIM files')
    # parser.add_argument('params', type=str,nargs='+',help='Parameters (i.e e um vm')
    # parser.add_argument('tstart', type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),help='start time')


    # args = parser.parse_args()
    

    # process(
    #     args.consfile,\
    #     args.tstart)


    process('/home/remy/Calypso/Projects/MarlboroughSounds/bnd/cookstrait_oceanum_tidecons.nc','2020-01-01')