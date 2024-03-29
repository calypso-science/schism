#!/usr/bin/env python3.7
import os
import numpy as np
from matplotlib.dates import date2num,num2date
import xarray as xr
#import cdms2
#from vacumm.misc.grid.regridding import fill2d
from scipy.interpolate import interp1d
from netCDF4 import Dataset
from atmospheric import rh2sh

file_sections={}
file_sections['air']=['uwind','vwind','prmsl','stmp','spfh']
file_sections['prc']=['prate']
file_sections['rad']=['dlwrf','dswrf']


def create_netcdf_file(filename,lon,lat,T,file_sections):

        T=T-1
        root_grp = Dataset(filename, 'w', format='NETCDF4',clobber=True)
        root_grp.Conventions = 'CF-1.0'
        # dimensions
        root_grp.createDimension('time', len(T))
        root_grp.createDimension('nx_grid', lon.shape[1])
        root_grp.createDimension('ny_grid', lon.shape[0])
        # variables
        times = root_grp.createVariable('time', 'f8', ('time',))
        base_date=np.double(num2date(T[0]).year),np.double(num2date(T[0]).month),np.double(num2date(T[0]).day),np.double(num2date(T[0]).hour),np.double(num2date(T[0]).minute),np.double(num2date(T[0]).second)
        times.base_date=base_date
        latitudes = root_grp.createVariable('lat', 'f4', ('ny_grid','nx_grid',))
        longitudes = root_grp.createVariable('lon', 'f4', ('ny_grid','nx_grid',))
        ## add the main var
        times[:]=np.ceil(abs(T-T[0])*1000)/1000
        latitudes[:,:]=lat
        longitudes[:,:]=lon
        ## save the rest
        temp={}
        for m in file_sections:
        	temp[m] = root_grp.createVariable(m, 'f4', ('time', 'ny_grid','nx_grid',))
        return temp,root_grp

class Meteo(object):
    """A class that manages atmospheric forcing fields
    """

    def __init__(self,in_dir,t0,t1,dt, logger=None):
        """ Constructor"""

        if logger:
            self.logger = logger
            self.logger.info("Processing atmospheric forcing")


        self._initialise(in_dir)
        self._make_sflux_txt()
        self.Tforcing=np.arange(date2num(t0),date2num(t1)+1,dt/(24.))

        self.dataset=[]
        self.lon=[]
        self.lat=[]

        self.rh2m=False



    def _initialise(self,input_dir):
        self.input_dir=input_dir
        if not os.path.exists(input_dir):
            os.makedirs(input_dir)    	

    def _make_sflux_txt(self):

        tmpl=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','template','sflux_inputs.tmpl')
        os.system('cp '+tmpl +' '+os.path.join(self.input_dir,os.path.expanduser('sflux_inputs.txt')))


    def add_dataset(self,filename,vardict):
        data=xr.open_dataset(filename)
        if 'longitude' in data.dims:
            lon_name='longitude'
            lat_name='latitude'
        else:
            lon_name='lon'
            lat_name='lat'   

        var={}
        for key in vardict.keys():
            if vardict[key] in data:
                arri=data[vardict[key]][:]
                var[key] =arri.interpolate_na(dim='time',method='nearest')
                if 'rh2m' in vardict[key]:
                    self.rh2m=True
            else:
                var[key]=vardict[key]




        self.dataset.append(var)

        if data[lat_name][0]>data[lat_name][-1]:
            self.need_sorting=True
            [lon,lat]=np.meshgrid(data[lon_name][:],np.sort(data[lat_name][:]))
        else:
            self.need_sorting=False
            [lon,lat]=np.meshgrid(data[lon_name][:],np.sort(data[lat_name][:]))

        self.lon.append(lon)
        self.lat.append(lat)





    def make_meteo(self):
    ## get time series daily
        dt=self.Tforcing[2]-self.Tforcing[1]
        unique_days=np.unique(np.floor(self.Tforcing))


        for section in file_sections:

            for k in range(0,len(unique_days)):
                t=np.arange(unique_days[k],unique_days[k]+1-dt,dt)
                tin=[np.datetime64(num2date(x)) for x in t]

               
                for n,dataset in enumerate(self.dataset):
                   
                    netcdf_name=('sflux_%s_%.f.%04.f.nc' % (section,n+1,k+1))
                   

                    temp,root_grp=create_netcdf_file(os.path.join(self.input_dir,netcdf_name),self.lon[n],self.lat[n],t+1,file_sections[section])
                   
                    for var in file_sections[section]:
                        if hasattr(dataset[var],'interp'):
                            tmp=dataset[var].interp(time=tin)
                            if self.rh2m and var is 'prmsl':
                              prmsl=tmp
                            if self.rh2m and var is 'stmp':
                           	  stmp=tmp
                            if self.rh2m and var is 'spfh':
                              tmp=rh2sh(tmp/100.,prmsl,stmp)
                            
                            if self.need_sorting:
                                tmp=tmp[:,::-1,:]

                            for nn in range(0,tmp.shape[0]):
                                temp[var][nn,:,:]=tmp[nn,:,:]
                        else:
                            temp[var][nn,:,:]=dataset[var]


                    root_grp.close()
                    str='ncks -O --fl_fmt=classic %s %s' %(os.path.join(self.input_dir,netcdf_name),os.path.join(self.input_dir,netcdf_name))
                    os.system(str)
                


