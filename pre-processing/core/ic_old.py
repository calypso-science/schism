#!/usr/bin/env python3.7
import copy
import numpy as np
import cdms2

from vacumm.misc.grid.regridding import fill2d
from vcmq import create_time,grid2xy,extend2d

from matplotlib.dates import date2num,num2date
import fiona
from matplotlib.path import Path

class InitialConditions(object):

    def __init__(self,fileout,hgrid,t0,value=None,shapefile=None,ncfile=None,var=None, logger=None):
        """ Constructor"""

        if logger:
            self.logger = logger
            self.logger.info("Processing initial conditions")


        self.hgrid=copy.deepcopy(hgrid.mesh)
        self.hgrid.values[:]=np.ones(self.hgrid.x.shape,dtype=np.float64)*0
        self.hgrid.longitude=hgrid.longitude
        self.hgrid.latitude=hgrid.latitude
        self.t0=t0
        self.fileout=fileout

            

        if ncfile:
            self._create_nc_gr3(ncfile,var)
        elif shapefile:
            self._create_shp_gr3(shapefile,value)
        else:
            if fileout.endswith('.prop'):
               self._create_constante_prop(value)
            else:
               self._create_constante_gr3(value)


    def _create_constante_prop(self,value):

        f = open(self.fileout, 'w')
        n_elems = len(self.hgrid.elements)
        attr_array=np.ones((n_elems,))
        attr_array[:]=value

        for i in range(n_elems):
            buf = "%d %d\n" % (i+1, attr_array[i])
            f.write(buf)
        f.close()

        self.logger.info("	%s exported"%self.fileout)
    
    def _create_constante_gr3(self,value):
        self.hgrid.values[:]=value*-1.
        self.hgrid.write(self.fileout)
        self.logger.info("	%s exported"%self.fileout)

    def _create_shp_gr3(self,shapefile,value):
        shapes=shapefile.split(',')
        self.hgrid.values[:]=value*-1.
        for i,file in enumerate(shapes):
            self._add_value_from_shapefile(file)

        self.hgrid.write(self.fileout)
        self.logger.info("  %s exported"%self.fileout)
    def _create_nc_gr3(self,ncfile,var):
        data=cdms2.open(ncfile)
        
        lon=self.hgrid.longitude
        lat=self.hgrid.latitude
        src=fill2d(data[var][:], method='carg')
        time0=[t.torelative('days since 1-1-1').value for t in data['time'].asRelativeTime()]
        tin=create_time(np.ones(len(lon))*date2num(self.t0)+1,units='days since 1-1-1')
        tb=grid2xy(src,xo=lon,yo=lat,method='linear',to=tin)

        if np.any(tb.mask==True):
            bad=(tb.mask==True).nonzero()[0]
            tin_bad=create_time(np.ones(len(bad))*date2num(self.t0)+1,units='days since 1-1-1')
            tb[bad]=grid2xy(src,xo=np.array(lon)[bad].tolist(),yo=np.array(lat)[bad].tolist(),method='nearest',to=tin_bad)

        self._create_constante_gr3(tb)



    def _add_value_from_shapefile(self,filename):
        shape = fiona.open(filename)
        for poly in shape:
            if poly['geometry']!=None :
                    nvert=len(poly['geometry']['coordinates'][0])
                    vert=[(poly['geometry']['coordinates'][0][x][0],poly['geometry']['coordinates'][0][x][1]) for x in range(0,nvert)]

                    Vin=poly['properties']['in']

                    

                    p=Path(vert)

                    flag = p.contains_points(self.hgrid.xy)
                    self.hgrid.values[flag]=Vin*-1.


                    if 'out' in poly['properties']:
                        Vout=poly['properties']['out']
                        if Vout:
                            self.hgrid.values[~flag]=Vout*-1.

