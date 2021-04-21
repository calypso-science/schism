#!/usr/bin/env python3.7
import copy
import numpy as np
import netCDF4
import numpy.matlib
from interp2D import mask_interp
#from vacumm.misc.grid.regridding import fill2d
#from vcmq import create_time,grid2xy,extend2d

from matplotlib.dates import date2num,num2date
import fiona
from matplotlib.path import Path

class InitialConditions(object):

    def __init__(self,fileout,hgrid,t0,value=None,shapefile=None,ncfile=None,var=None,bnd=None,distance=None,strength=None, logger=None):
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
        elif strength:
            self._create_dist_from_bnd_gr3(distance,bnd,strength)
        else:
            if fileout.endswith('.prop'):
               self._create_constante_prop(value)
            else:
               self._create_constante_gr3(value)

    def _create_dist_from_bnd_gr3(self,distance,bnds,strength):
        self.hgrid.values[:]=0
        lon=self.hgrid.longitude
        lat=self.hgrid.latitude
        
        bnd_nodes=[self.hgrid.boundaries[None][bnd]['indexes'] for bnd in self.hgrid.boundaries[None]]

        node_to_take=[]
        for x in bnds:
            node_to_take+=[int(y)-1 for y in bnd_nodes[x-1]]


        lons=self.hgrid.x
        lats=self.hgrid.y

        dist=np.ones((len(lons),))*np.inf
        for node in node_to_take:
            xs=np.abs(lons-self.hgrid.x[node])
            ys=np.abs(lats-self.hgrid.y[node])
            ds=np.sqrt(xs**2+ys**2)
            dist=np.minimum(ds,dist)


        gd_node=dist<distance


        rnu_max=1./strength/86400.
        import pdb;pdb.set_trace()
        self.hgrid.values[gd_node]=(1-(dist[gd_node]/distance))*rnu_max

        self.hgrid.write(self.fileout)
        self.logger.info("  %s exported"%self.fileout)

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
        data=netCDF4.Dataset(ncfile)
            # Read the grid file

        X=data.variables['longitude'][:]
        Y=data.variables['latitude'][:]
        if len(X.shape)==1:
            xx, yy = np.meshgrid(X, Y)

        lon=self.hgrid.longitude
        lat=self.hgrid.latitude

        

        time0=netCDF4.num2date(data['time'][:],data['time'].units)
        #time0=[np.datetime64(x) for x in time0]

        geo_idx = (np.abs(date2num(time0)-date2num(self.t0))).argmin() # closest timestep
       
        varin=data[var][geo_idx]
        if len(varin.shape)>2:
            varin=varin[0] # get surface level

        src=mask_interp(xx,yy,varin)

        tb=src(np.vstack((lon,lat)).T)

        if np.any(np.isnan(tb)):
            import pdb;pdb.set_trace()

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

