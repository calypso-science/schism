# from bctide import develop_bnd
from openbnd import  OpenBoundaries
import numpy as np
from tidal_tools import extract_HC,get_tide
from matplotlib.dates import num2date,date2num
import copy
# from matplotlib.dates import num2date,date2num
# from tidal_tools import extract_HC,get_tide
from res_tools import get_file_interpolator,vertical_extrapolation
from filetype import create_hotstartNC
# import xarray as xr
from interp2D import mask_interp
import interpvert
# import netCDF4
# import scipy.io
# from scipy.interpolate import griddata,interp1d
# import interpvert
# import time
# import numpy.matlib

INTERNAL_EDGE = 0
BOUNDARY_EDGE = 1
def fill_in_gap(tmp):
    for nlev in range(1,tmp.shape[1]):
        bad=np.isnan(tmp[:,nlev])
        if np.any(bad):
            tmp[bad,nlev]=tmp[bad,nlev-1]
    return tmp

class HotStart(object):
    def __init__(self,filename,config,hgrid,vgrid,t0,z0=0.001,logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.config=config
        self.hgrid=hgrid
        self.vgrid=vgrid
        self.t0=t0

        self.tidal=False
        self.residual=False
        self.i23d=2
        self.ivs=1
        self.z0=z0
        self.lat0=np.mean(self.hgrid.latitude)

        self.llat,self.llon,self.zz=self.get_all_nodes()
        self.zz=np.array(self.zz)

        self._node2elems=None
        #
        ne=len(self.hgrid.mesh.elements)
        self.build_centre_of_elems()
        self.build_edges_from_elems()
        self.xs,self.ys,self.zs=self.build_edgecenters()
        ns=len(self.xs)
        nvrt=int(self.vgrid._nvrt)
        nnodes= len(self.hgrid.mesh.nodes)
        ntracer = len(self.config['tracers'])
        self.logger.info("\tWriting:%s" %filename)
        self.dset,self.nc=create_hotstartNC(filename,ntracer,nnodes,ne,ns,nvrt)


    def build_edgecenters(self):
        """ Build centers of sides
            Returns
            -------
            numpy.array
                list of side centers
        """
        edges = np.array([edge[:2] for edge in self._edges])
        #nodes = self._nodes[edges]
        xs=np.mean(self.hgrid.longitude[edges],axis=1)
        ys=np.mean(self.hgrid.latitude[edges],axis=1)
        zs=np.ndarray((len(ys),self.zz.shape[1]))
        for edgi in range(0,edges.shape[0]):
            z1=self.vgrid.sigma_to_zlayer(edges[edgi,0],self.hgrid.h[edges[edgi,0]]*-1,0.,0.1)
            z2=self.vgrid.sigma_to_zlayer(edges[edgi,1],self.hgrid.h[edges[edgi,1]]*-1,0.,0.1)
            zs[edgi,:]=(np.array(z1)+np.array(z2))/2


        return xs,ys,zs

    def build_elem_balls(self):
        """ Build balls of elements around each node
        """
        if self._node2elems is not None:
            if self.logger is not None:
                self.logger.info("\tRemapping nodes to elements...")
            del self._node2elems
        if self.logger is not None:
            self.logger.info("\tMapping nodes to elements...")
        # build array for point->element lookup
        # Use set for later convenience
        self._node2elems = [set() for _ in range(len(self.hgrid.mesh.nodes))]
        for elem_i in range(len(self.hgrid.mesh.elements)):
            elem = self.hgrid.mesh.elements[elem_i]

            for node_index in elem:
                self._node2elems[int(node_index)-1].add(int(elem_i))

    def get_elems_i_from_node(self, node_i):
        """ Get the ball of elements around a node, node_i
            Parameters()
            ----------
            node_i: int
                node index (zero-based)
            Returns
            -------
            set
                set of element indexes
        """
        if self._node2elems is None:
            self.build_elem_balls()

        return self._node2elems[int(node_i)]
    def build_centre_of_elems(self):
        self.logger.info("\tBuilding center of elements")
        nvrt=int(self.vgrid._nvrt)
        self.xe=np.ndarray((len(self.hgrid.mesh.elements)))
        self.ye=np.ndarray((len(self.hgrid.mesh.elements)))
        self.ze=np.ndarray((len(self.hgrid.mesh.elements),nvrt))
        
        for elem_i, elem in enumerate(self.hgrid.mesh.elements):
            elem=[int(x)-1 for x in elem]
            self.xe[elem_i]=np.mean(self.hgrid.longitude[elem])
            self.ye[elem_i]=np.mean(self.hgrid.latitude[elem])

            ztmp=np.ndarray((3,nvrt))
            for n in range(0,3):
                ztmp[n]=self.vgrid.sigma_to_zlayer(elem[n],self.hgrid.h[elem[n]]*-1,0.,0.1)

            self.ze[elem_i]=np.mean(ztmp,axis=0)


    def build_edges_from_elems(self):
        """ This is a function copied and modified from TriGrid
        """
        self.logger.info("\tBuilding edges from elements")
        # iterate over elements, and for each element, if it's index
        # is smaller than a neighbor or if no neighbor exists,
        # write an edge record
        edges = []

        # this will get built on demand later.
        self._node2edges = None

        for elem_i, elem in enumerate(self.hgrid.mesh.elements):
            # find the neighbors:
            # the first neighbor: need another element that has
            # both self._elems[i,0] and self._elems[i,1] in its
            # list.
            elem=[int(x)-1 for x in elem]
            my_set = set([elem_i])
            for j in range(3):
                node_a = elem[j]
                node_b = elem[(j + 1) % 3]

                elem_ball_node_a = self.get_elems_i_from_node(node_a)
                elem_ball_node_b = self.get_elems_i_from_node(node_b)
                
                # the intersection is us and our neighbor
                # so difference out ourselves...
                adj_elem_of_edge = elem_ball_node_a.intersection(elem_ball_node_b).difference(my_set)
                # and maybe we get a neighbor, maybe not (we're a boundary)
                n_neighbors = len(adj_elem_of_edge)
                if n_neighbors == 1:
                    adj_elem_i = adj_elem_of_edge.pop()
                elif n_neighbors == 0:
                    adj_elem_i = -1
                else:
                    raise RuntimeError("Cannot have more than two neighbors "
                                       "of one edge.")
                if adj_elem_i == -1 or elem_i < adj_elem_i:
                    edges.append((node_a,
                                  node_b,
                                  BOUNDARY_EDGE if adj_elem_i == -1 else INTERNAL_EDGE,
                                  elem_i, adj_elem_i))

        self._edges = np.array(edges, dtype=np.int32)

    def get_all_nodes(self):
        Lat_ocean=[]
        Lon_ocean=[]
        Zlayer_ocean=[] 
        for node_i in range(0,len(self.hgrid.longitude)):
            Lat_ocean.append(self.hgrid.latitude[node_i])
            Lon_ocean.append(self.hgrid.longitude[node_i])
            Zlayer_ocean.append(self.vgrid.sigma_to_zlayer(node_i,self.hgrid.h[node_i]*-1,0.,0.1))

        return Lat_ocean,Lon_ocean,Zlayer_ocean 

    def get_elev(self,config):
        Time=date2num(self.t0)
        if 'tidal' in config:
            self.cons=config['tidal'].get('cons',None)
            OpenBoundaries.add_tide(self,config['tidal'])

        if 'residual' in config:
            OpenBoundaries.add_res(self,config['residual'])

        total=np.zeros(shape=(len(self.llon),1))

                    # get tide
        if self.tidal:
            var=list(self.HC.keys())[0]
            # horizontal interpolation
            tmp=get_tide(self.constidx,self.tfreq,self.HC[var],np.array(Time),np.mean(self.llat))                      
            total+=tmp
        
        if self.residual:
            if 'longitude' in self.res_file.dims:
                lon_name='longitude'
                lat_name='latitude'
            else:
                lon_name='lon'
                lat_name='lat'                

            xx,yy=np.meshgrid(self.res_file[lon_name][:],self.res_file[lat_name][:])

            var=self.res_vars[0]
            arri=self.res_file[var][:]
            arri_time=arri.interp(time=num2date(date2num(np.datetime64(num2date(Time)))).strftime('%Y-%m-%d %H:%M:%S'))
            arr=mask_interp(xx,yy,arri_time.to_masked_array())
            tb=arr(np.vstack((self.llon,self.llat)).T, nnear=6, p=2)
            total[:,0]=total[:,0]+tb

        return total
    def get_uv(self,config):
        Time=date2num(self.t0)
        self.llon=copy.deepcopy(self.xs)
        self.llat=copy.deepcopy(self.ys)
        Nlev=self.zz.shape[1]
        total=np.zeros(shape=(2,len(self.llon),Nlev))

        Time=date2num(self.t0)
        if 'tidal' in config:
            self.cons=config['tidal'].get('cons',None)
            OpenBoundaries.add_tide(self,config['tidal'])

        if 'residual' in config:
            OpenBoundaries.add_res(self,config['residual'])

        # get tide
        if self.tidal:
            var=self.HC.keys()
            for i,v in enumerate(sorted(var)):
                # horizontal interpolation
                tmp=get_tide(self.constidx,self.tfreq,self.HC[v],np.array(Time),np.mean(self.llat))      
                tmp=vertical_extrapolation(tmp,self.zs,z0=self.z0)   
                tpm=fill_in_gap(tmp)                   
                total[i,:,:]=total[i,:,:]+tmp

        if self.residual:
            if 'longitude' in self.res_file.dims:
                lon_name='longitude'
                lat_name='latitude'
            else:
                lon_name='lon'
                lat_name='lat'                

            xx,yy=np.meshgrid(self.res_file[lon_name][:],self.res_file[lat_name][:])
            var=self.res_vars

            zi=self.res_file['lev'][:].values
            if np.mean(zi)>0:
              zi=zi*-1

            for i,v in enumerate(sorted(var)):
                arri=self.res_file[v][:]
                arri_time=arri.interp(time=num2date(date2num(np.datetime64(num2date(Time)))).strftime('%Y-%m-%d %H:%M:%S'))
                
                tb=np.ndarray((len(self.llon),Nlev))
                tmp=np.ndarray((len(self.llon),arri_time.shape[0]))*np.nan
                for nlev in range(0,arri_time.shape[0]):
                    if np.any(arri_time[nlev].to_masked_array()):
                        arr=mask_interp(xx,yy,arri_time[nlev].to_masked_array())
                        if len(arr.z)>1:
                            tmp[:,nlev]=arr(np.vstack((self.llon,self.llat)).T, nnear=1, p=2)

                tpm=fill_in_gap(tmp) 

                if self.zs.shape[1]==2: # 2D
                    for p in range(0,tmp.shape[0]):
                        total_depth=self.zs[p,0]
                        bad=np.isnan(tmp[p,:])
                        depth=zi[~bad]
                        vel=tmp[p,~bad]
                        depth=np.insert(depth,0,0,axis=0)
                        ve=0
                        tot=0
                        for dep in range(0,len(depth)-1):
                            dz=depth[dep]-depth[dep+1]
                            dx=vel[dep]
                            ve+=dx*dz
                            tot+=dz
                        
                        tb[p,:]=ve/np.abs(tot)
                else: # 3D
                    varin=np.fliplr(tmp)#[:,~bad])
                    zin=np.fliplr(np.matlib.repmat(zi,varin.shape[0],1))

                    for nl in range(0,self.zs.shape[1]):
                        zout=self.zs[:,nl]
                        tb[:,nl]=np.squeeze(interpvert.interpz1d(varin,zin,zout,np=varin.shape[0],nzin=varin.shape[1],nzout=1,kz=1, null_value=-9.9e15))

                total[i,:,:]=total[i,:,:]+tb
        return total
    def get_tracer_at_elem(self,config):
        Time=date2num(self.t0)
        self.llon=copy.deepcopy(self.xe)
        self.llat=copy.deepcopy(self.ye)
        self.zz=copy.deepcopy(self.ze)
        Nlev=self.zz.shape[1]
        total=np.zeros(shape=(len(self.llon),Nlev))

        Time=date2num(self.t0)
        OpenBoundaries.add_res(self,config)

        if 'longitude' in self.res_file.dims:
            lon_name='longitude'
            lat_name='latitude'
        else:
            lon_name='lon'
            lat_name='lat'                

        xx,yy=np.meshgrid(self.res_file[lon_name][:],self.res_file[lat_name][:])

        zi=self.res_file['lev'][:].values
        if np.mean(zi)>0:
          zi=zi*-1

        var=self.res_vars[0]
        arri=self.res_file[var][:]
        arri_time=arri.interp(time=num2date(date2num(np.datetime64(num2date(Time)))).strftime('%Y-%m-%d %H:%M:%S'))
        
        tb=np.ndarray((len(self.llon),Nlev))
        tmp=np.ndarray((len(self.llon),arri_time.shape[0]))*np.nan
        for nlev in range(0,arri_time.shape[0]):
            if np.any(arri_time[nlev].to_masked_array()):
                arr=mask_interp(xx,yy,arri_time[nlev].to_masked_array())
                if len(arr.z)>1:
                    tmp[:,nlev]=arr(np.vstack((self.llon,self.llat)).T, nnear=1, p=2)

        tpm=fill_in_gap(tmp) 

        if self.zs.shape[1]==2: # 2D
            for p in range(0,tmp.shape[0]):
                total_depth=self.zz[p,0]
                bad=np.isnan(tmp[p,:])
                depth=zi[~bad]
                vel=tmp[p,~bad]
                depth=np.insert(depth,0,0,axis=0)
                ve=0
                tot=0
                for dep in range(0,len(depth)-1):
                    dz=depth[dep]-depth[dep+1]
                    dx=vel[dep]
                    ve+=dx*dz
                    tot+=dz
                
                tb[p,:]=ve/np.abs(tot)
        else: # 3D
            varin=np.fliplr(tmp)#[:,~bad])
            zin=np.fliplr(np.matlib.repmat(zi,varin.shape[0],1))

            for nl in range(0,self.zz.shape[1]):
                zout=self.zz[:,nl]
                tb[:,nl]=np.squeeze(interpvert.interpz1d(varin,zin,zout,np=varin.shape[0],nzin=varin.shape[1],nzout=1,kz=1, null_value=-9.9e15))


        return tb

    def get_tracet_at_node(self,config):
        Time=date2num(self.t0)
        self.llat,self.llon,self.zz=self.get_all_nodes()
        self.zz=np.array(self.zz)
        Nlev=self.zz.shape[1]
        total=np.zeros(shape=(len(self.llon),Nlev))

        Time=date2num(self.t0)
        OpenBoundaries.add_res(self,config)

        if 'longitude' in self.res_file.dims:
            lon_name='longitude'
            lat_name='latitude'
        else:
            lon_name='lon'
            lat_name='lat'                

        xx,yy=np.meshgrid(self.res_file[lon_name][:],self.res_file[lat_name][:])

        zi=self.res_file['lev'][:].values
        if np.mean(zi)>0:
          zi=zi*-1

        var=self.res_vars[0]
        arri=self.res_file[var][:]
        arri_time=arri.interp(time=num2date(date2num(np.datetime64(num2date(Time)))).strftime('%Y-%m-%d %H:%M:%S'))
        
        tb=np.ndarray((len(self.llon),Nlev))
        tmp=np.ndarray((len(self.llon),arri_time.shape[0]))*np.nan
        for nlev in range(0,arri_time.shape[0]):
            if np.any(arri_time[nlev].to_masked_array()):
                arr=mask_interp(xx,yy,arri_time[nlev].to_masked_array())
                if len(arr.z)>1:
                    tmp[:,nlev]=arr(np.vstack((self.llon,self.llat)).T, nnear=1, p=2)

        tpm=fill_in_gap(tmp) 

        if self.zs.shape[1]==2: # 2D
            for p in range(0,tmp.shape[0]):
                total_depth=self.zz[p,0]
                bad=np.isnan(tmp[p,:])
                depth=zi[~bad]
                vel=tmp[p,~bad]
                depth=np.insert(depth,0,0,axis=0)
                ve=0
                tot=0
                for dep in range(0,len(depth)-1):
                    dz=depth[dep]-depth[dep+1]
                    dx=vel[dep]
                    ve+=dx*dz
                    tot+=dz
                
                tb[p,:]=ve/np.abs(tot)
        else: # 3D
            varin=np.fliplr(tmp)#[:,~bad])
            zin=np.fliplr(np.matlib.repmat(zi,varin.shape[0],1))

            for nl in range(0,self.zz.shape[1]):
                zout=self.zz[:,nl]
                tb[:,nl]=np.squeeze(interpvert.interpz1d(varin,zin,zout,np=varin.shape[0],nzin=varin.shape[1],nzout=1,kz=1, null_value=-9.9e15))


        return tb

    def set_hotstart(self):
        
        # get ssh at nodes
        if 'ssh' in self.config:
            self.logger.info("\Adding elev to hostart")
            ssh=self.get_elev(self.config['ssh'])
            self.dset['eta2'][:]=ssh

        # get U and V at the side
        if 'uv' in self.config:
            self.logger.info("\Adding UV to hostart")
            uv=self.get_uv(self.config['uv'])
            self.dset['su2'][:]=uv[0]
            self.dset['sv2'][:]=uv[1]


        for key in self.config['tracers']:
            self.logger.info("\Adding %s to hostart" % key)
            tr_nd=self.get_tracet_at_node(self.config['tracers'][key])
            self.dset['tr_nd'][:,:,key-1]=tr_nd
            self.dset['tr_nd0'][:,:,key-1]=tr_nd
            tr_el=self.get_tracer_at_elem(self.config['tracers'][key])
            self.dset['tr_el'][:,:,key-1]=tr_el

        self.nc.close()
