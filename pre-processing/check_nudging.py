#!/usr/bin/env python2.7

import netCDF4
import sys,os
from pyschism.mesh import Hgrid
import matplotlib.pyplot as plt
import numpy

def build_edgecenters(gr):
    """ Build centers of sides
        Returns
        -------
        numpy.array
            list of side centers
    """
    edges = np.array([edge[:2] for edge in gr._edges])
    #nodes = self._nodes[edges]
    xs=np.mean(gr.hgrid.longitude[edges],axis=1)
    ys=np.mean(gr.hgrid.latitude[edges],axis=1)

    return xs,ys

def build_elem_balls(gr):
    """ Build balls of elements around each node
    """
    if gr._node2elems is not None:
        del gr._node2elems

    # build array for point->element lookup
    # Use set for later convenience
    gr._node2elems = [set() for _ in range(len(gr.hgrid.mesh.nodes))]
    for elem_i in range(len(gr.hgrid.mesh.elements)):
        elem = gr.hgrid.mesh.elements[elem_i]

        for node_index in elem:
            gr._node2elems[int(node_index)-1].add(int(elem_i))
    return gr
def get_elems_i_from_node(gr, node_i):
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
    if gr._node2elems is None:
        gr=build_elem_balls(gr)

    return gr._node2elems[int(node_i)]
def build_centre_of_elems(self):

    gr.xe=np.ndarray((len(gr.hgrid.mesh.elements)))
    gr.ye=np.ndarray((len(gr.hgrid.mesh.elements)))

    
    for elem_i, elem in enumerate(gr.hgrid.mesh.elements):
        elem=[int(x)-1 for x in elem]
        gr.xe[elem_i]=np.mean(gr.hgrid.longitude[elem])
        gr.ye[elem_i]=np.mean(gr.hgrid.latitude[elem])


    return gr


def build_edges_from_elems(gr):
    """ This is a function copied and modified from TriGrid
    """
    # iterate over elements, and for each element, if it's index
    # is smaller than a neighbor or if no neighbor exists,
    # write an edge record
    edges = []

    # this will get built on demand later.
    gr._node2edges = None

    for elem_i, elem in enumerate(gr.elements):
        # find the neighbors:
        # the first neighbor: need another element that has
        # both self._elems[i,0] and self._elems[i,1] in its
        # list.
        elem=[int(x)-1 for x in elem]
        my_set = set([elem_i])
        for j in range(3):
            node_a = elem[j]
            node_b = elem[(j + 1) % 3]

            elem_ball_node_a = gr.get_elems_i_from_node(node_a)
            elem_ball_node_b = gr.get_elems_i_from_node(node_b)
            
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

    gr._edges = np.array(edges, dtype=np.int32)
    return gr
def process(filein,hgrid_file,var,epsg,T,L):

    nd,nc,v=read_nudge(filein,var)
    gr=Hgrid.open(hgrid_file,crs="EPSG:%i" % epsg)

    #gr = load_g3(os.path.expanduser())
    Easting= gr.x
    Northing= gr.y
    elems=gr.elements.elements

    Face=numpy.ndarray((len(elems),3))
    for elem_i,elem in enumerate(elems.keys()):
        Face[elem_i]=[int(x)-1 for x in elems[elem][:-1]]


    # gr=build_edges_from_elems(gr)
    # gr.xs,gr.ys=build_edgecenters(gr)

    # Side= numpy.array([edge[:2] for edge in gr.mesh._edges])
    Z=v[0,:,L,-1]



    plot2D(Easting,Northing,Face,Z,var)



def read_nudge(filename,var):
    nc=netCDF4.Dataset(filename)
    nd=len(nc.dimensions['node'])
    v=nc.variables[var]
    return nd,nc,v

def plot2D(X,Y,F,Z,params):
    fig = plt.figure(figsize=(15,9))
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    tt1 = ax.text(.5, 1.05, '', transform = ax.transAxes, va='center',fontsize = 30)
    ax.tick_params(labelsize=15)

    Zmin=Z[Z>0].min()
    Zmax=Z.max()
    levels = numpy.linspace(Zmin,Zmax, 60)
    F=plt.tricontourf(X,Y,F,Z,vmin=Z.min(),vmax=Z.max(),cmap=plt.cm.Spectral_r,levels=levels)
    plt.clim(Zmin,Zmax)
    plt.xlabel('Easting (meters)',fontsize = 30)
    plt.ylabel('Northing (meters)',fontsize = 30)

    cbar=plt.colorbar(F,ax=ax)#numpy.linspace(Zmin,Zmax,10))
    cbar.set_label(r"%s" %(params), size=15)
    cbar.ax.tick_params(labelsize=10) 
    plt.draw()

    ax.set_xlim([X.min(), X.max()])
    ax.set_ylim([Y.min(), Y.max()])

    plt.show()



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog='check_hotstart.py', usage='%(prog)s filein hgrid var (i.e python check_hotstart.py hotstart.nc hgrid.gr3 tr_nd0 -tracer 1 -layer 1')
    ## main arguments
    
    parser.add_argument('filein', type=str,help='name of the input hot file')
    parser.add_argument('hgrid', type=str,help='name of the input grid file')
    parser.add_argument('var', type=str,help='name of the input variable to check')
    parser.add_argument('epsg', type=int,help='epsg')
    parser.add_argument('-tracer', type=int,help='tracer to plot (optional)',default=0)
    parser.add_argument('-layer', type=int,help='layer to plot (optional)',default=0)
    args = parser.parse_args()
    

    process(args.filein,args.hgrid,args.var,args.epsg,args.tracer,args.layer)

