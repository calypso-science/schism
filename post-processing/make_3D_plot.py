import netCDF4
from tvtk.api import tvtk
from mayavi import mlab
import numpy
from scipy.interpolate import griddata
from numpy import mgrid, empty, sin, pi, array, meshgrid, arange, prod
import rasterio
import os
	# some path to your local netcdf file
path = '/home/remy/Calypso/Projects/006_Southport/hydro/child/hg6_dt30_newD/outputs/schout_14.nc'

lim=[1241352,1245886,4825135,4830428]
quiver_res=100
# try and import the dataset, prefer netcdf4, you might want to use pydap also here if you access a dap server.

ds = netCDF4.Dataset(path)
X=ds.variables['SCHISM_hgrid_node_x'][:]
Y=ds.variables['SCHISM_hgrid_node_y'][:]
XY=numpy.array((X,Y)).T
 # read the variables (these are just links to the objects for now)
u = ds.variables['hvel'][0,:,-1,0] # velocity in u direction
v = ds.variables['hvel'][0,:,-1,1] # v direction

Xreg, Yreg, Zreg = numpy.meshgrid(numpy.linspace(lim[0],lim[1], quiver_res), numpy.linspace(lim[2], lim[3], quiver_res),0)

U = griddata(XY, u,(Xreg,Yreg),method='linear')
V = griddata(XY, v,(Xreg,Yreg),method='linear')

mag=numpy.sqrt(U**2+V**2)

sg = tvtk.StructuredGrid(dimensions=Xreg.shape)

# The actual points.
pts = empty(Xreg.shape + (3,), dtype=float)
pts[..., 0] = Xreg
pts[..., 1] = Yreg
pts[..., 2] = Zreg

pts = pts.transpose(2, 1, 0, 3).copy()
pts.shape = pts.size // 3, 3

sg.points = pts

vectors = empty(Xreg.shape + (3,), dtype=float)
vectors[..., 0] = U
vectors[..., 1] = V
vectors[..., 2] = 0


# We reorder the points, scalars and vectors so this is as per VTK's
# requirement of x first, y next and z last.


vectors = vectors.transpose(2, 1, 0, 3).copy()
vectors.shape = vectors.size // 3, 3

# Create the dataset.
sg.point_data.vectors = vectors
sg.point_data.vectors.name = 'velocity'


# Now visualize the data.
d = mlab.pipeline.add_dataset(sg)
gx = mlab.pipeline.grid_plane(d)
gx.actor.actor.visibility = False
gy = mlab.pipeline.grid_plane(d)
gy.actor.actor.visibility = False
gy.grid_plane.axis = 'y'
gz = mlab.pipeline.grid_plane(d)
gz.grid_plane.axis = 'z'
gz.actor.actor.visibility = False
#iso = mlab.pipeline.iso_surface(d)
#iso.contour.maximum_contour = 75.0
vec = mlab.pipeline.vectors(d)
vec.glyph.mask_input_points = True
vec.glyph.glyph.range = array([0. , 1.5])
vec.glyph.glyph.scale_factor = 300.0


d.scene.z_plus_view()
# file='/home/remy/Calypso/Projects/006_Southport/hydro/child/animation/chart_nztm.jpg'
# #img = rasterio.open(file)
# #lidar_dem_im = img.read(1, masked=True)

# bmp1 = tvtk.JPEGReader()
# bmp1.file_name=file #any jpeg file


# Ixlim=[1241221.9970655012875795,1249403.4132295222952962]
# Iylim=[4821406.8234612224623561,4832765.3438056204468012]
# Ximage_res=2267#6582
# Yimage_res=2647# 9138

# Ximg, Yimg = numpy.meshgrid(numpy.linspace(Ixlim[0],Ixlim[1], Ximage_res), numpy.linspace(Iylim[0], Iylim[1], Yimage_res))
# #texture = tvtk.Texture()
# #texture.interpolate=0
# #texture.set_input_data(0,bmp1.get_output())
# texture = tvtk.Texture(input_connection=bmp1.output_port, interpolate=0)


# data=numpy.random.random((Ximage_res,Yimage_res))*0

# surf = mlab.surf(Ximg,Yimg,data,color=(1,1,1))#,warp_scale=0.2) 
# surf.actor.enable_texture = True  
# surf.actor.tcoord_generator_mode = 'plane'  
# surf.actor.actor.texture = texture
fig=mlab.figure(size = (1024,768),\
            bgcolor = (1,1,1), fgcolor = (0.5, 0.5, 0.5))
import mplleaflet

gj = mplleaflet.fig_to_geojson(fig=fig)