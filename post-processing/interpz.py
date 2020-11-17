#!/usr/bin/python
from argparse import ArgumentParser
import numpy
import netCDF4
import datetime
import interpz

null_real = -9e15
kz = 1

parser = ArgumentParser(description='convert sigma residuals to z residuals')
parser.add_argument('-fz', '--filez',  required = True, help = 'file with zcor [../outncsig/s1_sig.nc]')
parser.add_argument('-fin','--filein', required = True, help = 'residuals file in sigma [../compute_res/res_s1_sig.nc]')
parser.add_argument('-fout','--fileout', required = True, help = 'residuals file in z [res_s1_z.nc]')
parser.add_argument('-d','--deplevels', default='0, 5, 10, 15, 20, 25, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000', help='csv list of depth levers') 
args = vars(parser.parse_args())

filez      = args['filez']
filein     = args['filein']
fileout    = args['fileout']
deplevels  = map(int, args['deplevels'].split(','))
nz         = len(deplevels)

print 'filez     = ', filez
print 'filein    = ', filein
print 'fileout   = ', fileout
print 'deplevels = ', deplevels
print 'nz = ', nz

zout = -1.*numpy.array(deplevels[::-1])
print 'zout = ', zout

## open filein ##
print 'Opening '+filein
f = netCDF4.Dataset(filein)
nt = len(f.variables['time'])
print 'nt     = ', nt

sigma = f.variables['sigma'][:]
nsig = len(sigma)
print 'sigma = ', sigma
print 'nsig = ', nsig

nnodes = f.variables['uo'].shape[2]
print 'nnodes = ', nnodes

## open filez ##
print 'Opening '+filez
f1 = netCDF4.Dataset(filez)
# check z variable is conformal with filein
tmp = f1.variables['z'].shape
print tmp
if tmp[0] != nt:
    print "filez = %s doesn't have same nt as filein = %s : %d != %d" % (filez, filein, tmp[0], nt)
    sys.exit(-2)

if tmp[1] != nsig:
    print "filez = %s doesn't have same nsig as filein = %s: %d != %d" % (filez, filein, tmp[1], nsig)
    sys.exit(-3)

if tmp[2] != nnodes:
    print "filez = %s doesn't have same nnodes as filein = %s: %d != %d" % (filez, filein, tmp[2], nnodes)
    sys.exit(-3)

####################
## create fileout ##
####################

print 'Creating '+fileout
fout = netCDF4.Dataset(fileout, 'w', format = 'NETCDF3_64BIT')
fout.createDimension('node', nnodes)
fout.createDimension('time', nt)
fout.createDimension('lev', nz)

print 'writing time '
x = fout.createVariable('time', 'float', 'time')
x.long_name = f.variables['time'].long_name
x.units     = f.variables['time'].units
x[:] = f.variables['time'][:]

print 'writing lev '
x = fout.createVariable('lev', 'float', 'lev')
x.long_name     = 'Level below sea surface'
x.units         = 'm'
x.positive      = 'down'
x[:] = deplevels

fout.sync()

 ## copying eo so that we don't have to use filein (with eo) anymore
misse = f.variables['eo']._FillValue
eo = fout.createVariable('eo', 'int16', ('time', 'node'), fill_value = misse)
eo.scale_factor = f.variables['eo'].scale_factor
eo.add_offset   = f.variables['eo'].add_offset
eo.long_name    = f.variables['eo'].long_name
eo.relative_to  = f.variables['eo'].relative_to
eo.units        = f.variables['eo'].units
eo.valid_range  = f.variables['eo'].valid_range
eo[:] = f.variables['eo'][:]
fout.sync()

##########################
###    interpolation   ###
##########################
## loop u and v together so that we dont have to read z variable twice

print 'creating uo(time,lev,nodes) at ', datetime.datetime.now()
missu = f.variables['uo']._FillValue
missv = f.variables['uo']._FillValue
print 'missu, missv = ', missu, missv

uo = fout.createVariable('uo', 'int16', ('time', 'lev', 'node'), fill_value = missu)
uo.scale_factor = f.variables['uo'].scale_factor
uo.add_offset   = f.variables['uo'].add_offset
uo.long_name    = f.variables['uo'].long_name
uo.units        = f.variables['uo'].units
uo.valid_range  = f.variables['uo'].valid_range

vo = fout.createVariable('vo', 'int16', ('time', 'lev', 'node'), fill_value = missv)
vo.scale_factor = f.variables['vo'].scale_factor
vo.add_offset   = f.variables['vo'].add_offset
vo.long_name    = f.variables['vo'].long_name
vo.units        = f.variables['vo'].units
vo.valid_range  = f.variables['vo'].valid_range

######################
# tide mask doesn't go into interp but its ok because if the node is maksed, it will not be used.
# tidal masked nodes have velocities of zero
# but nodes interpolated to below bathy must be masked with null_real

for it in range(nt):
    zin   = f1.variables['z'][it,:,:]

    varin = f.variables['uo'][it,:,:]
    # print 'varin = ', varin[:,inode]
    # print 'varin.data = ', varin[:,inode].data
    varout = interpz.interpz1d(varin,zin,zout,kz,null_real,np=nnodes,nzin=nsig,nzout=nz)
    varout = varout[::-1,:]
    varout = numpy.ma.array(varout)
    varout = numpy.ma.masked_where(varout == null_real,varout)  # mask below bathy
    # print 'varout = ', varout[:,inode]
    # print 'varout.data = ', varout[:,inode].data
    uo[it,:,:] = varout[:,:]

    varin = f.variables['vo'][it,:,:]
    # print 'varin = ', varin[:,inode]
    # print 'varin.data = ', varin[:,inode].data
    varout = interpz.interpz1d(varin,zin,zout,kz,null_real,np=nnodes,nzin=nsig,nzout=nz)
    varout = varout[::-1,:]
    varout = numpy.ma.array(varout)
    varout = numpy.ma.masked_where(varout == null_real,varout)  # mask below bathy
    # print 'varout = ', varout[:,inode]
    # print 'varout.data = ', varout[:,inode].data
    vo[it,:,:] = varout[:,:]

fout.sync()


print 'Closing  at ', datetime.datetime.now()
f1.close
f.close()
fout.close()
print 'Finished at ', datetime.datetime.now()



