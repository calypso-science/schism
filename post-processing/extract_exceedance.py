	#!/usr/bin/env python2.7

import os,sys
from glob import glob
import netCDF4
import numpy as np
import copy
	
					
						
def process_exce(ncout,all_files,params,threshs):
	

	for param in params:	
		for i,file in enumerate(all_files):
			print '%s => %s' % (param,file)
			ncin = netCDF4.Dataset(file,  'r')	


			D=ncin.variables[param][:]
			ncin.close()
			little_divier=np.sum(D.mask==False,0)
			if type(little_divier)==type(np.int64()):little_divier=D.shape[0]
			if i==0:
				divider=np.ones(shape=D.shape[1:])*0
				matrix=np.ones(shape=(len(threshs),)+D.shape[1:])*0
			divider=divider+little_divier
					
			for ithresh,thresh in enumerate(threshs):
				import pdb;pdb.set_trace()
				matrix[ithresh,:]=(D>thresh).sum(axis=0)


			divider[divider==0]=np.nan
		
			
			ncout.variables[param+'_'+stat][:]=np.squeeze(matrix[istat,:]*100./divider)

	
	return ncout

def create_output_file(filout,first_file,threshs,params):

	ncin = netCDF4.Dataset(first_file,  'r')
 	y = ncin.variables['SCHISM_hgrid_node_y'][:]
	x= ncin.variables['SCHISM_hgrid_node_x'][:]
	nnodes=len(x)
	ele=ncin.variables['SCHISM_hgrid_face_nodes'][:]
	sigma=ncin.dimensions['nSCHISM_vgrid_layers'].size

	ncout = netCDF4.Dataset(filout,  'w')
        
	dnode = ncout.createDimension('nSCHISM_hgrid_node', len(x))
	dele = ncout.createDimension('nSCHISM_hgrid_face', len(ele))
	dface = ncout.createDimension('nMaxSCHISM_hgrid_face_nodes', 4)
	dtime =ncout.createDimension('time', 1)
	dsigma =ncout.createDimension('nSCHISM_vgrid_layers', sigma)	
	done =ncout.createDimension('one', 1)
	dtwo =ncout.createDimension('two', 2)

	vtime = ncout.createVariable('time',np.float64,dimensions=('time'))
	vele = ncout.createVariable('SCHISM_hgrid_face_nodes',np.int32,dimensions=('nSCHISM_hgrid_face','nMaxSCHISM_hgrid_face_nodes'))
	vx = ncout.createVariable('SCHISM_hgrid_node_x',np.float64,dimensions=('nSCHISM_hgrid_node'))
	vy = ncout.createVariable('SCHISM_hgrid_node_y',np.float64,dimensions=('nSCHISM_hgrid_node'))
	vdepth = ncout.createVariable('depth',np.float64,dimensions=('nSCHISM_hgrid_node'))
	vsigma = ncout.createVariable('zcor',np.float64,dimensions=('time','nSCHISM_hgrid_node','nSCHISM_vgrid_layers'))



	vtime[:]=1
	vele[:]=ncin.variables['SCHISM_hgrid_face_nodes'][:]
	vx[:]=ncin.variables['SCHISM_hgrid_node_x'][:]
	vy[:]=ncin.variables['SCHISM_hgrid_node_y'][:]
	vdepth[:]=ncin.variables['depth'][:]

	for j,param in enumerate(params):
		i23d=ncin.variables[param].i23d
		ivs=ncin.variables[param].ivs
		for S in threshs:
			print param+'_'+S
			if ivs==1 and i23d==1: # elev
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node'))
			if ivs==2 and i23d==1: # dahv
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node','two'))
			if ivs==1 and i23d==2: # zcor
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node','nSCHISM_vgrid_layers'))
			if ivs==2 and i23d==2: # hvel
				ncout.createVariable(param+'_'+S,np.float64,dimensions=('time','nSCHISM_hgrid_node','nSCHISM_vgrid_layers','two'))					

	
	ncin.close()
	return ncout,nnodes
	
	

def process(fileout,dirin,prefix,Istart,Iend,params,thresh):


	all_files_tmp = [y for x in os.walk(dirin) for y in glob(os.path.join(x[0], prefix+'_*.nc'))]
	all_files=[]
	for file in all_files_tmp:
		[tmp,filenum]=file.replace('.nc','').split('_')
		if int(filenum)>=Istart and int(filenum)<=Iend:
			all_files.append(file)


	print "%.f files found" % len(all_files)
	ncin,nnodes=create_output_file(fileout,all_files[0],thresh,params)


	ncout=process_exce(ncin,all_files,params,threhs)

	
	ncout.close()

            
    
    
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(prog='extract_exceedance.py', usage='%(prog)s fileout dirout params thresholds')
	## main arguments
	
	parser.add_argument('fileout', type=str,help='name of the output file (without the extension)')
	parser.add_argument('dirout', type=str,help='name of the output where are the SCHISM files')
	parser.add_argument('parameter', type=str,nargs='+',help='name of the parameter to plot')
	parser.add_argument('thresholds', type=float,nargs='+',help='tresholds')	
	parser.add_argument('-prefix', type=str,help='prefix default:schout_',default='schout')
	parser.add_argument('-start', type=int,help='First file to take',default=1)
	parser.add_argument('-end', type=int,help='Last file to take',default=10000)

	args = parser.parse_args()
	

	### PRINT IT ALL
	print 'output name : %s' % (args.fileout)
	print 'Direcory : %s' % (args.dirout)
	print 'From file #%i and #%i' % (args.start,args.end)
	print 'Do parameters : %s' % (args.params)
	
	
	process(args.fileout,args.dirout,args.prefix,args.start,args.end,args.params,args.thresholds)
