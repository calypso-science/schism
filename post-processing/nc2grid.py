import os
from os.path import join 
import sys


def nc2grid(ncfile,):
    grbfile=ncfile[:-3]+'.grb'
    tmp_folder=os.path.split(grbfile)[0]
    tmp_folder=join(tmp_folder,'in')
    os.system('mkdir %s' % tmp_folder)
    for v in ['u','v']:
        os.system('ncks -v %s -d lev,0 %s %s' %(v,
            ncfile,
            join(tmp_folder,v+'.nc')))
        os.system('ncwa -O -a lev %s %s' %(
            join(tmp_folder,v+'.nc'),
            join(tmp_folder,v+'.nc')))
        os.system('cdo splitname %s %s' % (
            join(tmp_folder,v+'.nc'),
            join(tmp_folder,'split_')))
        os.system('ncrename -O -v %s,%s %s %s' % (v,v+'curr',
            join(tmp_folder,'split_'+v+'.nc'),
            join(tmp_folder,'split_'+v+'.nc'))) 
        os.system('cdo -f grb copy %s %s' % (
            join(tmp_folder,'split_'+v+'.nc'),
            join(tmp_folder,v+'.grb'))) 


    os.system('cdo -f grb chparam,1,049 %s %s' % (
        join(tmp_folder,'u.grb'),
        join(tmp_folder,'u1.grb'))) 

    os.system('cdo -f grb chparam,1,050 %s %s' % (
        join(tmp_folder,'v.grb'),
        join(tmp_folder,'v1.grb'))) 

    os.system('cdo merge %s %s %s' % (
        join(tmp_folder,'u1.grb'),
        join(tmp_folder,'v1.grb'),
        join(tmp_folder,'current0.grb')))

    os.system('cdo griddes %s > %s' % (
        ncfile,
        join(tmp_folder,'truegrid.txt')))

    os.system('cdo setgrid,%s %s %s' % (
        join(tmp_folder,'truegrid.txt'),
        join(tmp_folder,'current0.grb'),
        grbfile))

    os.system('rm -rf %s' % tmp_folder)