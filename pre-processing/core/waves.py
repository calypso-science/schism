#!/usr/bin/env python3.7
from mako.template import Template
from mako.runtime import Context
from io import StringIO
from mako import exceptions
import ast
from matplotlib.dates import date2num
import numpy as np
import os,copy

class Wave(object):
    
    def __init__(self,root, hgrid,hydro,t0,t1,
        logger=None):
        '''Constructor'''

        if logger:
            self.logger = logger
            
        self.hgrid=hgrid
        self.config={}
        self.userconfig={}
        self.userconfig['hydro'] = hydro
        self.path=root

    def _make_bnd(self):
        values=np.zeros((len(self.hgrid.hgrid.x)))
        for bnd in self.hgrid.mesh.boundaries[None]:
            tmp=self.hgrid.mesh.boundaries[None][bnd]['indexes']
            tmp=[int(y)-1 for y in tmp]
            values[tmp]=2


        hgridWWM=copy.deepcopy(self.hgrid.mesh)
        hgridWWM.values[:]=values
        hgridWWM.write(os.path.join(self.path,'wwmbnd.gr3'))

    def _make_grid(self):
        hgridWWM=copy.deepcopy(self.hgrid.mesh)
        hgridWWM.write(os.path.join(self.path,'hgrid_WWM.gr3'))
        
    def make_wave(self):
        '''Docstring''' 
        if not os.path.isfile(os.path.join(self.path,'hgrid_WWM.gr3')):
            self.logger.info("\tWriting:%s" %'hgrid_WWM.gr3')
            self._make_grid()

        if not os.path.isfile(os.path.join(self.path,'wwmbnd.gr3')):
            self.logger.info("\tWriting:%s" %'wwmbnd.gr3')
            self._make_bnd()

        # cmdblank = Template(filename=cmdfile,input_encoding='utf-8',output_encoding='utf-8')
        # buf = StringIO()
        # stri='Context(buf,'
        
        # for key in self.config[module]:
        #     stri=stri+key+'='+str(self.config[module][key])+','

        # stri=stri[0:-1]+')'
        # ctx = eval(stri.encode('utf-8'))
        

        # #cmdblank.render_context(ctx)
        # try:
        #     cmdblank.render_context(ctx)
        # except:
        #     self.logger.info(exceptions.text_error_template().render())

	
        # a=open(paramfile,'w')
        # a.write(buf.getvalue())
        # a.close()



