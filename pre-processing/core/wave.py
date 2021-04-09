#!/usr/bin/env python3.7
from mako.template import Template
from mako.runtime import Context
from io import StringIO
from mako import exceptions
import ast
from matplotlib.dates import date2num

class Wave(object):
    
    def __init__(self, grid,hydro,t0,t1,
        logger=None):
        '''Constructor'''

        if logger:
            self.logger = logger
            

        self.config={}
        self.userconfig={}
        self.userconfig['hydro'] = hydro
       

    def make_wave(self):
        '''Docstring''' 

        self.logger.info("\tWriting:%s" %'hgrid_WWM.gr3')
        import pdb;pdb.set_trace()
        self.logger.info("\tWriting:%s" %'wwmbnd.gr3')

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



