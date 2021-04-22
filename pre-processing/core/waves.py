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
    
    def __init__(self,root, hgrid,wwm,t0,t1,deltc,h0,
        logger=None):
        '''Constructor'''

        if logger:
            self.logger = logger
            
        self.hgrid=hgrid
        self.config={}
        self.userconfig={}
        self.userconfig['wwm'] = wwm
        self.path=root

        self.config['wwm']=self._get_default_params()
        self._replace_default_params('wwm')

        self._set_params_relations(t0.strftime('%Y%m%d.%H%M%S'),t1.strftime('%Y%m%d.%H%M%S'))

        self.config['wwm']['DELTC']= deltc#MUST match dt*nstep_wwm in SCHISM!
        self.config['wwm']['DMIN']= h0


    def _make_bnd(self,fileout):
        values=np.zeros((len(self.hgrid.hgrid.x)))
        for bnd in self.hgrid.mesh.boundaries[None]:
            tmp=self.hgrid.mesh.boundaries[None][bnd]['indexes']
            tmp=[int(y)-1 for y in tmp]
            values[tmp]=2


        hgridWWM=copy.deepcopy(self.hgrid.mesh)
        hgridWWM.values[:]=values*-1
        hgridWWM.write(fileout)

    def _make_grid(self):
        hgridWWM=copy.deepcopy(self.hgrid.mesh)
        hgridWWM.write(os.path.join(self.path,'hgrid_WWM.gr3'))
        
    def _make_bndfile(self,ww3,filename):
        with open(filename,'w') as file: 
            ww3_order=['dir','fp','hs','spr','t02']
            for key in ww3_order:
                if key not in ww3:
                    self.logger.info("\t!!!!Need Key:%s!!!!"%key)
                    sys.exit(-1)
                file.write(ww3[key]+'\n')

    def make_wave(self,cmdfile,paramfile):
        '''Docstring'''        
        self.make_config( cmdfile,paramfile)

    def make_forcing(self,filewave,ww3=None):
        if not os.path.isfile(os.path.join(self.path,'hgrid_WWM.gr3')):
            self.logger.info("\tWriting:%s" %'hgrid_WWM.gr3')
            self._make_grid()

        if not os.path.isfile(os.path.join(self.path,'wwmbnd.gr3')):
            self.logger.info("\tWriting:%s" %'wwmbnd.gr3')
            self._make_bnd(os.path.join(self.path,'wwmbnd.gr3'))

        if not os.path.isfile(os.path.join(self.path,filewave)) and ww3:
            self.logger.info("\tWriting:%s" %filewave)
            self._make_bndfile(ww3,os.path.join(self.path,filewave))

    def _get_default_params(self):
        ''' Docstring'''
        default={}
        default['LQSTEA']=False
        default['LSPHE']=True
        default['MDC']=36
        default['MSC']=36
        default['LSLOP']=False
        default['SLMAX']=0.2
        default['LBCSE']=True
        default['LBINTER']=True
        default['LBCWA']=True
        default['LBCSP']=False
        default['LINHOM']=True
        default['LBSP1D']=False
        default['LBSP2D']=False
        default['DELTC2']=1
        default['IBOUNDFORMAT']=3
        default['LINDSPRDEG']=False
        default['LPARMDIR']=False
        default['MESNL']=1
        default['MESIN']=1
        default['IFRIC']=1
        default['MESBF']=1
        default['FRICC']=.067
        default['MESBR']=1
        default['IBREAK']=1
        default['ICRIT']=1
        default['BRCR']=0.78
        default['a_BRCR']=0.76
        default['b_BRCR']=0.29
        default['min_BRCR']=0.25
        default['max_BRCR']=0.8
        default['a_BIPH']=0.2
        default['BR_COEF_METHOD']=1
        default['B_ALP']=0.5
        default['ZPROF_BREAK']=2
        default['BC_BREAK']=1
        default['IROLLER']=0
        default['ALPROL']=0.85
        default['MEVEG']=0
        default['LMAXETOT']=True
        default['MESDS']=1
        default['MESTR']=1
        default['TRICO']=0.1
        default['TRIRA']=5
        default['TRIURS']=0.1
        default['ICOMP']=3
        default['AMETHOD']=7
        default['SMETHOD']=1
        default['ROLMETHOD']=2
        default['DMETHOD']=2
        default['RTHETA']=0.5
        default['LITERSPLIT']=False
        default['LFILTERTH']=False
        default['MAXCFLTH']=1.0
        default['FMETHOD']=1
        default['LFILTERSIG']=False
        default['MAXCFLSIG']=1.0
        default['MELIM']=1
        default['LIMFAK']=0.1
        default['LDIFR']=False
        default['IDIFFR']=1
        default['LCONV']=False
        default['LCFL']=False
        default['NQSITER']=1
        default['QSCONV1']=0.98
        default['QSCONV2']=0.98
        default['QSCONV3']=0.98
        default['QSCONV4']=0.98
        default['QSCONV5']=0.98
        default['LEXPIMP']=False
        default['FREQEXP']=0.1
        default['EPSH1']=0.01
        default['EPSH2']=0.01
        default['EPSH3']=0.01
        default['EPSH4']=0.01
        default['EPSH5']=0.01
        default['LVECTOR']=False
        default['IVECTOR']=2
        default['LADVTEST']=False
        default['LCHKCONV']=False
        default['DTMIN_DYN']=1
        default['NDYNITER']=100
        default['DTMIN_SIN']=1
        default['DTMIN_SNL4']=1
        default['DTMIN_SDS']=1
        default['DTMIN_SNL3']=1
        default['DTMIN_SBR']=.1
        default['DTMIN_SBF']=1
        default['NDYNITER_SIN']=10
        default['NDYNITER_SNL4']=10
        default['NDYNITER_SDS']=10
        default['NDYNITER_SBR']=10
        default['NDYNITER_SNL3']=10
        default['NDYNITER_SBF']=10
        default['LSOUBOUND']=False
        default['WAE_SOLVERTHR']=1e-6
        default['MAXITER']=1000
        default['PMIN']=1
        default['LNANINFCHK']=False
        default['LZETA_SETUP']=False
        default['ZETA_METH']=0
        default['BLOCK_GAUSS_SEIDEL']=True
        default['LNONL']=False
        default['ASPAR_LOCAL_LEVEL']=0
        default['L_SOLVER_NORM']=False
        default['LACCEL']=False
        default['DELTC3']=86400
        default['OUTSTYLE']="'NO'"
        default['FILEWAVE']="'bndfiles.dat'"
        return default
    def _replace_default_params(self,module):
        for key in self.config[module]:
            if key in self.userconfig[module]:
                self.config[module][key]=self.userconfig[module][key]   
    
    def _set_params_relations(self,t0,t1):
        self._update_timings(t0,t1)
        self._choose_mode()
    def _choose_mode(self):
        mode=self.userconfig['wwm'].get('mode','ww3')
        if mode=='ww3':
            pass
        elif mode=='spectra':
            self.config['wwm']['LBCWA']=False
            self.config['wwm']['LBCSP']=True
            self.config['wwm']['LBSP2D']=True
            self.config['wwm']['IBOUNDFORMAT']=6
            self.config['wwm']['FILEWAVE']="'"+self.config['wwm']['FILEWAVE']+"'"
        else:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('!!!!!!!!Mode must be WW3 or Spectra !!!!!!!!!!!!!')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            import sys;sys.exit(-1)


    def _update_timings(self,t0,t1):
        self.config['wwm']['BEGTC']="'"+t0+"'"#'"'+"'"+t0+"'"+'"'
        self.config['wwm']['ENDTC']="'"+t1+"'"#'"'+"'"+t1+"'"+'"'
        


    def make_config(self, cmdfile,paramfile):
        '''Docstring''' 

        self.logger.info("\tWriting:%s" %paramfile)
        cmdblank = Template(filename=cmdfile,input_encoding='utf-8',output_encoding='utf-8')
        buf = StringIO()
        stri='Context(buf,'
        
        for key in self.config['wwm']:
            if self.config['wwm'][key] is True:
                ss="'T'"
            elif self.config['wwm'][key] is False:
                ss="'F'"
            else:
                ss=self.config['wwm'][key]
                if type(ss) != int and type(ss) != float:
                    ss='"'+ss+'"'



            stri=stri+key+'='+str(ss)+','

        stri=stri[0:-1]+')'
        ctx = eval(stri.encode('utf-8'))
        
        try:
            cmdblank.render_context(ctx)
        except:
            self.logger.info(exceptions.text_error_template().render())

    
        a=open(paramfile,'w')
        a.write(buf.getvalue())
        a.close()





