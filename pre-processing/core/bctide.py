from matplotlib.dates import date2num
import numpy as np
from ttide.t_getconsts import t_getconsts
from ttide.t_vuf import t_vuf
import copy
from tidal_tools import extract_HC
from ttide.t_getconsts import t_getconsts

TRACERS_MODULE=['GEN','AGE','SED3D','EcoSim','ICM','CoSiNE','FIB','TIMOR']

def t_equilib(freq,doodsonamp,doodsonspecies,lat0):
  g=9.81;            # m/s^2;
  erad=6365;         # km
  earthmoond=3.84e5; # km
  Mmoon=7.38e22;     # kg
  Mearth=5.977e24;   # kg
  Gravconst=6.658e-11;  # m^3/kg/s^2

  # There appears to be a typo in Godin's text, and this
  # should likely be checked against Doodson's original.
  # This is what I *think* it should be.
  G=3/4*Mmoon*(erad*1e3)**3/(earthmoond*1e3)**2/Mearth;

  # The 4/3 is to correct for the 3/4 in G
  gfac=Gravconst*Mearth/(erad*1e3)**2*(4/3);
    
  slat=np.sin(lat0*np.pi/180);
  clat=np.cos(lat0*np.pi/180);

  G1=np.zeros((6,1));

  # Latitude dependence of amplitude for various species -
  # + for A, -for B (from Godin, 1972).

  G1[3+0-1,0]=    0.5*G*(1-3*slat**2);
  G1[3-1-1,0]=      2*G*slat*clat;
  G1[3+1-1,0]= .72618*G*clat*(1-5*slat**2);
  G1[3-2-1,0]=2.59808*G*slat*clat**2;
  G1[3+2-1,0]=        G*clat**2;
  G1[3+3-1,0]=        G*clat**3;

  idx=[int(idx) for idx in doodsonspecies+3-1]
  #amp=np.abs(doodsonamp/gfac*G1[np.ix_(idx)].T)[0]
  amp=doodsonamp/3.747394476289734
  return amp
   

def Calculate(lat, t0, cons,typ=0):

    const = t_getconsts(np.array([]))
    Const= [con.decode('UTF-8') for con in const[0]['name']] 


    consindex = [Const.index(con.ljust(4)) for con in cons]

    # V: astronomical phase, U: nodal phase modulation, F: nodal amplitude correction
    v,u,f = t_vuf('nodal', np.array(t0), consindex, lat)
    tear = 360.*(v+u)
    tfreq = (2*np.pi)*const[0]['freq'][consindex]/3600.
    talpha = const[0]['name'][consindex]

    if typ==1:

      tpspec=np.abs(const[0]['doodsonspecies'][consindex])
      tpamp=t_equilib(const[0]['freq'][consindex],const[0]['doodsonamp'][consindex],const[0]['doodsonspecies'][consindex],lat)
      return talpha,tpspec,tpamp, tfreq, tear, f

    else:
      return talpha, tfreq, tear, f

def develop_bnd(obc):
    # check if we have a loop in the dictionnary
    if 'cons' in obc:
        obc.pop('cons')
    if 'tip_dp' in obc:
        obc.pop('tip_dp')
    if 'tp cons' in obc:
        obc.pop('tp cons')

    btype=dict()
    
    for nk in obc.keys():
        if type(nk) is not int:
            [aa,bb]=nk.split('-')
            for nkk in range(int(aa)-1,int(bb)):
                btype[int(nkk)]=obc[nk]
                
        else:
    
            btype[int(nk)-1]=obc[nk]

    return btype

class BCinputs(object):
    """
    Class that manages input file generation for a model simulation.
    """

    def __init__(self,obc,hgrid,lat0,t0, logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.lat0 = lat0
        self.t0= t0
        self.t0dec=date2num(t0)
        self.obc=obc
        self.hgrid=hgrid
        self.nnode=hgrid.nnode




    def update_bctide(self,filin,rnday,soft='./tide_fac'):
        f=open('date_param.txt', 'w')
        # 48-character start time info string (only used for visualization with xmvis6)
        f.write('%s%i\n' % (self.t0.strftime('%H,%d,%m,%Y,'),rnday))
        f.flush()
        f.close()
        os.system(soft)
        


    
    def _write_iettype(self,btypes,iettype):
            

            ## TAKE CARE OF iettype
            if iettype==2: #this boundary is forced by a constant elevation
                ethconst =btypes['iettype']['const']
                self.bctides.write("%.2f\n" % (ethconst)) 

            elif iettype==3 or iettype==5: #forced in frequency domain
                for cons in self.iecons.keys():
                  self.bctides.write("%s\n" % (cons))
                  keys=self.iecons[cons].items()
                  eamp=[ k for k,v in keys if 'amp' in k][0]
                  epha=[ k for k,v in keys if 'pha' in k][0]

                  for n in range(0,len(self.iecons[cons][eamp])):
                    self.bctides.write("%.4f %.4f\n" % (self.iecons[cons][eamp][n],self.iecons[cons][epha][n]))


    def _write_ifltype(self,btypes,ifltype):

            ## TAKE CARE OF ifltype
        if ifltype==2: #forced in frequency domain
            const=btypes['ifltype']['const']
            self.bctides.write("%.4f\n" % (const))
        elif ifltype==3 or ifltype==5: #this boundary is forced by a constant elevatio
            for cons in self.ifcons.keys():
              self.bctides.write("%s\n" % (cons))
              keys=self.ifcons[cons].items()
              uamp=[ k for k,v in keys if 'amp' in k and 'u' in k][0]
              upha=[ k for k,v in keys if 'pha' in k and 'u' in k][0]
              vamp=[ k for k,v in keys if 'amp' in k and 'v' in k][0]
              vpha=[ k for k,v in keys if 'pha' in k and 'v' in k][0]
              for n in range(0,len(self.ifcons[cons][uamp])):
                self.bctides.write("%.4f %.4f %.4f %.4f\n" % (self.ifcons[cons][uamp][n],\
                                                              self.ifcons[cons][upha][n],\
                                                              self.ifcons[cons][vamp][n],\
                                                              self.ifcons[cons][vpha][n]))
        elif ifltype==-4: #time history of velocity 
            inflow = btypes['ifltype']['inflow']
            outflow = btypes['ifltype']['outflow']
            self.bctides.write("%.2f %.2f\n" % (inflow,outflow))
        elif ifltype==-1: #Flanther type radiation b.c.
            eta_m0 = btypes['ifltype']['eta_m0']
            qthcon = btypes['ifltype']['qthcon']
            self.bctides.write("%.2f %.2f\n" % (eta_m0,qthcon))


    def _write_tracers(self,btypes,module,flag):

        if flag==1 or flag==3 or flag==4: #
            tobc = btypes[module]['tobc']
            self.bctides.write("%.2f\n" % (tobc))
        elif flag==2: 
            const = btypes[module]['const']
            tobc = btypes[module]['tobc']
            self.bctides.write("%.2f\n" % (const))
            self.bctides.write("%.2f\n" % (tobc))

    def _write_bctides(self,filename,tpalpha=[],tpspec=[],tpamp=[],tpfreq=[],tpear=[],tpnf=[], talpha=[], tfreq=[], tear=[], tnf=[]):
        '''
           ----------------------------------- Elevation b.c. section ------
            iettype(j): (DEFAULT: 5)

                0: elevations are not specified for this boundary (in this case the velocity must be specified).

                1: time history of elevation on this boundary
                  > no input in bctides.in; time history of elevation is read in from elev.th (ASCII);

                2: forced by a constant elevation
                  > ethconst (constant elevation value for this segment)

                3: forced by tides
                  > for k=1, nbfr
                  >   alpha(k) !tidal constituent name
                  >   for i=1, nond(j) !loop over all open boundary nodes on this segment
                  >       emo((j,i,k) efa (j,i,k) !amplitude and phase for each node on this open boundary
                  >   end for i
                  > end for k

                4: space- and time-varying input
                  > no input in this file; time history of elevation is read in from elev2D.th(binary)

                5: combination of (3) and (4):
                   > time history of elevation is read in from elev2D.th, and then added to tidal B.C. specified below
                   > for k = 1, nbfr
                   >     alpha(k) !tidal constituent name
                   >     for i = 1, nond(j) !loop over all open boundary nodes on this segment
                   >         emo((j, i, k) efa(j, i, k) !amplitude and phase for each node on this open boundary
                   >     end for i
                   > end for k

            ----------------------------------- Velocity b.c. section -----------------------
            ifltype(j): (DEFAULT: 5)

                0: velocity not specified (no input needed)

                1: time history of discharge on this boundary
                  > no input in this file; time history of discharge is read in from flux.th (ASCII)

                2: this boundary is forced by a constant discharge
                  > vthconst (constant discharge; note that a negative number means inflow)

                3: vel. (not discharge!) is forced in frequency domain
                  > for k=1, nbfr
                  >    vmo(j,k) vfa(j,k) (uniform amplitude and phase along each boundary segment)
                  >    end for;

                4 or -4: 3D input
                        > time history of velocity (not discharge!) is read in from uv3D.th (binary)

                5: combination of (4) and tides
                  > time history of velocity (not discharge!) is read in from uv3D.th (binary) and then added to tidal velocity specified below
                  >    for k=1, nbfr
                  >       alpha(k) (tidal constituent name)
                  >       for i=1, nond(j) (loop over all open boundary nodes on this segment)
                  >           umo(j,i,k) ufa(j,i,k) vmo(j,i,k) vfa(j,i,k) !amplitude and phase for (u,v) at each node on this open boundary
                  >       end for i
                  >    end for k

                -1: Flanther type radiation b.c. (iettype must be 0 in this case)
                   > eta_mean (mean elevation below)
                   > for i=1,nond(j) (loop over all nodes)
                   >     eta_m0(i) !mean elev at each node
                   > end for i
                   > vn_mean (mean normal velocity)
                   > for i=1,nond(j)
                   >    qthcon(1:Nz,i,j) (mean normal velocity at the node; at all levels)
                   > end for i

            ----------------------------- Temperature b.c. section ---------------------------
            itetype: (DEFAULT: ?)

                0: temperature not specified (no input needed)

                1: time history of temperature on this boundary
                  > tobc !nudging factor (between 0 and 1 with 1 being strongest nudging) for inflow; time history of temperature will be read in from TEM_1.th (ASCII)

                2: this boundary is forced by a constant temperature
                  > tthconst (constant temperature on this segment)
                  > tobc (nudging factor (between 0 and 1) for inflow)

                3: initial temperature profile for inflow
                  > tobc (nudging factor (between 0 and 1) for inflow)

                4: 3D input
                  > tobc (nudging factor (between 0 and 1); time history of temperature is read in from TEM_3D.th (binary)

            ---------------------------- Salintiy B.C. section ------------------------------
            isatype: (DEFAULT: ?)

                the structure is similar to temperature

            ---------------------------- Tracers B.C. section -------------------------------
            If any tracer module is invoked, you also need the corresponding B.C. part for each tracer module

            isatype: (DEFAULT: 0)

                the structure is similar to temperature and salinity

            end for !j: open boundary segment
            *** Notes: the tidal amplitudes and phases can be generated using utility scripts shown on the web.'''


        self.bctides = open(filename, 'w')
        # 48-character start time info string (only used for visualization with xmvis6)
        self.bctides.write(self.t0.strftime('%d/%m/%Y %H:%M:%S\n'))
        # ntip tip_dp > # of constituents used in earth tidal potential; cut-off depth for applying tidal potential (i.e., it is not calculated when depth < tip_dp).
        if len(tpalpha)>0:
          self.bctides.write("%i %.f\n"%(len(tpalpha),self.obc.get('tip_dp',1000)))
          for k in range(0, len(tpalpha)):
              # tidal constituent name
            self.bctides.write("%s\n"%(str(tpalpha[k].decode('UTF-8')))) 
              # angular frequency (rad/s), nodal factor, earth equilibrium argument (deg)
            self.bctides.write("%i %.6f %.6f %.6f %.6f\n" % (tpspec[k],tpamp[k],tpfreq[k], tpnf[k], np.mod(tpear[k],360)))
        else:
          self.bctides.write("%i %.f\n"%(len(tpalpha),self.obc.get('tip_dp',1000)))


        # nbfr > # of tidal boundary forcing frequencies
        self.bctides.write("%d\n"%len(talpha))
        for k in range(0, len(talpha)):
            # tidal constituent name
        	self.bctides.write("%s\n"%(str(talpha[k].decode('UTF-8')))) 
            # angular frequency (rad/s), nodal factor, earth equilibrium argument (deg)
        	self.bctides.write("%.6f %.6f %.6f\n" % (tfreq[k], tnf[k], np.mod(tear[k],360)))

        


        #number of open boundary segments
        nopen=len(self.nnode)
        self.bctides.write('%s\n'%nopen)


        btypes=develop_bnd(self.obc)



        
        for k in range(0,nopen):
    # open boundary option
            iettype = btypes[k]['iettype']['value']
            ifltype = btypes[k]['ifltype']['value']
            itetype = btypes[k]['itetype']['value']
            isatype = btypes[k]['isatype']['value']

            
            self.bctides.write("%.f %.f %.f %.f %.f" % (len(self.nnode[k]),iettype,ifltype,itetype,isatype))
            ## add the tracers
            for modules in TRACERS_MODULE:
                if modules in btypes[k]:
                    self.bctides.write(" %.f" % (btypes[k][modules]['value']))
                else:
                    self.bctides.write(" %.f" % (0))
            
            self.bctides.write("\n" )


            self._write_iettype(btypes[k],iettype)
            self._write_ifltype(btypes[k],ifltype)


            self._write_tracers(btypes[k],'itetype',itetype) # temperature
            self._write_tracers(btypes[k],'isatype',isatype) # salinity

            for modules in TRACERS_MODULE:
                if modules in btypes[k]:
                    self._write_tracers(btypes[k],modules,btypes[k][modules]['value']) # salinity
                            
                        
        self.bctides.flush()
        self.bctides.close()


    def make_bctides(self,filename):
        '''Docstring'''
        
        if self.logger:
            self.logger.info("  Writing %s" %filename)


        talpha, tfreq, tear, tnf = Calculate(self.lat0, self.t0dec, self.obc['cons'].split(' '))
        if 'tp cons' in self.obc:
          tpalpha,tpspec,tpamp, tpfreq, tpear, tpnf = Calculate(self.lat0, self.t0dec, self.obc['tp cons'].split(' '),1)
        else:
          tpalpha=[]
          tpspec=[]
          tpamp=[]
          tpfreq=[]
          tpear=[]
          tpnf=[]


        

        btypes=develop_bnd(copy.deepcopy(self.obc))
        iet=[btypes[k]['iettype']['value'] for k in btypes.keys()]
        ifl=[btypes[k]['ifltype']['value'] for k in btypes.keys()]
        if (3 in iet) or (5 in iet):
          self.get_e_cons(btypes)
          

        if (3 in ifl) or (5 in ifl):
          self.get_uv_cons(btypes)

        self._write_bctides(filename,tpalpha,tpspec,tpamp,tpfreq,tpear,tpnf,talpha, tfreq, tear, tnf)

    def get_e_cons(self,btypes):
          Lat=[]
          Lon=[]  
          for b in range(0,len(btypes)):
            if (btypes[b]['iettype']['value']==3) or (btypes[b]['iettype']['value']==5):
              file=btypes[b]['iettype']['filename']
              var=btypes[b]['iettype']['vars']
           
              for node_i in self.nnode[b]:
                Lat.append(self.hgrid.latitude[int(node_i)])
                Lon.append(self.hgrid.longitude[int(node_i)])
          
              
          HC,tfreq,constidx=extract_HC(file,var,Lon,Lat, logger=self.logger)
          const = t_getconsts(np.array([]))
          self.iecons={}
          
          for i,n in enumerate(constidx):
            self.iecons[const[0]['name'][n].decode('UTF-8')]={}
            for k in HC.keys():
              k=k.split('_')[0]
              self.iecons[const[0]['name'][n].decode('UTF-8')][k+'_amp']=HC[k][i,0,:]
              self.iecons[const[0]['name'][n].decode('UTF-8')][k+'_pha']=np.mod(HC[k][i,2,:],360)


    def get_uv_cons(self,btypes):
          Lat=[]
          Lon=[]  
          for b in range(0,len(btypes)):
            if (btypes[b]['ifltype']['value']==3) or (btypes[b]['ifltype']['value']==5):
              file=btypes[b]['ifltype']['filename']
              var=btypes[b]['ifltype']['vars']
           
              for node_i in self.nnode[b]:
                Lat.append(self.hgrid.latitude[int(node_i)])
                Lon.append(self.hgrid.longitude[int(node_i)])
          
              
          HC,tfreq,constidx=extract_HC(file,var,Lon,Lat, logger=self.logger)
          const = t_getconsts(np.array([]))
          self.ifcons={}
          
          for i,n in enumerate(constidx):
            self.ifcons[const[0]['name'][n].decode('UTF-8')]={}
            for k in var:
              k=k.split('_')[0]
              self.ifcons[const[0]['name'][n].decode('UTF-8')][k+'_amp']=HC[k][i,0,:]
              self.ifcons[const[0]['name'][n].decode('UTF-8')][k+'_pha']=np.mod(HC[k][i,2,:],360)
