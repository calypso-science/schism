import os
import datetime
import get_opendap
import glob
import copy
import xarray as xr
import cdsapi

def daterange(tstart, tend, delta=86400,typ='seconds'):
    days = []
    d = tstart
    if typ=='seconds':
        delta = datetime.timedelta(seconds=delta)
    elif typ=='months':
        delta = relativedelta(months=delta)#seconds=delta)

    while d <= tend:
        days.append(d)
        d += delta
    return days

class download_data(object):

    def __init__(self,t0,t1,logger=None):
        '''Docstring'''  

        if logger:
            self.logger = logger

        self.t0=t0
        self.t1=t1

    def clean_pw(self,filein):
        os.system('mv %s %s' % (filein,filein+'.grb'))
        try:
            ds=xr.open_dataset(filein+'.grb', engine="cfgrib")
            ds.to_netcdf(filein)
        except:
            for v in ['10u','10v','msl']:
                ds=xr.open_dataset(filein+'.grb', engine="cfgrib",backend_kwargs={'filter_by_keys':{'shortName': '%s' % v}})
                ds.to_netcdf(os.path.join(os.path.split(filein)[0],'tmp_%s.nc' % v))
            os.system("mv %s %s" % (os.path.join(os.path.split(filein)[0],'tmp_msl.nc'),filein))
            os.system('ncks -A -v u10 %s %s' % (os.path.join(os.path.split(filein)[0],'tmp_%s.nc' % '10u'),filein))
            os.system('ncks -A -v v10 %s %s' % (os.path.join(os.path.split(filein)[0],'tmp_%s.nc' % '10v'),filein))
            os.system('rm %s' % os.path.join(os.path.split(filein)[0],'tmp_*.nc'))

        os.system("ncks -O -C -x -v step %s %s"%(filein, filein))
        os.system("ncks -O -C -x -v time %s %s"%(filein, filein))
        os.system("ncrename -O -v valid_time,time %s %s"%(filein, filein))
        os.system("ncks -O -3 %s %s"%(filein, filein))
        os.system("ncrename -O -d step,time %s %s"%(filein, filein))
        os.system("ncks -O -4 %s %s"%(filein, filein))
#        os.system("ncks -O --mk_rec_dmn time %s %s" %(filein, filein)) 
        os.system("ncpdq -O -U %s %s" %(filein, filein)) 
        os.system("rm %s" %(filein+'.grb*'))

        #os.system("ncatted -O -a _FillValue,,o,f,9.96920996838687e+36 %s %s" %(filein, filein))
    
    def download_olympics(self,fileout,source,t0,t1):

        url='wget -O '+fileout+' "'+source.get('url')

        url+='&username='+source.get('user')
        url+='&password='+source.get('pass')

        
        root,filename=os.path.split(fileout)

        for itry in range(0,10):
            self.logger.info('Try #%i for %s' % (itry,'olympics'))
            os.system(url)
            if os.path.isfile(os.path.join(root,filename)):
                break


    def download_pw(self,fileout,source,t0,t1):


        url='wget -O '+fileout+' "'+source.get('url')


        url+='username='+source.get('user')
        url+='&password='+source.get('pass')




        if (source.get('Grid')['y2']==source.get('Grid')['y']) & (source.get('Grid')['x2']==source.get('Grid')['x']):
            url+='&lat='+str(source.get('Grid')['y'])
            url+='&lon='+str(source.get('Grid')['x'])   
            url+='&timestep='+str(source.get('dt'))  
            url+='&source='+source.get('product')   
            url+='&compress=false&forecast=all&Z=10'
            url+='&resolution='+str(source.get('Grid')['dx']) 
            var='all'  

        else:
            url+='&nlat='+str(source.get('Grid')['y2'])
            url+='&slat='+str(source.get('Grid')['y'])
            url+='&elon='+str(source.get('Grid')['x2'])
            url+='&wlon='+str(source.get('Grid')['x'])
            url+='&time='+str(source.get('dt'))
            url+='&model='+source.get('product')
            url+='&res='+str(source.get('Grid')['dx'])

            #url+='&time=10'
            nvar=copy.deepcopy(source.get('vars'))
            url+='&variables='
            varname=[]
            if 'u10' in nvar:
                varname.append('wind')
            if 'msl' in nvar:
                varname.append('pressure')

            for var in varname:
                url+=var+','

            url=url[:-1]
        
        

        url+='"'
        print(url)
        root,filename=os.path.split(fileout)

        for itry in range(0,10):
            self.logger.info('Try #%i for %s' % (itry,var))
            os.system(url)
            if os.path.isfile(os.path.join(root,filename)):
                break

    def clean_mercator(self,filein):
        os.system("ncrename -O -d .depth,lev %s %s" %(filein, filein))
        os.system("ncrename -O -v depth,lev %s %s" %(filein, filein))
        os.system("ncks -O --mk_rec_dmn time %s %s" %(filein, filein)) 
        os.system("ncpdq -O -U %s %s" %(filein, filein)) 
        os.system("ncatted -O -a _FillValue,,o,f,9.96920996838687e+36 %s %s" %(filein, filein))
    def download_mercator(self,fileout,source,t0,t1):
        
        service=source.get('service')
        product=source.get('product')
        xmin=source.get('Grid')['x']
        xmax=source.get('Grid')['x2']
        ymin=source.get('Grid')['y']
        ymax=source.get('Grid')['y2']
        nvar=source.get('vars')
        pwd=source.get('pass')
        user=source.get('user')
        root,filename=os.path.split(fileout)
        add_url=' --depth-min '+str(source.get('Grid').get('z',0))+' --depth-max '+str(source.get('Grid').get('z2',6000))


        url='python3 -m motuclient --motu '+source.get('url')+' '+\
        '--service-id '+service+\
        ' --product-id '+product+\
        ' --longitude-min '+str(xmin)+' --longitude-max '+str(xmax)+' '+\
        '--latitude-min '+str(ymin)+' --latitude-max '+str(ymax)+' '+\
        '--date-min "'+(t0-datetime.timedelta(hours=12)).strftime('%Y-%m-%d %H:%M:00')+\
        '" --date-max "'+(t1+datetime.timedelta(hours=12)).strftime('%Y-%m-%d %H:%M:00')+'"'

        var_url=''
        for var in nvar:
            var_url+=' --variable '+var

        url+=var_url
        url+=add_url
        url+=' --out-dir '+root+\
        ' --out-name '+filename+\
        ' --user '+user+' --pwd '+pwd
        print(url)
        for itry in range(0,10):
            self.logger.info('Try #%i for %s' % (itry,var))
            os.system(url)
            if os.path.isfile(os.path.join(root,filename)):
                break
    def download_ecmwf(self,fileout,source,t0,t1):
        c = cdsapi.Client()

        service=source.get('service')
        product=source.get('product')
        xmin=source.get('Grid')['x']
        xmax=source.get('Grid')['x2']
        ymin=source.get('Grid')['y']
        ymax=source.get('Grid')['y2']
        nvar=source.get('vars')

        dwnl_opt={}
        dwnl_opt['product_type']=product
        dwnl_opt['variable']=nvar
        dwnl_opt['format']= 'netcdf'
        dwnl_opt['year']=t0.strftime('%Y')
        dwnl_opt['month']=t0.strftime('%m')
        dwnl_opt['day']=t0.strftime('%d')
        dwnl_opt['time']=[
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ]
        dwnl_opt['area']=[source.get('Grid')['y2'],source.get('Grid')['x'],source.get('Grid')['y'],source.get('Grid')['x2']]
        
        c.retrieve(service,dwnl_opt,fileout)
    def clean_ecmwf(self,filein):

        os.system("ncks -O --mk_rec_dmn time %s %s" %(filein, filein)) 
        os.system("ncpdq -O -U %s %s" %(filein, filein)) 
        os.system("ncatted -O -a _FillValue,,o,f,9.96920996838687e+36 %s %s" %(filein, filein))

    def download_hycom(self,fileout,source,t0,t1):

        # try:

            get_opendap.main(fileout,t0.strftime('%Y-%m-%d %H:%M:%S'),\
                t1.strftime('%Y-%m-%d %H:%M:%S'),\
                source.get('Grid')['y'],source.get('Grid')['y2'],\
                source.get('Grid')['x'],source.get('Grid')['x2'],\
                var=source.get('vars'),\
                z=source.get('Grid').get('z',0),Z=source.get('Grid').get('z2',0),\
                Source=source.get('dset'))
                
        # except:
        #     print "HYCOM request FAILED for %s" %(t0)
    def clean_hycom(self,filein):
        os.system("ncrename -O -d .depth,lev %s %s" %(filein, filein))

        os.system("ncks -O --mk_rec_dmn time %s %s" %(filein, filein)) 
        os.system("ncpdq -O -U %s %s" %(filein, filein)) 
        os.system("ncatted -O -a _FillValue,,o,f,9.96920996838687e+36 %s %s" %(filein, filein))
    
    def clean_uds(self,filein):
        os.system("ncks -O --mk_rec_dmn time %s %s" %(filein, filein))


    def concat_files(self,indir,mergefile,fileid,filetype):
        filelist = glob.glob( indir + "/%s*.nc"% fileid ) 

        file_tmp=os.path.join(indir,'%s*' % fileid) 
        self.logger.info( " Merging daily files into %s" %(mergefile))

        if filetype=='tide':
            os.system("mv %s %s" % (filelist[0],mergefile))
        elif fileid=='predictwind':
            os.system("mv %s %s" % (filelist[0],mergefile))
        else:
            os.system("ncrcat --fl_fmt=classic %s %s" %(file_tmp, mergefile))
            os.system("ncks -O --no_rec_dmn time %s %s" %(mergefile, mergefile))
            for f in filelist: os.remove(f)


    def download_uds(self,fileout,source,day,tend):


        from pymo_stuff import Data,DataList,Grid, Times

        src=copy.deepcopy(source)

        grid = Grid(**src.pop('Grid'))
        idd=src.pop('id')
        vaars = src.pop('vars')
        url = src.pop('url')
        filout=src.pop('filename')
        try:
            var=src.pop('var')
        except:
            None
        data = Data( idd, grid, vaars, url,**src) 
        dtimes = Times(t0=day, t1=tend)

        data.get(dtimes, fileout, self.logger, [])

        #if not data.get(dtimes, fileout, self.logger, []):
        #    self.logger.info( " UDS request FAILED for %s" %(day))

        

    def get_input_data(self,source,tsleep=30):
        """
        Request input data 
        """
        filename=source['filename']
        self.logger.info('  Sourcing external data files')
        if os.path.isfile(filename):
            self.logger.info('  file %s already exists' % filename)
            return


        rootdir=os.path.dirname(os.path.abspath(filename))
        if not os.path.exists(os.path.join(rootdir, "in")):
            os.makedirs(os.path.join(rootdir, "in"))


        # Download as daily file
        delta = source.get('dt')
        days = daterange(self.t0, self.t1+ datetime.timedelta(seconds=delta*3600))

        if self.t1>datetime.datetime.now() and source['id'].lower()=='predictwind':
            self.logger.info('  FORECAST')
            days = daterange(self.t0, self.t0+ datetime.timedelta(days=1))

        date_str=[]

        for day in days[:-1]:
            self.logger.info( " Sourcing data for %s-%s-%s" %(day.year, day.month, day.day))
            tend = day + datetime.timedelta(seconds=(24-delta)*3600)

            try:

                filetmp = "%s%s.000000.nc" %(source.get('id'), day.strftime("%Y%m%d"))
                filetmp=os.path.join(rootdir, "in", filetmp)

                if source['id'].lower()=='hycom':
                    self.download_hycom(filetmp,source,day,tend)
                    self.clean_hycom(filetmp)
                elif source['id'].lower()=='mercator':
                    self.download_mercator(filetmp,source,day,tend)
                    self.clean_mercator(filetmp)
                elif source['id'].lower()=='ecmwf':
                    self.download_ecmwf(filetmp,source,day,tend)
                    self.clean_ecmwf(filetmp)
                elif source['id'].lower()=='predictwind':
                    self.download_pw(filetmp,source,day,tend)
                    self.clean_pw(filetmp) 
                elif source['id'].lower()=='olympics':
                    self.download_olympics(filetmp,source,day,tend)
                    self.clean_pw(filetmp)               
                elif source['id'].lower()=='uds':
                    self.download_uds(filetmp,source,day,tend)
                    self.clean_uds(filetmp)
                    if source.get('type','')=='tide':
                        break

                else:
                    self.logger.info('  Source not understood')

            except:
                self.logger.info( " !!!! Not Found: %s-%s-%s !!!!" %(day.year, day.month, day.day))


        self.concat_files(os.path.join(rootdir, "in"),filename,source.get('id'),source.get('type'))

