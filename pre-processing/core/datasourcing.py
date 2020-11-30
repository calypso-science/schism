import os
import datetime
import get_opendap
import glob
import copy


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

    def clean_mercator(filetmp):
        import pdb;pdb.set_trace()
    def download_mercator(self,fileout,source,t0,t1):
        url=source.get('url')
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
        if 'z' in source:
            add_url=' --depth-min '+str(source.get('Grid').get('z'))+' --depth-max '+str(source.get('Grid').get('z2'))
        else:
            add_url=''

        TRY=10

        for var in nvar:
            url='python -m motuclient --motu '+url+' '+\
            '--service-id '+service+\
            ' --product-id '+product+\
            ' --longitude-min '+str(xmin)+' --longitude-max '+str(xmax)+' '+\
            '--latitude-min '+str(ymin)+' --latitude-max '+str(ymax)+' '+\
            '--date-min "'+t0.strftime('%Y-%m-%d %H:%M:00')+\
            '" --date-max "'+t1.strftime('%Y-%m-%d %H:%M:00')+\
            '" --variable '+var

            url+=add_url
            url+=' --out-dir '+root+\
            ' --out-name MERC_'+var+'_'+t0.strftime('%Y%m%d%H%M00')+'.nc'+\
            ' --user '+user+' --pwd '+pwd

            for itry in range(0,TRY):
                if logging:
                    logging.info('Try #%i for %s' % (itry,var))
                os.system(url)
                if os.path.isfile(os.path.join(root,prf+'_'+t0.strftime('%Y%m%d%H%M00')+'.nc')):
                    break

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
        days = daterange(self.t0, self.t1+ datetime.timedelta(days=1))
        date_str=[]
        for day in days[:-1]:
            self.logger.info( " Sourcing data for %s-%s-%s" %(day.year, day.month, day.day))
            if day == days[-2]: 
                # if it's the end of a month, extrapolate to future
                tend = day + datetime.timedelta(seconds=24*3600)
            else: # it not, we do  dt not want to extrpolate to future
                delta = source.get('dt')
                # to avoid UDS getting time record in the future,
                #     we need to hardwire for a valid last time record
                tend = day + datetime.timedelta(seconds=(24-delta)*3600)


            filetmp = "%s%s.000000.nc" %(source.get('id'), day.strftime("%Y%m%d"))
            filetmp=os.path.join(rootdir, "in", filetmp)

            if source['id'].lower()=='hycom':
                self.download_hycom(filetmp,source,day,tend)
                self.clean_hycom(filetmp)
            if source['id'].lower()=='mercator':
                self.download_mercator(filetmp,source,day,tend)
                self.clean_mercator(filetmp)
            elif source['id'].lower()=='uds':
                self.download_uds(filetmp,source,day,tend)
                self.clean_uds(filetmp)
                if source.get('type','')=='tide':
                    break
            else:
                self.logger.info('  Source not understood')



        self.concat_files(os.path.join(rootdir, "in"),filename,source.get('id'),source.get('type'))

