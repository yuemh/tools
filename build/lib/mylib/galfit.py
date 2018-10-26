import os
from astropy.io import fits
from subprocess import *
import pkg_resources

#xcard = '_XCENT'
#ycard = '_YCENT'
# galfit3
xcard = '_XC'
ycard = '_YC'

dir_data = pkg_resources.resource_filename('mylib', 'data')
exe_dir = dir_data + '/executable'

class galfitpars(object):
    def __init__(self):
        self.objtype = 'none'
        self.parameters = []
        for i in range(10):
            self.parameters.append(dict(val=None,fit=False,comment=''))
            self.paramNum = None
            self.subtract = True
    def getobjecttype(self):
        return self.objtype
    def _getParamIndex(self,paramName):
        return self.paramNum[paramName] - 1
    def setfit(self,paramName,dofit):
        pi = self._getParamIndex(paramName)
        self.parameters[pi]['fit'] = dofit
    def setxy(self,x,y,dofit=True):
        self.parameters[0]['val'] = (x,y)
        self.parameters[0]['fit'] = dofit
    def getxy(self,x,y):
        return self.parameters[0]['val']
    def setsubtract(self,dosubtract):
        self.subtract = dosubtract
    def dosubtraction(self):
        return self.subtract
    def setparam(self,paramName,val,dofit=True,comment='',default=None):
        pi = self._getParamIndex(paramName)
        if default==None:
            self.parameters[pi]['val'] = val
            self.parameters[pi]['fit'] = dofit
        else:
            self.parameters[pi]['val'] = default
            self.parameters[pi]['fit'] = False
        self.parameters[pi]['comment'] = comment

    def getparam(self,paramName):
        pi = self._getParamIndex(paramName)
        return self.parameters[pi]['val']
    def writeparams(self,parfile):
        for i,par in enumerate(self.parameters):
            pi = i+1
            if pi==2: continue # because of position, there is no #2
            if (par['val'] == None):
                parfile.write(' %1d) 0.0   0.0              # \n' % pi)
            elif (isinstance(par['val'],tuple)):
                # for the position (x,y) entry
                dat = (pi,)+par['val']+(par['fit'],par['fit'])+(par['comment'],)
                parfile.write(' %1d) %.4f  %.4f  %1d  %1d   # %s\n' % dat)
            else:
                dat = (pi,par['val'],par['fit'],par['comment'])
                parfile.write(' %1d) %.4f   %1d              # %s\n' % dat)

class galfit_sky(galfitpars):
    def __init__(self,skyval,freelist=[True]):
        galfitpars.__init__(self)
        self.objtype = 'sky'
        self.paramNum = {'sky':1}
        self.setparam('sky',skyval,freelist[0],'sky value')

class galfit_psf(galfitpars):
    def __init__(self,x,y,mag,pa=0.0,freelist=[True,True]):
        galfitpars.__init__(self)
        self.objtype = 'psf'
        self.paramNum = {'position':1,'mag':3,}
        self.setxy(x,y,freelist[0])
        self.setparam('mag',mag,freelist[1],'total magnitude')

class galfit_gaussian(galfitpars):
    def __init__(self,x,y,mag,fwhm,axis_ratio=None,pa=None,freelist=[True,True,True,True,True]):
        galfitpars.__init__(self)
        self.objtype = 'gaussian'
        self.paramNum = {'position':1,'mag':3,'fwhm':4,'axis_ratio':8,'pa':9}
        self.setxy(x,y,freelist[0])
        self.setparam('mag',mag,freelist[1],'total magnitude')
        self.setparam('fwhm',fwhm,freelist[2],'Gaussian FWHM')
        self.setparam('axis_ratio',axis_ratio,freelist[3],'axis ratio',default=1.0)
        self.setparam('pa',pa,freelist[4],'position angle',default=0.0)

class galfit_moffat(galfitpars):
    def __init__(self,x,y,mag,fwhm,n=1.0,axis_ratio=None,pa=None,freelist=[True,True,True,True,True,True]):
        galfitpars.__init__(self)
        self.objtype = 'moffat'
        self.paramNum = {'position':1,'mag':3,'fwhm':4,'n':5,'axis_ratio':8,'pa':9}
        self.setxy(x,y,freelist[0])
        self.setparam('mag',mag,freelist[1],'total magnitude')
        self.setparam('fwhm',fwhm,freelist[2],'Gaussian FWHM')
        self.setparam('n',n,freelist[3],'concentration index')
        self.setparam('axis_ratio',axis_ratio,freelist[4],'axis ratio',default=1.0)
        self.setparam('pa',pa,freelist[5],'position angle',default=0.0)

class galfit_sersic(galfitpars):
    def __init__(self,x,y,mag,Re,n=4.0,axis_ratio=None,pa=None,freelist=[True,True,True,True,True,True]):
        galfitpars.__init__(self)
        self.paramNum = {'position':1,'mag':3,'Re':4,'n':5,'axis_ratio':9,'pa':10}
        self.objtype = 'sersic'
        self.setxy(x,y,freelist[0])
        self.setparam('mag',mag,freelist[1],'total magnitude')
        self.setparam('Re',Re,freelist[2],'Reff')
        self.setparam('n',n,freelist[3],'concentration index')
        self.setparam('axis_ratio',axis_ratio,freelist[4],'axis ratio')
        self.setparam('pa',pa,freelist[5],'position angle')

class galfit(object):
    def __init__(self,image,fitname,fitdir='./'):
        self.fitname = fitname
        self.image = image
        self.output = fitdir+fitname+'.fits'
        self.noise = 'none'
        self.psf = 'none'
        self.psfovrsamp = 1
        self.bpmask = 'none'
        self.constraints = 'none'
        self.zeropt = 20.0
        self.displaytype = 'regular'
        self.outonly = False
        self.interactive = False
        self.objects = []
        self.fitdir = fitdir
        self.parfile = fitdir+self.fitname+'.par'
        self.lastparfile = ''
        self.fitresults = None
        self.galfit_output = None

    def setfitname(self,fitname):
        self.fitname = fitname
        self.parfile = self.fitdir+self.fitname+'.par'
        self.output = self.fitdir+fitname+'.fits'
    def setimagename(self,image):
        self.image = image
    def addobject(self,objpars):
        self.objects.append(objpars)
    def addobjects(self,objparlist):
        self.objects.extend(objparlist)
#	def delobject(self,objnum):
#		pass
    def setpsf(self,psffile):
        self.psf = psffile

    def setzeropoint(self,zeropoint):
        self.zeropt = zeropoint

    def setimageregion(self,region):
        self.region = '%d  %d  %d  %d' % region
    def setboxsize(self,boxsize):
        if (isinstance(boxsize,int)):
            self.boxsize = '%d  %d' % (boxsize,boxsize)
        else:
            self.boxsize = '%d  %d' % boxsize
    def setplatescale(self,platescale):
        if (isinstance(platescale,float)):
            self.platescale = '%.4f  %.4f' % (platescale,platescale)
        else:
            self.platescale = '%.4f  %.4f' % platescale
    def setbadpixmask(self,badpix):
        self.bpmask = badpix
    def setconstraints(self,constraintf):
        self.constraints = constraintf
    def writeparams(self):
        f = open(self.parfile,'w')
        f.write(
"""#
A) %(image)s
B) %(fitname)s.fits
C) %(noise)s
D) %(psf)s
E) %(psfovrsamp)d
F) %(bpmask)s
G) %(constraints)s
H) %(region)s
I) %(boxsize)s
J) %(zeropt).3f
K) %(platescale)s
O) %(displaytype)s
P) %(outonly)d
S) %(interactive)d
""" % vars(self))
        objnum = 1
        for obj in self.objects:
            f.write('\n # Object number %d\n' % objnum)
            f.write(' 0) %s\n' % obj.getobjecttype())
            obj.writeparams(f)
            #f.write(' Z) ',obj.getoutputimagetype())
            f.write(' Z) %1d\n\n' % (not obj.dosubtraction()))
            objnum = objnum + 1

        f.close()

    def fit(self):
        self.writeparams()
        failure = os.system(exe_dir +'/galfit '+self.parfile)
#		p = Popen(['/home/minghao/software/executable/galfit',self.parfile],stdout=PIPE)
#		self.galfit_output = p.communicate()

    def readfit(self):
        self.fitresults = read_fit_results(self.output)
#		return fitresults
#	def refit(self):



def read_fit_results(fitfile):
    fitim = fits.open(fitfile)
    fith = fitim[2].header
    fitdat = {}
    ntype = {}
    #print fith.keys()
    fitdat['chisqr'] = fith['CHI2NU']
    #for i in range(1,len(self.objects)+1):
    for i in range(1,10):
        istr = '%1d' % i
        #fittype = fith['COMP_'+istr]
        try:
            fittype = fith['COMP_'+istr]
        except KeyError:
            break
        try:
            ntype[fittype] += 1
        except KeyError:
            ntype[fittype] = 1
        ns = '%1d' % ntype[fittype]
        if (fittype == 'sky'):
            sdat = read_card(fith[istr+'_SKY'])
#           xdat = fith[istr+'_DSDX'].split('+/-')
#           ydat = fith[istr+'_DSDY'].split('+/-')
            fitdat['sky'] = dict(component=i,\
                                sky=sdat[0],\
                                skyerr=sdat[1])
#			                     grad=(float(xdat[0]),float(ydat[0])),
#			                     grad_err=(float(xdat[1]),float(ydat[1])))
        elif (fittype == 'psf'):
            xdat = read_card(fith[istr+xcard])
            ydat = read_card(fith[istr+ycard])
            mdat = read_card(fith[istr+'_MAG'])
            fn = 'psf'+ns
            fitdat[fn] = dict(component=i,
                              position=(float(xdat[0]),float(ydat[0])),
                              position_err=(float(xdat[1]),float(ydat[1])),
                              mag=float(mdat[0]),
                              mag_err=float(mdat[1]))
        elif (fittype == 'gaussian'):
            xdat = read_card(fith[istr+xcard])
            ydat = read_card(fith[istr+ycard])
            mdat = read_card(fith[istr+'_MAG'])
            fwhmdat = read_card([istr+'_FWHM'])
            ardat = read_card(fith[istr+'_AR'])
            padat = read_card(fith[istr+'_PA'])
            fn = 'gaussian'+ns
            fitdat[fn] = dict(component=i,
                              position=(float(xdat[0]),float(ydat[0])),
                              position_err=(float(xdat[1]),float(ydat[1])),
                              mag=float(mdat[0]),
                              mag_err=float(mdat[1]),
                              fwhm=float(fwhmdat[0]),
                              fwhm_err=float(fwhmdat[1]),
                              axis_ratio=float(ardat[0]),
                              axis_ratio_err=float(ardat[1]),
                              pa=float(padat[0]),
                              pa_err=float(padat[1]))
        elif (fittype == 'sersic'):
            xdat = read_card(fith[istr+xcard])
            ydat = read_card(fith[istr+ycard])
            mdat = read_card(fith[istr+'_MAG'])
            Redat = read_card(fith[istr+'_RE'])
            ndat = read_card(fith[istr+'_N'])
            ardat = read_card(fith[istr+'_AR'])
            padat = read_card(fith[istr+'_PA'])
            fn = 'sersic'+ns
            fitdat[fn] = dict(component=i,
                              position=(float(xdat[0]),float(ydat[0])),
                              position_err=(float(xdat[1]),float(ydat[1])),
                              mag=float(mdat[0]),
                              mag_err=float(mdat[1]),
                              Re=float(Redat[0]),
                              Re_err=float(Redat[1]),
                              n=float(ndat[0]),
                              n_err=float(ndat[1]),
                              axis_ratio=float(ardat[0]),
                              axis_ratio_err=float(ardat[1]),
                              pa=float(padat[0]),
                              pa_err=float(padat[1]))
    return fitdat




def read_card(string):
    strlist=string.split('+/-')
    if len(strlist)==1:
        return (float(strlist[0][1:-1]),-1)
    elif len(strlist)==2:
        if strlist[0][0]=='*':
            return (float(strlist[0][1:-2]),float(strlist[1][2:-1]))
        else:
            return (float(strlist[0]),float(strlist[1]))
