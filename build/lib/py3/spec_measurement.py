import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from math import *
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import pkg_resources

dir_data = pkg_resources.resource_filename('mylib', 'data')
filter_dir = dir_data + '/filter'

WFC_filter_dict=\
        {'F435W':filter_dir+'/HST/ACS/WFC/wfc_F435W.dat',\
         'F475W':filter_dir+'/HST/ACS/WFC/wfc_F475W.dat',\
         'F555W':filter_dir+'/HST/ACS/WFC/wfc_F555W.dat',\
         'F606W':filter_dir+'/HST/ACS/WFC/wfc_F606W.dat',\
         'F625W':filter_dir+'/HST/ACS/WFC/wfc_F625W.dat',\
         'F775W':filter_dir+'/HST/ACS/WFC/wfc_F775W.dat',\
         'F814W':filter_dir+'/HST/ACS/WFC/wfc_F814W.dat',\
         'F850LP':filter_dir+'/HST/ACS/WFC/wfc_F850LP.dat'}

SDSS_filter_dict=\
        {'SDSSu':filter_dir+'/SDSS/sdssu',\
         'SDSSg':filter_dir+'/SDSS/sdssg',\
         'SDSSr':filter_dir+'/SDSS/sdssr',\
         'SDSSi':filter_dir+'/SDSS/sdssi',\
         'SDSSz':filter_dir+'/SDSS/sdssz'}

IRAC_filter_dict=\
        {'CH1':filter_dir+'/IRAC/irac_ch1.res',\
         'CH2':filter_dir+'/IRAC/irac_ch2.res',\
         'CH3':filter_dir+'/IRAC/irac_ch3.res',\
         'CH4':filter_dir+'/IRAC/irac_ch4.res'}

PanStaars_filter_dict=\
        {'PSg':filter_dir + '/PanStarrs/PAN-STARRS_PS1.g.dat',\
        'PSi':filter_dir + '/PanStarrs/PAN-STARRS_PS1.i.dat'}

MIPS_filter_dict=\
        {'MIPS24':filter_dir+'/linhua/mips24.res'}

WFCAM_filter_dir=\
        {'H':filter_dir+'/WFCAM/WFCAM_H.res',\
         'K':filter_dir+'/WFCAM/WFCAM_K.res'}

WISE_filter_dict=\
        {'W1':filter_dir+'/WISE/w1.dat',
         'W2':filter_dir+'/WISE/w2.dat',
         'W3':filter_dir+'/WISE/w3.dat',
         'W4':filter_dir+'/WISE/w4.dat'}

MegaCam_filter_dict=\
        {'MegaCam_g': filter_dir + '/MegaCam/gS.dat',\
        'MegaCam_i': filter_dir + '/MegaCam/iS.dat'}

instru_filter_dict=\
        {'SDSS':SDSS_filter_dict,\
         'ACSWFC':WFC_filter_dict,\
         'IRAC':IRAC_filter_dict,\
         'MIPS':MIPS_filter_dict,\
         'WFCAM':WFCAM_filter_dir,\
         'PanStaars':PanStaars_filter_dict,\
         'WISE':WISE_filter_dict,\
         'MegaCam': MegaCam_filter_dict}


class lambdafunc(object):
    def update(self,wavelength,value,units):
        if not len(wavelength) == len(value):
            raise ValueError('The length of wavelength list\
                             does not equal to the length of value list')
        else:
            self.wavelength =\
                    np.array(wavelength)*float(1*units[0].cgs/u.Angstrom)
            self.value = np.array(value)
            self.units = [u.Angstrom.cgs,units[1].cgs]
            self.interval = \
                    (np.append(self.wavelength,0)\
                     - np.append(0,self.wavelength))[1:-1]
            self.midvalue = \
                    (np.append(self.value,0) + \
                     np.append(0,self.value))[1:-1] / 2.0

    def getvalue(self,x):
        return np.interp(x,self.wavelength,self.value)

    def plot(self):
        plt.plot(self.wavelength,self.value)
        plt.xlabel('Wavelength / $\AA$',fontsize=16)
        plt.ylabel('Arbitrary Unit',fontsize=16)
        plt.show()


class FiltCurve(lambdafunc):
    def __init__(self, wavelength, value, units=[u.Angstrom,u.Quantity(1)]):
        self.update(wavelength,value,units)


def read_filter(instrument, filtername, index=[0,1], uflag=0):
    filename = instru_filter_dict[instrument][filtername]
    data = np.loadtxt(filename)
    wavelength = data[:,index[0]]
    throughput = data[:,index[1]]

    if uflag == 1:
        filt = FiltCurve(wavelength,throughput,units=[u.nm,u.Quantity(1)])
    elif uflag == 2:
        filt = FiltCurve(wavelength,throughput,units=[u.um,u.Quantity(1)])
    elif uflag == 0:
        filt = FiltCurve(wavelength,throughput)

    return filt

class Spectrum(lambdafunc):
    def __init__(self, wavelength, value,\
                 units=[u.Angstrom,u.erg/u.s/u.cm/u.cm/u.Angstrom],\
                 mode='OBS'):
        self.update(wavelength, value, units)
        self.mode = mode

    def flux(self, filt):
        newflux = self.getvalue(\
                np.array(filt.wavelength*filt.units[0]/self.units[0]))
        filtered_flux = newflux*filt.value
        midflux = (np.append(filtered_flux,0) + \
                   np.append(0,filtered_flux))[1:-1]/2.0
        interval=np.abs(filt.interval)

        return np.sum(interval*midflux)*filt.units[0]*self.units[1]

    def magnitude(self, filt, style='AB'):
        '''
        Calculate the magnitude in the input filter.
        '''
        ### Now, AB mag only ###

        objflux = self.flux(filt)
        stdspec = ((3631*u.Jy).cgs * 3e10 * u.cm / u.s \
                / ((filt.wavelength*filt.units[0]).cgs)**2).to(self.units[1])

        midstdspec = (np.append(stdspec,0) + np.append(0,stdspec))[1:-1]/2.0

        interval = np.abs(filt.interval)
        stdflux = np.sum(interval * midstdspec * filt.midvalue)\
                * filt.units[0] * self.units[1]

        mag = -2.5*np.log10(float(objflux/stdflux))
        return mag

    def to_obs(self, redshift=1e-5):
        if self.mode == 'OBS':
            print('This Spectrum is Already an Observed Spectrum')

        elif self.mode == 'LUMI':
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            lumi_dist = (cosmo.luminosity_distance(redshift)).cgs
            newvalue = self.value/(1+redshift)/4/pi/lumi_dist**2
            newwavelength = self.wavelength*(1+redshift)
            newunits = [self.units[0],self.units[1]/u.cm/u.cm]

            self.update(newwavelength, newvalue, newunits)

        elif self.mode == 'ABS':
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            lumi_dist = (cosmo.luminosity_distance(redshift)).cgs
            std_dist = (10 * u.pc).cgs
            newvalue = self.value /\
                    float((1+redshift)*(lumi_dist/std_dist)**2)
            newwavelength = self.wavelength * (1+redshift)
            newunits = [self.units[0], self.units[1]]

            self.update(newwavelength, newvalue, newunits)

        self.mode = 'OBS'

    def to_abs(self, redshift=1e-5):

        if self.mode == 'ABS':
            print('This Spectrum is Already an Absolute Spectrum')

        elif self.mode == 'OBS':
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            lumi_dist = (cosmo.luminosity_distance(redshift)).cgs
            std_dist = (10*u.pc).cgs
            newvalue = self.value*float((1+redshift)*(lumi_dist/std_dist)**2)
            newwavelength = self.wavelength/(1+redshift)
            newunits = [self.units[0], self.units[1]]

            self.update(newwavelength, newvalue, newunits)

        elif self.mode == 'LUMI':
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            std_dist = float((10*u.pc).cgs/u.cm)
            newvalue = self.value/float(4/pi/std_dist**2)
            newwavelength = self.wavelength
            newunits = [self.units[0],self.units[1]/u.cm/u.cm]

            self.update(newwavelength, newvalue, newunits)

        self.mode = 'ABS'

    def to_luminosity(self, redshift=1e-5):
        if self.mode == 'LUMI':
            print('This Spectrum is Already in Luminosity Mode')

        elif self.mode == 'OBS':
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            lumi_dist = float(((cosmo.luminosity_distance(redshift)).cgs)/u.cm)
            newvalue = self.value * (1 + redshift) * 4 * pi * lumi_dist**2
            newwavelength = self.wavelength / (1 + redshift)
            newunits = [self.units[0], self.units[1]*u.cm*u.cm]

            self.update(newwavelength, newvalue, newunits)

        self.mode='LUMI'

    def normalize(self, filt, mag, style='AB'):
        mag0 = self.magnitude(filt, style)
        scale = 10**(0.4*(mag0 - mag))

        self.value = self.value * scale


