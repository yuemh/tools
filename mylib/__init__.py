#mylib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from math import *
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord,Angle
from astropy.cosmology import FlatLambdaCDM
import pkg_resources

dir_data = pkg_resources.resource_filename('mylib', 'data')
ext_dir = dir_data + '/extinction'
exe_dir = dir_data + '/executable'


#crossmatch two catalogs by astropy
def crossmatch(RA_list1, DEC_list1, RA_list2, DEC_list2, radius=2/3600.):
    if not len(RA_list1) == len(DEC_list1):
        raise ValueError(\
                'The length of RA_list1 and DEC_list1 is not the same')
    if not len(RA_list2) == len(DEC_list2):
        raise ValueError(\
                'The length of RA_list2 and DEC_list2 is not the same')

    RA_list1 = np.array(RA_list1)
    DEC_list1 = np.array(DEC_list1)

    RA_list2 = np.array(RA_list2)
    DEC_list2 = np.array(DEC_list2)

    catalog1 = SkyCoord(ra = RA_list1*u.degree, dec = DEC_list1*u.degree)
    catalog2 = SkyCoord(ra = RA_list2*u.degree, dec = DEC_list2*u.degree)

    idx, d2d, d3d = catalog1.match_to_catalog_sky(catalog2)

    distance = np.array(d2d/u.degree)

    allindex_i = np.array(range(len(RA_list1)), dtype=int)
    output_i = allindex_i[ distance < radius ]
    output_j = idx[distance < radius ]

    '''
    output_i = []
    output_j = []

    for i in range(len(RA_list1)):
        if (distance[i]<radius):
            output_i.append(i)
            output_j.append(idx[i])
    '''

    return [output_i,output_j]


def name_to_coord(name):
    J_index = name.find('J')
    start_index = J_index + 1

    mid_index = np.max([name.find('+'), name.find('-')])

    ra_str = name[start_index:mid_index]
    dec_str = name[mid_index:]

    ra_str = ra_str[:2] + ':'\
            + ra_str[2:4] + ':'\
            +ra_str[4:]
    dec_str = dec_str[:3] + ':'\
            + dec_str[3:5] + ':'\
            +dec_str[5:]

    return (ra_str, dec_str)


def coord_to_name(coord, style='full', header=''):
    dummy = 1

    ra_str, dec_str = coord
    new_ra_str = ra_str.replace(':','')
    new_dec_str = dec_str.replace(':','')

    if not new_dec_str[0] in '+-':
        new_dec_str = '+' + new_dec_str

    if style=='full':
        name = header + 'J' + new_ra_str + new_dec_str

    elif style=='abbr':
        name = header + 'J' + new_ra_str[:4] + new_dec_str[:5]

    else:
        print('Unknown style in coord_to_name(). Use style==full')
        name = header + 'J' + new_ra_str + new_dec_str

    return name


#trim a image by IDL commands
def trim_image(image_name, IDL_dir = '', python_dir = '',\
               operate_dir = '', output_dir = '',\
               paras = [], suffix = '_trim.fits'):

    #IDL and python code directory
    if len(IDL_dir) == 0:
        IDL_dir = '/home/minghao/Desktop/research/commonuse/idl_files'
    if len(python_dir) == 0:
        python_dir = os.getcwd()

    f = open(IDL_dir + '/script_filename.txt', 'w')
    ftrim = open(IDL_dir + '/script_filename_trimmed.txt', 'w')
    string = image_name
    string_trimmed = image_name+suffix

    if len(paras) == 5:
        string_param = 'nextension=' + str(paras[0]) + '\n' + \
                      'x0=' + str(paras[1]) + '\n' + \
                      'x1=' + str(paras[2]) + '\n' + \
                      'y0=' + str(paras[3]) + '\n' + \
                      'y1=' + str(paras[4]) + '\n'
    else:
        raise ValueError('Number of parameters must be 5, but '\
                          + str(len(paras)) + ' input')

    rest_part=\
'''
OPENR, lun,'script_filename.txt', /GET_LUN  ; read file name
filename = '' ; Read one line at a time, saving the result into array
line = ''
row=0
WHILE NOT EOF(lun) DO BEGIN 
  READF, lun, line 
  filename = [filename, line] 
  row=row+1
ENDWHILE
FREE_LUN, lun ; Close the file and free the file unit

OPENR, lun,'script_filename_trimmed.txt', /GET_LUN  ; read trimed file name
filename_trimed = '' ; Read one line at a time, saving the result into array
line = ''
row=0
WHILE NOT EOF(lun) DO BEGIN 
  READF, lun, line 
  filename_trimed = [filename_trimed, line] 
  row=row+1
ENDWHILE
FREE_LUN, lun ; Close the file and free the file unit

print, row
print, filename
print, filename_trimed

for num=1,row do begin
    image=mrdfits(filename[num],nextension,header)
    hextract,image,header,x0,x1,y0,y1
    writefits,filename_trimed[num],image,header
endfor
end
'''
    string_IDL_command = string_param + rest_part
    fcommand = open(IDL_dir+'/IDL_journal.pro','w')
    fcommand.write(string_IDL_command)
    fcommand.close()

    f.write(operate_dir + '/' + string)
    ftrim.write(output_dir + '/' + string_trimmed)
    f.close()
    ftrim.close()

    os.chdir(IDL_dir)

    os.system('idl ' + IDL_dir + '/idlcommand.pro')

    os.system('rm ' + IDL_dir + '/script_filename.txt')
    os.system('rm ' + IDL_dir + '/script_filename_trimmed.txt')

    os.chdir(python_dir)
    return string_trimmed


def extinction_curve(wavelength, Rv=3.1):
    ### Wavelength in unit AA ###

    x = 10000. / wavelength
    if 0.3 < x <= 1.1:
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61

    elif 1.1 < x <= 3.3:
        y = x - 1.82
        a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + \
                0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4\
                - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7

    elif 3.3 < x <= 8:
        if 3.3 < x <= 5.9:
            Fa = 0
            Fb = 0
        elif 5.9 < x <= 8:
            Fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
            Fb = 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3

        a = 1.752-0.316*x - 0.104/((x-4.67)**2+0.341) + Fa
        b = -3.09 + 1.825*x + 1.206/((x-4.62)**2+0.263) + Fb

    elif 8 < x <= 10:
        a = -1.073 -0.628*(x-8) + 0.137*(x-8)**2 - 0.07*(x-8)**3
        b = 13.67 + 4.257*(x-8) - 0.420*(x-8)**2 + 0.374*(x-8)**3

    else:
        raise ValueError('Wavelength Out of Range')

    return a + b / Rv


def read_ebv(ra,dec,dustfile=ext_dir + '/lambda_zea_sfd_ebv.fits'):
    c = SkyCoord(ra=ra, dec=dec, unit='deg', frame='fk5')
    galactic = c.galactic

    l, b = galactic.l / u.deg, galactic.b / u.deg

#    print(l,b)
    if b > 0:
        ext = 1
    else:
        ext = 2

    dust = fits.open(dustfile)[ext]
    w = wcs.WCS(dust.header)

    pixcrd = w.all_world2pix([[l,b]],1)

    x, y = pixcrd[0]
    print(x, y)
    xmax, ymax = dust.data.shape
    x0, y0 = int(x), int(y)
    x1, y1 = x0+1, y0+1

    e00, e01, e10, e11 = \
            dust.data[y0,x0], dust.data[y1,x0], dust.data[y0,x1], dust.data[y1,x1]
    ky = e01 - e00
    kx = e10 - e00

    e = e00 + (x-x0)*kx + (y-y0)*ky

    return e


def read_catalog_sextractor(filename):
    f=open(filename)

    namelist=[]

    headflag=1
    while headflag:
        string=f.readline()
        if len(string)==0:
            break
        if not string[0]=='#':
            headflag=0
        else:
            name=string.split()[2]
            namelist.append(name)

    f.close()

    df=pd.read_csv(filename,delim_whitespace=1,comment='#',names=namelist)
    return df.copy()


def DS9(filename, extension=-1, options=''):
    if not type(extension)==int:
        raise ValueError('The extension number must be a non-negative integer')

    if extension > 0:
        extension_str = '['+str(extension)+']'
    else:
        extension_str=''

    os.system(exe_dir + '/ds9 ' + filename + extension_str + options)
    os.system('exit')

