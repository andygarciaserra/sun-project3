#Packages and constants
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.special import exp10
from tqdm import tqdm
DATADIR = 'data/'
FITSDIR = '/home/andy/Downloads/prac3solar/data/'





#Booleans to execute code sections:
bool1a = False
bool1b = False
bool1c = False  
bool2a = False
bool2b = False
bool3a = False
bool3bc = True






##---TASKS---##

## 1A ##
if (bool1a):

        #Defining temperature and electron density
    temp = np.geomspace(5e4,5e7,46)
    N_E = 3e9
    edens = np.full(46,N_E)


        #Calculating gain function
    Fe_ix = ch.ion('fe_9', temperature=temp, eDensity=edens, \
                abundance='sun_coronal_1992_feldman_ext')
    Fe_ix.intensity()
    index_ord = np.argsort(Fe_ix.Intensity['wvl'])
    wvl_ord = Fe_ix.Intensity['wvl'][index_ord]
    intensities_ord = Fe_ix.Intensity['intensity'][:, index_ord]
        
        #Finding maximums
    pos_T, pos_wvl = np.unravel_index(intensities_ord.argmax(), \
                                    intensities_ord.shape)

        #Saving data in .txt to plot later
    full_array = np.stack([intensities_ord[pos_T,:], wvl_ord], axis=1)
    np.savetxt(DATADIR+'INT_WVL_1a.txt', full_array, header='Intensity\tWavelength', \
            delimiter='\t', comments='')
    np.savetxt(DATADIR+'T.txt', temp, header='Temperature', \
            delimiter='\t', comments='') 
    np.savetxt(DATADIR+'pos_T_WVL_1a.txt', [pos_T, pos_wvl], header='Intensity\tWavelength', \
            delimiter='\t', comments='')

## 1B ##
if (bool1b):
        #Defining temperature and electron density
    temp = np.geomspace(5e4,5e7,46)
    edens = np.full(46,N_E)

        #Calculating gain function
    Fe_ix = ch.ion('fe_9', temperature=temp, eDensity=edens, \
                abundance='sun_coronal_1992_feldman_ext')
    Fe_ix.intensity()
    index_ord = np.argsort(Fe_ix.Intensity['wvl'])
    wvl_ord = Fe_ix.Intensity['wvl'][index_ord]
    intensities_ord = Fe_ix.Intensity['intensity'][:, index_ord]

        #Calculating the aggregated gain function
    g_total = np.sum(intensities_ord, axis=1)

        #Saving data in .txt to plot later
    full_array = np.stack([g_total, temp], axis=1)
    np.savetxt(DATADIR+'gtotal_T_1b.txt', full_array, header='g_total\tTemperature', \
            delimiter='\t', comments='')

## 1C ##
if (bool1c):
        #Defining temperature and electron density
    temp = np.geomspace(5e4,5e7,46)
    N_E = [3e7, 3e8, 3e9, 3e10]
    
    for i in range(4):
        edens = np.full(46,N_E[i])

            #Calculating gain function
        Fe_ix = ch.ion('fe_9', temperature=temp, eDensity=edens, \
                    abundance='sun_coronal_1992_feldman_ext')
        Fe_ix.intensity()
        index_ord = np.argsort(Fe_ix.Intensity['wvl'])
        wvl_ord = Fe_ix.Intensity['wvl'][index_ord]
        intensities_ord = Fe_ix.Intensity['intensity'][:, index_ord]

            #Calculating the aggregated gain function
        g_total = np.sum(intensities_ord, axis=1)

            #Saving data in .txt to plot later
        full_array = np.stack([g_total, temp], axis=1)
        np.savetxt(DATADIR+'gtotal_T_'+'{:.1e}'.format(N_E[i])+'_1c.txt', \
                   full_array, header='g_total\tTemperature', \
                   delimiter='\t', comments='')
    np.savetxt(DATADIR+'N_E.txt', N_E, header='Electron density', delimiter='\t', comments='')

## 2A ##
if (bool2a):

        #Importing from .fits
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data
    z = hdu[1].data

        #Mean values for each z
    mean_lgT = np.zeros(768)
    stdev = np.zeros(768)

    for i in range(768):
        mean_lgT[i] = lgT[i,:,:].mean(axis=(0,1))
        stdev[i] = lgT[i,:,:].std(axis=(0,1))

        #Saving data in .txt to plot later:
    full_array = np.stack([z,mean_lgT,stdev], axis=1)
    np.savetxt(DATADIR+'z_lgT_2a.txt', full_array, header='z [Mm]\tlog(T) [K]', \
            delimiter='\t', comments='')

## 2B ##
if (bool2b):

        #Importing from .fits
    hdu = fits.open(FITSDIR+'lgne_750.fits')
    lgne = hdu[0].data
    z = hdu[1].data
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data

        #Mean values for each z
    mean_lgT = np.zeros(768)
    stdev_T = np.zeros(768)
    mean_ne = np.zeros(768)
    stdev_ne = np.zeros(768)

    for i in range(768):
        mean_lgT[i] = lgT[i,:,:].mean(axis=(0,1))
        stdev_T[i] = lgT[i,:,:].std(axis=(0,1))
        mean_ne[i] = lgT[i,:,:].mean(axis=(0,1))
        stdev_ne[i] = lgT[i,:,:].std(axis=(0,1))

        #Saving data in .txt to plot later:
    full_array = np.stack([z,mean_lgT,stdev_T,mean_ne,stdev_ne], axis=1)
    np.savetxt(DATADIR+'z_lgT_ne_2b.txt', full_array, header='z [Mm]\t<log(T)> \
               [K]\tSTD(<log(T)>) [K]\t<Ne> [cm^-3]\tSTD(<Ne> [cm^-3])', \
               delimiter='\t', comments='')

## 3A ##
if (bool3a):
        #Importing from .fits
    hdu = fits.open(FITSDIR+'lgne_750.fits')
    lgne = hdu[0].data
    z = hdu[1].data*1e8                 #Mm to cm
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data
    edens = exp10(lgne[:,100,100]-6)    #m^-3 to cm^-3
    temp = exp10(lgT[:,100,100])

        #Calculating gain function via Chianti
    Fe_xiv = ch.ion('fe_14', temperature=temp, eDensity=edens, \
                abundance='sun_coronal_1992_feldman_ext')
    Fe_xiv.intensity()
    index_ord = np.argsort(Fe_xiv.Intensity['wvl'])
    wvl_ord = Fe_xiv.Intensity['wvl'][index_ord]
    intensities_ord = Fe_xiv.Intensity['intensity'][:, index_ord]

        #Finding maximums
    pos_T, pos_wvl = np.unravel_index(intensities_ord.argmax(), \
                                    intensities_ord.shape)

        #Calculating the aggregated gain function at lam = 211A
    int_reduced = intensities_ord[:,pos_wvl-5:pos_wvl+5]
    g211 = np.sum(int_reduced, axis=1)
    for i in range(len(temp)):
        if (temp[i]<1e5 or edens[i]>1e11):
            g211[i] = 0
    integrand = (1/1.2)*np.multiply(g211,np.multiply(edens,edens))
    int = np.trapz(integrand,z)
    print(int)


## 3B ##
if (bool3bc):
    
        #Importing from .fits
    hdu = fits.open(FITSDIR+'lgne_750.fits')
    lgne = hdu[0].data
    z = hdu[1].data*1e8                 #Mm to cm
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data

        #Equally spacing 2d grid for N=4
    N = 16                               #number of columns (NxN grid)
    pos = np.empty(shape=[N])
    I_2d = np.empty(shape=[N,N])
    step = int(767/(N-1))
    for i in range(N):
        pos[i] = i*step
    pos = pos.astype(int)

        #Calculating each column and creating 2D intensity array
    for l in range(N):                  #loop for x-axis
        for j in range(N):              #loop for y-axis
            edens = exp10(lgne[:,pos[l],pos[j]]-6)    #m^-3 to cm^-3
            temp = exp10(lgT[:,pos[l],pos[j]])

                #Calculating gain function via Chianti
            Fe_xiv = ch.ion('fe_14', temperature=temp, eDensity=edens, \
                        abundance='sun_coronal_1992_feldman_ext')
            Fe_xiv.intensity()
            index_ord = np.argsort(Fe_xiv.Intensity['wvl'])
            wvl_ord = Fe_xiv.Intensity['wvl'][index_ord]
            intensities_ord = Fe_xiv.Intensity['intensity'][:, index_ord]

                #Finding maximums
            pos_T, pos_wvl = np.unravel_index(intensities_ord.argmax(), \
                                            intensities_ord.shape)

                #Calculating the aggregated gain function at lam = 211A
            int_reduced = intensities_ord[:,pos_wvl-5:pos_wvl+5]
            g211 = np.sum(int_reduced, axis=1)
            for k in range(len(temp)):
                if (temp[k]<1e5 or edens[k]>1e11):
                    g211[k] = 0
            integrand = (1/1.2)*np.multiply(g211,np.multiply(edens,edens))
            int = np.trapz(integrand,z)
            I_2d[l,j] = int
    
        #Printing to double check I_2d dimensions
    print(I_2d)

        #Saving to .fits to plot later
    hdu = fits.PrimaryHDU(I_2d)
    hdu.writeto(DATADIR+'I_2d_N'+str(N)+'.fits')