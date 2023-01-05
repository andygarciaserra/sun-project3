#Packages and constants
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
DATADIR = 'data/'





#Booleans to execute code sections:
bool1a = False
bool1b = False
bool1c = False  
bool2a = True





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
    N_E = [3e7, 3e8, 3e10]
    
    for i in range(3):
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
    print('doing 2a')