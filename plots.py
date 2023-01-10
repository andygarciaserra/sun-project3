#Packages and constants
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.special import exp10
import mpl_scatter_density
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import LogStretch
from matplotlib import colors
DATADIR = 'data/'
FITSDIR = '/home/andy/Downloads/prac3solar/data/'
FIGDIR = 'figures/'





#Booleans to execute parts of code
plot1a = False
plot1b = False
plot1c = False
plot2a = False
plot2a_paper = False
plot2b = False
plot3b = True





##---TASKS---##

## 1A ##
if(plot1a):
    #Loading data from .txt
    data = np.loadtxt(DATADIR+'INT_WVL_1a.txt',delimiter='\t',skiprows=1)
    int_Tmax = data[:,0]
    wvl_ord = data[:,1]
    temp = np.loadtxt(DATADIR+'T.txt',delimiter='\t',skiprows=1)
    data = np.loadtxt(DATADIR+'pos_T_WVL_1a.txt',delimiter='\t',skiprows=1)
    pos_T = int(data[0])
    pos_wvl = int(data[1])

    #Plotting
    plt.figure(figsize=(10,6))
    m1, st1, b1 = plt.stem(wvl_ord, int_Tmax, basefmt='', markerfmt='.', \
                           label='g ('+r'$\lambda$'+', '+r'$T_{max}$'+')')
    m2, st2, b2 = plt.stem(wvl_ord[pos_wvl], int_Tmax[pos_wvl], \
                           label=r'$\lambda$= '+str('{:.2f}'.format(wvl_ord[pos_wvl]))+ \
                            ' '+r'$\AA$'+'\n'+r'$T_{max}=$'+'{:.2e}'.format(temp[pos_T])+ \
                            ' K', linefmt='-', markerfmt='.')
        #formatting
    plt.setp(m1, markersize=2, color='tab:blue')
    plt.setp(st1, linewidth=0.4, color='tab:blue')
    plt.setp(b1, color='tab:blue')
    plt.setp(m2, markersize=5, color='tab:orange', marker='D')
    plt.setp(st2, linewidth=1.5, color='tab:orange')
    plt.setp(b2, color='tab:orange')
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=15)
    plt.ylabel('g ('+r'$T_{max}$'+', 'r'$\lambda$'+', '+r'$n_e$'+')' \
               +' ['+r'$erg$'+' '+ r'$cm^3$'+' '+'/'+' '+'s '+r'$ \ st$'+']', fontsize=15)
    plt.legend(loc='upper left', fontsize=15)
    plt.xlim((160,180))
    plt.ylim((1e-31,2e-22))
    ax = plt.gca()
    ax.set_yscale('log')
    plt.title('Gain function, Fe IX, T = '+'{:.2e}'.format(temp[pos_T])+' K', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig(FIGDIR+'gain_function_max.png', dpi=200)
    plt.show()

## 1B ##
if(plot1b):
    #Loading data from .txt
    data = np.loadtxt(DATADIR+'gtotal_T_1b.txt',delimiter='\t',skiprows=1)
    g_total = data[:,0]
    temp = data[:,1]
    data = np.loadtxt(DATADIR+'pos_T_WVL_1a.txt',delimiter='\t',skiprows=1)
    pos_T = int(data[0])
    
    #Plotting
    plt.figure(figsize=(10,6))
    plt.plot(temp, g_total, label='Fe IX')

        #formatting
    plt.xlabel('T'+' [K]', fontsize=15)
    plt.ylabel(r'$g_{total}$'+' ('+r'$\lambda$'+', '+r'$n_e$'+')' \
               +' ['+r'$erg$'+' '+ r'$cm^3$'+' '+'/'+' '+'s '+r'$ \ st$'+']', fontsize=15)
    plt.vlines(temp[pos_T], 1e-43, 1e-23, colors='r', linestyles='--')
    plt.plot(temp[pos_T],g_total[pos_T],color='r',marker='D', label='g = '+ \
             '{:.2e}'.format(g_total[pos_T])+' '+r'$erg$'+' '+r'$cm^3$'+'/'+ \
             's '+r'$st$'+'\n'+'T = '+'{:.2e}'.format(temp[pos_T])+ \
             ' K')
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.ylim((1e-39,1e-21))
    plt.title('Aggregated gain function, Fe IX', fontsize=15)
    plt.legend(fontsize=13)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig(FIGDIR+'gain_total_feix.png', dpi=300)
    plt.show()

## 1C ##
if(plot1c):
    #Loading data from .txt
    N_E = np.loadtxt(DATADIR+'N_E.txt',delimiter='\t',skiprows=1)

    #Plotting
    plt.figure(figsize=(10,6))
    for i in range(4):
        data = np.loadtxt(DATADIR+'gtotal_T_'+'{:.1e}'.format(N_E[i])+'_1c.txt',delimiter='\t',skiprows=1)
        g_total = data[:,0]
        temp = data[:,1]
        plt.plot(temp, g_total, label=r'$n_e$'+'= '+'{:.1e}'.format(N_E[i])+' '+ \
                 r'$cm^{-3}$'+'\t'+'')

        #Formatting
    plt.xlabel('T'+' [K]', fontsize=15)
    plt.ylabel(r'$g_{total}$'+' ('+r'$\lambda$'+', '+r'$n_e$'+')' \
            +' ['+r'$erg$'+' '+ r'$cm^3$'+' '+'/'+' '+'s '+r'$ \ st$'+']', fontsize=15)
    plt.legend(fontsize = 15)
    plt.xlim((1e5,1e7))
    plt.ylim((5e-24,1.5e-23))
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.title('Aggregated gain function, Fe IX', fontsize=15)
    plt.savefig(FIGDIR+'gain_total_ne.png', dpi=300)
    plt.show()

## 2A ##
if(plot2a):
    #Loading data from .txt
    data = np.loadtxt(DATADIR+'z_lgT_2a.txt',delimiter='\t',skiprows=1)
    z = data[:,0]
    mean_lgT = data[:,1]
    meanplus = data[:,1] + data[:,2]
    meanminus = data[:,1] - data[:,2]
    meanplus2 = data[:,1] + 2*data[:,2]
    meanminus2 = data[:,1] - 2*data[:,2]
    data = np.loadtxt(DATADIR+'paper.txt',delimiter='\t',skiprows=1)
    T_paper = data[:,1]
    z_paper = data[:,0]

    #Loading data from .fits
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data
    z = hdu[1].data

    #Creating random distribution
    N_RAND = 10000                  #N of random 2D plane points
    xrandint = np.random.randint(0,767,N_RAND)
    yrandint = np.random.randint(0,767,N_RAND)
    T_all = np.empty(shape=[768,N_RAND])
    for i in range(len(z)):
        rand_z = np.zeros(N_RAND)
        for j in range(len(xrandint)):
            rand_z[j] = lgT[i,xrandint[j],yrandint[j]]
        T_all[i] = rand_z
    
        #Plotting
    plt.figure()
    plt.plot(z[0], exp10(T_all[0,0]), color='tab:blue', label='T (z)')
    plt.plot(z, exp10(T_all), color='tab:blue', alpha=0.01)
    plt.plot(z_paper, T_paper, color='tab:orange', label='Gudiksen et al',ls='--',lw=1)
    plt.plot(z, exp10(mean_lgT), 'k', lw=2, label='< T > (z) [K]')
    plt.plot(z, exp10(meanplus), 'r-.', lw=1)
    plt.plot(z, exp10(meanminus), 'r-.', lw=1)
    plt.plot(z, exp10(meanplus2), 'r--', lw=1, label='2'+r'$\sigma$'+' (< T >) (z)')
    plt.plot(z, exp10(meanminus2), 'r--', lw=1)
    plt.fill_between(z, exp10(meanplus), exp10(meanminus), alpha=0.5, \
                     facecolor='red', label=r'$\sigma$'+' (< T >) (z)')
        #formatting
    plt.xlabel('z [Mm]', fontsize=15)
    plt.ylabel('T [K]', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc='lower right', fontsize=12)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.savefig(FIGDIR+'plt2a.png',dpi=300)
    plt.show()

## 2A_paper ##
if (plot2a_paper):
    #Loading data from .txt
    data = np.loadtxt(DATADIR+'z_lgT_2a.txt',delimiter='\t',skiprows=1)
    z = data[:,0]
    mean_lgT = data[:,1]
    meanplus = data[:,1] + data[:,2]
    meanminus = data[:,1] - data[:,2]
    meanplus2 = data[:,1] + 2*data[:,2]
    meanminus2 = data[:,1] - 2*data[:,2]
    data = np.loadtxt(DATADIR+'paper.txt',delimiter='\t',skiprows=1)
    T_paper = data[:,1]
    z_paper = data[:,0]

    #Loading data from .fits
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data
    z = hdu[1].data

    #Creating random distribution
    N_RAND = 10000                  #N of random 2D plane points
    xrandint = np.random.randint(0,767,N_RAND)
    yrandint = np.random.randint(0,767,N_RAND)
    T_all = np.empty(shape=[768,N_RAND])
    for i in range(len(z)):
        rand_z = np.zeros(N_RAND)
        for j in range(len(xrandint)):
            rand_z[j] = lgT[i,xrandint[j],yrandint[j]]
        T_all[i] = rand_z
    
    #Plotting
    plt.figure()
    plt.plot(z[0], exp10(T_all[0,0]), color='tab:blue', label='T (z)')
    plt.plot(z, exp10(T_all), color='tab:blue', alpha=0.01)
    plt.plot(z_paper, T_paper, color='magenta', label='Gudiksen et al',ls='--')
    plt.plot(z, exp10(mean_lgT), 'k', lw=2, label='< T > (z) [K]')
        #formatting
    plt.xlabel('z [Mm]', fontsize=15)
    plt.ylabel('T [K]', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(loc='lower right', fontsize=12)
    ax = plt.gca()
    ax.set_yscale('log')
    plt.savefig(FIGDIR+'plt2a_paper.png',dpi=300)
    plt.show()

## 2B ##
if (plot2b):

    #Loading data from .txt
    data = np.loadtxt(DATADIR+'z_lgT_ne_2b.txt',delimiter='\t',skiprows=1)
    z = data[:,0]
    mean_lgT = data[:,1]
    meanTplus = data[:,1] + data[:,2]
    meanTminus = data[:,1] - data[:,2]
    meanTplus2 = data[:,1] + 2*data[:,2]
    meanTminus2 = data[:,1] - 2*data[:,2]
    meanNEplus = data[:,3] + data[:,4]
    meanNEminus = data[:,3] - data[:,4]
    meanNEplus2 = data[:,3] + 2*data[:,4]
    meanNEminus2 = data[:,3] - 2*data[:,4]

    #Loading data from .fits
    hdu = fits.open(FITSDIR+'lgtg_750.fits')
    lgT = hdu[0].data
    z = hdu[1].data
    hdu = fits.open(FITSDIR+'lgne_750.fits')
    lgne = hdu[0].data

    #Creating random distribution
    N_RAND = 10000                  #N of random 2D plane points
    xrandint = np.random.randint(0,767,N_RAND)
    yrandint = np.random.randint(0,767,N_RAND)
    T_all = np.empty(shape=[768,N_RAND])
    ne_all = np.empty(shape=[768,N_RAND])
    for i in range(len(z)):
        rand_z1 = np.zeros(N_RAND)
        rand_z2 = np.zeros(N_RAND)
        for j in range(len(xrandint)):
            rand_z1[j] = lgT[i,xrandint[j],yrandint[j]]
            rand_z2[j] = lgne[i,xrandint[j],yrandint[j]]
        T_all[i] = rand_z1
        ne_all[i] = rand_z2
    ne_all = ne_all-6

    #Plotting
    norm = ImageNormalize(vmin=0,vmax=1000, stretch=LogStretch())
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='scatter_density')
    dens = ax.scatter_density(exp10(ne_all),exp10(T_all), norm=norm, cmap='Greys', dpi=None)
        #formatting
    plt.xlabel(r'$n_e$'+' ['+r'$cm^{-3}$'+']', fontsize=15)
    plt.ylabel('T [K]', fontsize=15)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig(FIGDIR+'plt2b.png', dpi=300)
    plt.show()

if (plot3b):

    #Importing from fits
    N = 16
    data = fits.getdata(DATADIR+'I_2d_N'+str(N)+'.fits')
    plt.figure()
    plt.imshow(data, cmap='plasma', norm=colors.LogNorm())
    plt.xlabel('X [Mm]',fontsize=15)
    plt.ylabel('Y [Mm]',fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    cb = plt.colorbar(label='I /'+r'$I_{max}$')
    cb.ax.tick_params(labelsize=15)
    plt.tight_layout()
    plt.savefig(FIGDIR+'I_2D_N'+str(N)+'.png',dpi=300)
    plt.show()
