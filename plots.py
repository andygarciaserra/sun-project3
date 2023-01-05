#Packages and constants
import numpy as np
import matplotlib.pyplot as plt
DATADIR = 'data/'





#Booleans to execute parts of code
plot1a = False
plot1b = False
plot1c = False





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
                           label='g ('+r'$\lambda$'+')')
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
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']')
    plt.ylabel('g ('+r'$T_{max}$'+', 'r'$\lambda$'+', '+r'$n_e$'+')' \
               +' ['+r'$erg$'+' '+ r'$cm^3$'+' '+'/'+' '+'s '+r'$ \ st$'+']')
    plt.legend()
    plt.xlim((160,180))
    plt.ylim((1e-31,2e-22))
    ax = plt.gca()
    ax.set_yscale('log')
    plt.title('Gain function, Fe IX, T = '+'{:.2e}'.format(temp[pos_T])+' K')
    plt.show()

## 1B ##
if(plot1b):
    #Loading data from .txt
    data = np.loadtxt(DATADIR+'gtotal_T_1b.txt',delimiter='\t',skiprows=1)
    g_total = data[:,0]
    temp = data[:,1]
    #Plotting
    plt.figure(figsize=(10,6))
    plt.plot(temp, g_total, label='Fe IX')

        #formatting
    plt.xlabel('T'+' [K]')
    plt.ylabel(r'$g_{total}$'+' ('+r'$\lambda$'+', '+r'$n_e$'+')' \
               +' ['+r'$erg$'+' '+ r'$cm^3$'+' '+'/'+' '+'s '+r'$ \ st$'+']')
    plt.legend()
    #plt.xlim((160,180))
    #plt.ylim((1e-31,2e-22))
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.title('Aggregated gain function, Fe IX')
    plt.show()

## 1C ##
if(plot1c):
    #Loading data from .txt
    N_E = np.loadtxt(DATADIR+'N_E.txt',delimiter='\t',skiprows=1)

        #Plotting
    plt.figure(figsize=(10,6))
    for i in range(3):
        data = np.loadtxt(DATADIR+'gtotal_T_'+'{:.1e}'.format(N_E[i])+'_1c.txt',delimiter='\t',skiprows=1)
        g_total = data[:,0]
        temp = data[:,1]
        plt.plot(temp, g_total, label=r'$n_e$'+'= '+'{:.1e}'.format(N_E[i])+' '+r'$cm^{-3}$')

        #formatting
    plt.xlabel('T'+' [K]')
    plt.ylabel(r'$g_{total}$'+' ('+r'$\lambda$'+', '+r'$n_e$'+')' \
            +' ['+r'$erg$'+' '+ r'$cm^3$'+' '+'/'+' '+'s '+r'$ \ st$'+']')
    plt.legend()
    #plt.xlim((160,180))
    plt.ylim((5e-24,2e-23))
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.title('Aggregated gain function, Fe IX')
    plt.show()

