#! /usr/bin/env python3
summary = """plot_tools.py

Description:
  Module for plot various concepts in a built dataframe with dframe_tools.py

Prerequisite:
  numpy
  pandas
  matplotlib

Optional files:
  include/notebook.mplstyle
  include/aps.mplstyle
"""

# phase_tool.py
# Sang Wook Kim
# 06.02.2022

import numpy as np
from numpy import pi as π
import math
import pandas as pd
import matplotlib.pyplot as plt
import re,glob,os
from collections import defaultdict
import dgutils

# if os.path.exists('include/notebook.mplstyle') and os.path.exists('include/aps.mplstyle'):
#     # plot style
#     plot_style = {'notebook':'include/notebook.mplstyle','aps':'include/aps.mplstyle'}
#     plt.style.reload_library()
#     plt.style.use(plot_style['aps'])
#     figsize = plt.rcParams['figure.figsize']
#     plt.rcParams['text.latex.preamble'] = f'\input{{{os.getcwd()}/include/texheader}}'

#     colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#     print("plot style is loaded")
# else:
#     print("plot style don't exist")

if os.path.exists('../include/sans_NC.mplstyle') and os.path.exists('../include/notebook.mplstyle'):
    ### plot style
    plot_style = {'notebook':'../include/notebook.mplstyle','sans':'../include/sans_NC.mplstyle'}
    plt.style.reload_library()
    plt.style.use(plot_style['sans'])
    figsize = plt.rcParams['figure.figsize']

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # setup a possible custom font path
    from matplotlib import font_manager
    def what_font_path (filename):
        for p in font_manager.findSystemFonts(fontpaths=None, fontext="ttf"):
            if p.find(filename) != -1:
                return p

    font_path,bold_font_path = '.','.'
    if 'LOCAL_FONT_PATH' in os.environ:
        font_path = os.environ['LOCAL_FONT_PATH'] + os.path.sep + 'HelveticaNeue/HelveticaNeue-Light-08.ttf'
        bold_font_path = os.environ['LOCAL_FONT_PATH'] + os.path.sep + 'HelveticaNeue/HelveticaNeue-Bold-02.ttf'
    else:
        # local path (custom)
        font_path = what_font_path('HelveticaNeue-Light-08')
        bold_font_path = what_font_path('HelveticaNeue-Bold-02')
else:
    print("plot style don't exist")

# -----------------------------------------------------------------------------
def Σ(σ,q):
    '''Compute the Σ function needed for linear fits.'''
    return np.sum(q/σ**2)
# -----------------------------------------------------------------------------
def get_a(x,y,σ):
    '''Get the χ^2 best fit value of a0 and a1.'''
    # This fails when σ has zero.
    epsilon = np.finfo(np.float32).eps
    σ = np.where(σ == 0, epsilon, σ)
    
    # Get the individual Σ values
    Σy,Σx,Σx2,Σ1,Σxy = Σ(σ,y),Σ(σ,x),Σ(σ,x**2),Σ(σ,np.ones(x.size)),Σ(σ,x*y)

    # the denominator
    D = Σ1*Σx2 - Σx**2

    # compute the best fit coefficients
    a = np.array([Σy*Σx2 - Σx*Σxy,Σ1*Σxy - Σx*Σy])/D

    # Compute the error in a
    aErr = np.array([np.sqrt(Σx2/D),np.sqrt(Σ1/D)])

    return a,aErr
# -----------------------------------------------------------------------------
def linear(x,a):
    '''Return a polynomial of order'''
    return a[0] + a[1]*x
# -----------------------------------------------------------------------------
def datadic(df):
    '''Make subset dataframes for each strain value and chemical potential'''
    sliced_dict = {}

    for strain in df['strain'].unique():
        for mu in df['mu'].unique():
            subset = df[(df['strain'] == strain)&(df['mu'] == mu)]
            if subset.empty:
                continue
            subset = subset.sort_values('T', ascending=False)
            sliced_dict[str(strain)+','+str(mu)] = subset
    
    print(sliced_dict.keys())
    return sliced_dict
# -----------------------------------------------------------------------------
def esti_array(subdf):
    '''Put dataframe and return dictionary of features below'''
    '''Features list:
        'Tset': list of unique temperature values
        'totN': total number of adsorption sites based on graphene
        'n': number of particles
        'nerr': standard error of n
        'yarray': superfluid fraction
        'yerrarray': maximum error of binning error analysis
        'rhosarray': superfluid density (fraction*n/area)
        'rhoserrarray': yerrarray*n/area
    '''
    result = {}
    Tsets = sorted(list(set(subdf['T'])))
    yarray = []
    yerrarray = []
    narray = []
    nerrarray = []
    Narray = []
    for Ti in Tsets:
        sublst = subdf[subdf['T'] == Ti].sort_values('totN',ascending=False)
        yarray.append(sublst['rhos'])
        yerrarray.append(sublst['rhoserr'])
        narray.append(sublst['n'])
        nerrarray.append(sublst['nerr'])
        Narray.append(sublst['totN'])
    yarray = np.asarray(yarray)
    yerrarray = np.asarray(yerrarray)
    narray = np.asarray(narray)
    nerrarray = np.asarray(nerrarray)
    Narray = np.asarray(Narray)

    rhoarray = narray/Narray/(math.sqrt(3)/4*6*1.42)*10*6.7 # 10^15 (cm^-2) * 10^-24 (g) = 10* 10^(-10)
    rhosarray = yarray*rhoarray
    rhoserrarray = yerrarray*rhoarray
    
    result['Tset'] = Tsets
    result['totN'] = Narray
    result['n'] = narray
    result['nerr'] = nerrarray
    result['yarray'] = yarray
    result['yerrarray'] = yerrarray
    result['rhosarray'] = rhosarray
    result['rhoserrarray'] = rhoserrarray
    
    return result
# -----------------------------------------------------------------------------
def esti_array_multi(subdf):
    '''same with esti_array, but choose mean value for estimator and max value for error
    when we have more than one data row with the same configurations'''
    result = {}
    Tsets = sorted(list(set(subdf['T'])))
    yarray = []
    yerrarray = []
    narray = []
    nerrarray = []
    Narray = []
    for Ti in Tsets:
        sublst = subdf[subdf['T'] == Ti].sort_values('totN',ascending=False)
        yarray.append(sublst.groupby('totN')['rhos'].mean())
        yerrarray.append(sublst.groupby('totN')['rhoserr'].max())
        narray.append(sublst.groupby('totN')['n'].mean())
        nerrarray.append(sublst.groupby('totN')['nerr'].max())
        Narray.append(sublst.groupby('totN')['totN'].mean())
        
    yarray = np.asarray(yarray)
    yerrarray = np.asarray(yerrarray)
    narray = np.asarray(narray)
    nerrarray = np.asarray(nerrarray)
    Narray = np.asarray(Narray)

    rhoarray = narray/Narray/(math.sqrt(3)/4*6*1.42)*10*6.7 # 10^15 (cm^-2) * 10^-24 (g) = 10* 10^(-10)
    rhosarray = yarray*rhoarray
    rhoserrarray = yerrarray*rhoarray
    
    result['Tset'] = Tsets
    result['totN'] = Narray
    result['n'] = narray
    result['nerr'] = nerrarray
    result['yarray'] = yarray
    result['yerrarray'] = yerrarray
    result['rhosarray'] = rhosarray
    result['rhoserrarray'] = rhoserrarray
    
    return result
# -----------------------------------------------------------------------------
def nparallel(subdf):
    '''return number of copies with same configurations in the subset dataframes
    if the total # of data is n times total # of (T, totN) combination.
    return None if it's something else'''
    Tsets = subdf['T'].unique()
    Nsets = subdf['totN'].unique()
    intg = len(subdf)/(len(Tsets)*len(Nsets))
#     return intg
    if intg.is_integer():
        return int(intg)
    else:
        return None
# -----------------------------------------------------------------------------
def plot_frac(result):
    '''using result of esti_array, plot n/totN vs totN for each temperature'''
    '''also returns mean of it and maximum of error'''
    x = result['totN'][0]
    x = np.asarray([1/(item) for item in x])
#     x = np.asarray([item for item in x])
#     x = np.asarray([1/np.sqrt(item) for item in x])
#     alist = []
#     aerrlist = []
    ylist = []
    yerrlst = []
    nfig = len(result['Tset'])
    nrow = (nfig+1)//2
    fitplot, axs = plt.subplots(nrows=nrow, ncols=2, 
                                sharex=False, figsize = [9,2*nrow])
    
    # Defining custom 'xlim' and 'ylim' values.
    custom_xlim = (0, max(x)*1.1)
    custom_ylim = (0, np.amax(result['n']/result['totN'])*1.1)

    # Setting the values for all axes.
    plt.setp(axs, xlim=custom_xlim, ylim=custom_ylim)
    
    if len(result['Tset'])<=2:
        
        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x = np.asarray([1/(item) for item in x])
            y = result['n'][i]/result['totN'][i]
            σ = result['nerr'][i]

            # peform the fits
            a1,a1_err = get_a(x,y,σ)
            
            # plot the data
            axs[i].errorbar(x, y, yerr = σ,color=colors[1], fmt='--o',label=f'T={tag}')

            # plot the fit results
            fx = np.linspace(0,max(x)*1.1,20)
            axs[i].plot(fx,a1[0]+a1[1]*fx, color=colors[1], linewidth=1.5, zorder=0, label=f'T={tag} fit')
            
            axs[i].set_title(f'T={tag}')
            axs[i].set_xlabel('1/N')
            axs[i].set_ylabel('filling')

#             ylist.append(y)
#             yerrlst.append(σ)

            ylist.append(a1)
            yerrlst.append(a1_err)
        
    else:

        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x = np.asarray([1/(item) for item in x])
            y = result['n'][i]/result['totN'][i]
            σ = result['nerr'][i]

            # peform the fits
            a1,a1_err = get_a(x,y,σ)
            
            # plot the data
            axs[i//2,i%2].errorbar(x, y, yerr = σ, color=colors[1], fmt='--o',label=f'T={tag}')

            #plot the fit results
            fx = np.linspace(0,max(x)*1.1,20)
            axs[i//2,i%2].plot(fx,a1[0]+a1[1]*fx, color=colors[1], linewidth=1.5, zorder=0, label=f'T={tag} fit')
    
            axs[i//2,i%2].set_title(f'T={tag}')
            axs[i//2,i%2].set_xlabel('1/N')
            axs[i//2,i%2].set_ylabel('filling')

#             ylist.append(y)
#             yerrlst.append(σ)

            ylist.append(a1)
            yerrlst.append(a1_err)

    fitplot.tight_layout()
    
#     return np.mean(ylist), np.amax(yerrlst)
    return ylist, yerrlst, fitplot
# -----------------------------------------------------------------------------
def non_plot_frac(result):
    '''using result of esti_array, plot n/totN vs totN for each temperature'''
    '''also returns mean of it and maximum of error'''
    x = result['totN'][0]
    x = np.asarray([1/(item) for item in x])
#     x = np.asarray([item for item in x])
#     x = np.asarray([1/np.sqrt(item) for item in x])
#     alist = []
#     aerrlist = []
    ylist = []
    yerrlst = []
    epsilon = np.finfo(np.float32).eps

    for i in range(len(result['Tset'])):
        tag = str(result['Tset'][i])
        x = result['totN'][i]
        x = np.asarray([1/(item) for item in x])
        y = result['n'][i]/result['totN'][i]
        σ = result['nerr'][i]

        # peform the fits
        a1,a1_err = get_a(x,y,σ)

#             ylist.append(y)
#             yerrlst.append(σ)

        ylist.append(a1)
        yerrlst.append(a1_err)

    
#     return np.mean(ylist), np.amax(yerrlst)
    return ylist, yerrlst
# -----------------------------------------------------------------------------
def plot_superfrac(result):
    '''using result of esti_array, plot sf fraction vs 1/sqrt(totN) for each temperature'''
    '''also returns fitting parameters'''
    x = result['totN'][0]
    x = np.asarray([1/np.sqrt(item) for item in x])
    alist = []
    aerrlist = []
    nfig = len(result['Tset'])
    nrow = (nfig+1)//2
    fitplot, axs = plt.subplots(nrows=nrow, ncols=2, 
                                sharex=False, figsize = [9,2*nrow])
    
    # Defining custom 'xlim' and 'ylim' values.
    custom_xlim = (0, max(x)*1.1)
    custom_ylim = (0, np.amax(result['yarray'])*1.1)

    # Setting the values for all axes.
    plt.setp(axs, xlim=custom_xlim, ylim=custom_ylim)
    
    if len(result['Tset'])<=2:
        
        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x = np.asarray([1/np.sqrt(item) for item in x])
            y = result['yarray'][i]
            σ = result['yerrarray'][i]

            # peform the fits
            a1,a1_err = get_a(x,y,σ)

            # plot the data
            axs[i].errorbar(x,y,yerr = σ,color=colors[0], fmt='--o',label=f'T={tag}')

            # plot the fit results
            fx = np.linspace(0,max(x)*1.1,20)

            axs[i].plot(fx,a1[0]+a1[1]*fx, color=colors[0], linewidth=1.5, zorder=0, label=f'T={tag} fit')
            axs[i].set_title(f'T={tag}')
            axs[i].set_xlabel('sqrt(1/N)')
            axs[i].set_ylabel('SF fraction')

            alist.append(a1)
            aerrlist.append(a1_err)
        
    else:

        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x = np.asarray([1/np.sqrt(item) for item in x])
            y = result['yarray'][i]
            σ = result['yerrarray'][i]

            # peform the fits
            a1,a1_err = get_a(x,y,σ)

            # plot the data
            axs[i//2,i%2].errorbar(x,y,yerr = σ,color=colors[0], fmt='--o',label=f'T={tag}')

            # plot the fit results
            fx = np.linspace(0,max(x)*1.1,20)

            axs[i//2,i%2].plot(fx,a1[0]+a1[1]*fx, color=colors[0], linewidth=1.5, zorder=0, label=f'T={tag} fit')
            axs[i//2,i%2].set_title(f'T={tag}')
            axs[i//2,i%2].set_xlabel('sqrt(1/N)')
            axs[i//2,i%2].set_ylabel('SF fraction')

            alist.append(a1)
            aerrlist.append(a1_err)

    fitplot.tight_layout()
    
    return alist, aerrlist, fitplot
# -----------------------------------------------------------------------------
def non_plot_superfrac(result):
    '''using result of esti_array, plot sf fraction vs 1/sqrt(totN) for each temperature'''
    '''also returns fitting parameters'''
    x = result['totN'][0]
    x = np.asarray([1/np.sqrt(item) for item in x])
    alist = []
    aerrlist = []
        
    for i in range(len(result['Tset'])):
        tag = str(result['Tset'][i])
        x = result['totN'][i]
        x = np.asarray([1/np.sqrt(item) for item in x])
        y = result['yarray'][i]
        σ = result['yerrarray'][i]

        # peform the fits
        a1,a1_err = get_a(x,y,σ)

        alist.append(a1)
        aerrlist.append(a1_err)
    
    return alist, aerrlist
# -----------------------------------------------------------------------------
def plot_only(result):
    '''using result of esti_array, plot n/totN vs totN for each temperature'''
    '''also returns mean of it and maximum of error'''
    _x = result['totN'][0]
    _x1 = np.asarray([1/(item) for item in _x])
    _x2 = np.asarray([1/np.sqrt(item) for item in _x])
#     x = np.asarray([item for item in x])
#     x = np.asarray([1/np.sqrt(item) for item in x])
    alist = []
    aerrlist = []
    ylist = []
    yerrlst = []
    nfig = len(result['Tset'])
    if nfig <= 1:
        print("empty or single data")
        return None
    nrow = (nfig)//2+1
    plt.style.reload_library()
    with plt.style.context(plot_style['sans']):
        fig = plt.figure(figsize = (3.4646*2, 2.14122*nrow))

        ax1 = fig.add_subplot(nrow, 2, 1)
        axs = []
        
#         fitplot, axs = plt.subplots(nrows=nrow, ncols=2, 
#                                     sharex=True, figsize = (3.4646*2, 2.14122*nrow) )
        dgutils.fonts.set_custom_font(font_path)

        # Defining custom 'xlim' and 'ylim' values.
        custom_xlim1 = (0, max(_x1)*1.1)
        custom_ylim1 = (0, np.amax(result['n']/result['totN'])*1.1)

            # Defining custom 'xlim' and 'ylim' values.
        custom_xlim2 = (0, max(_x2)*1.1)
        custom_ylim2 = (0, np.amax(result['yarray'])*1.1)

        # Setting the values for all axes.
    #     plt.setp(axs, xlim=custom_xlim2, ylim=custom_ylim2)


        ax1.set_xlabel('1/N')
        ax1.set_ylabel('filling fraction')
        ax1.set_xlim(custom_xlim1)
        ax1.set_ylim(custom_ylim1)

        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x1 = np.asarray([1/(item) for item in x])
            x2 = np.asarray([1/np.sqrt(item) for item in x])
            y1 = result['n'][i]/result['totN'][i]
            y2 = result['yarray'][i]
            σ1 = result['nerr'][i]
            σ2 = result['yerrarray'][i]

            # peform the fits
            a1,a1_err = get_a(x1,y1,σ1)
            a1s,a1s_err = get_a(x2,y2,σ2)

            # plot the data
            ax1.errorbar(x1, y1, yerr = σ1, color=colors[i], fmt='.', ms = 5)

            # plot the fit results
            fx1 = np.linspace(0,max(x1)*1.1,20)
            ax1.plot(fx1,a1[0]+a1[1]*fx1, linewidth=1.5, color=colors[i], zorder=0, label=f'T={tag} fit')

            # plot the data
            if i<=1:
                axs.append( fig.add_subplot(nrow, 2, i+2) )
            else:
                axs.append( fig.add_subplot(nrow, 2, i+2, ) ) #sharex = axs[i-2]
                           
            
                
            axs[i].errorbar(x2,y2,yerr = σ2,color=colors[i], fmt='.',label=f'T={tag}', ms = 5)

            # plot the fit results
            fx2 = np.linspace(0,max(x2)*1.1,20)

            axs[i].plot(fx2,a1s[0]+a1s[1]*fx2, color=colors[i], linewidth=1.5, zorder=0, label=f'T={tag} fit')
#             axs[i].set_xlabel('sqrt(1/N)')
#             axs[i].set_ylabel('SF fraction')
            axs[i].set_xlim(custom_xlim2)
            axs[i].set_ylim(custom_ylim2)
            axs[i].legend()
            
            if i%2 == 1 or i==0:
                axs[i].set_ylabel('SF fraction')
            else:
                axs[i].axes.yaxis.set_ticklabels([])
            if i < len(result['Tset'])-2:
                axs[i].axes.xaxis.set_ticklabels([])

        axs[i-1].set_xlabel('sqrt(1/N)')
        axs[i].set_xlabel('sqrt(1/N)')

        ax1.legend()
        fig.tight_layout()

    return fig
# -----------------------------------------------------------------------------
def plot_superdens(result):
    '''using result of esti_array, plot sf density vs 1/sqrt(totN) for each temperature'''
    '''also returns fitting parameters'''
    x = result['totN'][0]
    x = np.asarray([1/np.sqrt(item) for item in x])
    alist = []
    aerrlist = []
    nfig = len(result['Tset'])
    nrow = (nfig+1)//2
    fitplot, axs = plt.subplots(nrows=nrow, ncols=2, 
                                sharex=False, figsize = [9,2*nrow])
    
    # Defining custom 'xlim' and 'ylim' values.
    custom_xlim = (0, max(x)*1.1)
    custom_ylim = (0, np.amax(result['rhosarray'])*1.1)

    # Setting the values for all axes.
    plt.setp(axs, xlim=custom_xlim, ylim=custom_ylim)
    
    if len(result['Tset'])<=2:
        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x = np.asarray([1/np.sqrt(item) for item in x])
            y = result['rhosarray'][i]
            σ = result['rhoserrarray'][i]

            # peform the fits
            a1,a1_err = get_a(x,y,σ)

            # plot the data
            axs[i].errorbar(x,y,yerr = σ,color=colors[0], fmt='--o',label=f'T={tag}')

            # plot the fit results
            fx = np.linspace(0,max(x)*1.1,20)

            axs[i].plot(fx,a1[0]+a1[1]*fx, color=colors[0], linewidth=1.5, zorder=0, label=f'T={tag} fit')
            axs[i].set_title(f'T={tag}')
            axs[i].set_xlabel('sqrt(1/N)')
            axs[i].set_ylabel(r'SF density$\times 10^{9} \,$ ($\text{g cm}^{-2}$)')

            alist.append(a1)
            aerrlist.append(a1_err)
        
    else:

        for i in range(len(result['Tset'])):
            tag = str(result['Tset'][i])
            x = result['totN'][i]
            x = np.asarray([1/np.sqrt(item) for item in x])
            y = result['rhosarray'][i]
            σ = result['rhoserrarray'][i]

            # peform the fits
            a1,a1_err = get_a(x,y,σ)

            # plot the data
            axs[i//2,i%2].errorbar(x,y,yerr = σ,color=colors[0], fmt='--o',label=f'T={tag}')

            # plot the fit results
            fx = np.linspace(0,max(x)*1.1,20)

            axs[i//2,i%2].plot(fx,a1[0]+a1[1]*fx, color=colors[0], linewidth=1.5, zorder=0, label=f'T={tag} fit')
            axs[i//2,i%2].set_title(f'T={tag}')
            axs[i//2,i%2].set_xlabel('sqrt(1/N)')
            axs[i//2,i%2].set_ylabel(r'SF density$\times 10^{9} \,$ ($\text{g cm}^{-2}$)')

            alist.append(a1)
            aerrlist.append(a1_err)

    fitplot.tight_layout()
    
    return alist, aerrlist

# -----------------------------------------------------------------------------
usage = '''How to use:
1. Make classified data dictionary with datadic(dataframe)
2. The dictionary values are subdf and take a look its keys
3. use esti_array_multi(subdf) to get dictionary of features, [result]
4. plot the [result] with
   mean, max = plot_frac(result)
   alist, aerrlist = plot_superfrac(result)
   alist, aerrlist = plot_superdens(result)
'''
# -----------------------------------------------------------------------------
def help():
    print(summary)
    print(usage)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("This is module")
    print(summary)
    print(usage)
