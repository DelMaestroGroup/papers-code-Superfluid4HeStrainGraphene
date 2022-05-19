#! /usr/bin/env python3
summary = """phase_tool.py

Description:
  Module for reading and building dataframe to analyze phase of atoms

Prerequisite:
  graphenetools-py : pip install graphenetools-py
  erroranalysis-py : pip install erroranalysis-py
  pandas
  scipy
"""

# phase_tool.py
# Sang Wook Kim
# 10.22.2021

import numpy as np
from graphenetools import gt
import sys,importlib
# from dgutils import colors as colortools
from numpy import pi as π
import math
# import pimcscripts.pimchelp as pimchelp
import pandas as pd
import os
import glob
import erroranalysis as ea
from scipy.constants import k as k_b

# -----------------------------------------------------------------------------
def bin(MC):
    
    # minimum number of MC bins required
    min_bin = 32

    # initialize B to MC data
    B = np.copy(MC)

    # Resize if 1D array
    if B.ndim == 1:
        B.resize(B.shape[0],1)
    
    if len(MC) < 32:
        return 0
    
    try:
        # Define number of binning levels
        Nl = int(math.floor(math.log(B.shape[0]/min_bin,2))+1)

    except:
        print("Not enough bins. Need {} bins and have {} bins. Setting binning level to 1.".format(min_bin,B.shape[0]))
        Nl = 1

    # initialize D
    D = np.zeros((Nl,B.shape[1]))
    
    # First level of binning is raw data
    D[0,:] = np.std(MC,0)/math.sqrt(B.shape[0]-1)
    
    # Binning loop over levels l
    for l in range(1,Nl):
        
        # Bin pairs of bins: if odd # of bins, truncate first bin
        if ((B.shape[0] % 2) == 0):
            B = (B[::2,:]+ B[1::2,:])/2
        else:
            B = (B[1::2,:]+ B[2::2,:])/2
        
        # Error at level l
        D[l,:] = np.std(B,0)/math.sqrt(B.shape[0]-1)
    
    return D
# -----------------------------------------------------------------------------
# !!!volume and k_B need to be review carefully!!!(physical values are not reflected yet)
def compressibility(N,N2,volume=1.0,k_B=1.0,T=1.0):
    pf = volume/(k_B*T)
    return pf*(((N2 - (N**2)))/(N**2))
# -----------------------------------------------------------------------------
def stats(data):
    '''Return the average and standard error in data. '''
    ave = np.average(data)
    ave2 = np.average(data*data)
    err = np.sqrt(np.abs(ave2-ave**2)/(1.0*data.size-1.0) ) 
    return ave,err
# -----------------------------------------------------------------------------
def find_total_Nsites(Lx,Ly,strain):
    '''find total # of absorption sites based on box size'''
    # if you didn't generate the right box size, this can be wrong
    a = 1.42*strain
    nx = round(Lx/(a*math.sqrt(3)))
    ny = round(Ly/(1.5*a))
    return nx*ny
# ----------------------------------------------------------------------
def get_estimator_names(base_dir,pimcid,verbose=False):
    '''Return a list of estimator and log file names.'''
    est_name = ['log', 'estimator', 'super',
                'linedensity', 'planeavedensity',
                'energy','position','scom', 
                'lineardensity']

    file_names = {}
    for est in est_name:
        name = f"{base_dir}/*-{est}-*-{pimcid}.dat" 
        file_name = glob.glob(name)
        if file_name:
            file_names[est] = file_name[0]
        else:
            if verbose:
                print(f"{name} doesn't exist")
    return file_names
# -----------------------------------------------------------------------------
class Filewrongformat(Exception):
    def __init__(self,ref):
        Exception.__init__(self,"There is something wrong in "+ref)
# -----------------------------------------------------------------------------
def df_pimc(base_dir, skip = 0, noenergy = 0, replace = 0):
    '''
    base_dir: base directory
    skip: if integer, skip the (skip) number of lines. if 0=<skip<1, skip the ratio. Default: 0
    noenergy: if target estimator file has no energy coloumns (K,V,V_ext,V_int,E,E_mu,K/N,V/N,E/N), set True
    '''
    
    """output data columns - 
    'id': pimc id
    'strain': iso-strain
    'mu': chemical potential
    'T': temperature
    'filedic': file_names dictionary ex) 'log', 'estimator', 'super'
    'n': number of particle
    'nerr': standard error of n
    'totN': total number of adsorption sites based on graphene
    'kap': compressibility
    'kaperr': standard error of kap
    'rhos': superfulid fraction
    'rhoserr': maximum error of binning error analysis
    'estsize': number of rows in estimator file
    'supsize': number of rows in superfulid file
    'boxdims': size of box,
    'com_ave': average of order parameter,
    'com_err': error of order parameter,
    """
    
    # iterate through all file
    file_ids = []
    file_strain = []
    file_mu = []
    file_T = []

    for file in os.listdir(base_dir):
        # Check whether file is in right format
        if file.endswith(".dat") and file.find("log")!=-1:
            file_path = f"{base_dir}/{file}"
            #Open it and go thorugh line by line
            with open(file_path) as infile:
                for line in infile:
                    if line.startswith("# PIMCID"):
                        file_ids.append(line[10:-1])
                    if line.startswith("Command String"):
                        try:
                            idx = line.find("graphene_isotropic_")
                            file_strain.append(float(line[idx+19:idx+23]))
                            idx = line.find("-u")
                            file_mu.append(float(line[idx+2:idx+8]))
                            idx = line.find("-T")
                            file_T.append(float(line[idx+2:idx+8]))
                        except:
                            print("There is problem with your 'log' or 'estimator' files")
                            raise Filewrongformat

    file_names = [get_estimator_names(base_dir,ids) for ids in file_ids]
    
    header_lines = 1
    n = []
    nerr = []
    totN = []
    boxdims = []
    kap = []
    kaperr = []
    rhos = []
    rhoserr = []
    estsize = []
    supsize = []
    com_avelst = []
    com_errlst = []
    
    for i, fnames in enumerate(file_names):
        if noenergy:
            estData = np.genfromtxt(fnames['estimator'],\
                                names=True,\
                                skip_header=header_lines,\
                                deletechars="",\
                                usecols = (0,1) ) # Only store N, N^2 column
        else:
            estData = np.genfromtxt(fnames['estimator'],\
                                names=True,\
                                skip_header=header_lines,\
                                deletechars="",\
                                usecols = (9,10) ) # Only store N, N^2 column
        
        numsize = estData.size
        estsize.append(numsize) #save
        
        if isinstance(skip, float) and skip <= 1.0:
            skipfrom = int(numsize*skip)
        elif isinstance(skip, int):
            skipfrom = skip
        else:
            raise Filewrongformat('skip')

        # Simple average and error for N
        ave,err = stats(estData['N'][skipfrom:])   #skip first 'skipfrom' rows
        
        n.append(ave) # save
        nerr.append(err) # save

        # Total number of sites
        # first, read box shape
        with open(fnames['log']) as thefile:
            for line in thefile:
                if line.startswith("Container Dimensions"):
                    boxdim = list(line[len("Container Dimensions    	:	"):-1].split('x'))
                    boxdim = [float(item.strip()) for item in boxdim]
            #
            strain_local = file_strain[i]+1
            totN_local = find_total_Nsites(boxdim[0],boxdim[1],strain_local)

        totN.append(totN_local) # save
        boxdims.append(boxdim)

        #Calculate average compressibility and error using jackknife method
        box_slice_area = int(totN_local)* math.sqrt(3)/4*6  # \AA^2
        temperature = file_T[i] # K
        κ_avg, κ_err = ea.jackknife_on_function( compressibility,
                                                 estData['N'][skipfrom:],
                                                 estData['N^2'][skipfrom:],
                                                 volume=box_slice_area,
                                                 k_B=1.0,
                                                 T=temperature )      #skip first 'skipfrom' rows
        kap.append(κ_avg) # save
        kaperr.append(κ_err) # save

        # Calculate superfluid density fraction
        if 'super' in fnames:
            supData = np.genfromtxt(fnames['super'],\
                                    names=True,\
                                    skip_header=header_lines,\
                                    deletechars="",\
                                    usecols = (0,1,2,3) )  #only store rhos/rho,Wx,Wy,Wz column
            supsize.append(supData.size) #save
            
            if replace:
                supData['rho_s/rho'] = np.nan_to_num(supData['rho_s/rho'])

            s_ave,s_err = stats(supData['rho_s/rho'][skipfrom:])      #skip first 'skipfrom' rows

            # binning analysis
            s_errbin = bin(supData['rho_s/rho'][skipfrom:])        #skip first 'skipfrom' rows

            rhos.append(s_ave) # save
            rhoserr.append(np.amax(s_errbin)) # save
            
        else:
            rhos.append('N/A') # save
            rhoserr.append('N/A') # save
            supsize.append('N/A') # save
            
        # Calculate superfluid density fraction
        if 'scom' in fnames:
            opData = np.genfromtxt(fnames['scom'], skip_header=2,deletechars="")
            com_ave, com_err = stats(opData[skipfrom:])   #skip first 'skipfrom' rows
            com_avelst.append(com_ave)
            com_errlst.append(com_err)
            
        else:
            com_avelst.append('N/A') # save
            com_errlst.append('N/A') # save
            
    dbase = pd.DataFrame({'id': file_ids,
                          'strain': file_strain,
                          'mu': file_mu,
                          'T': file_T,
                          'filedic': file_names,
                          'n': n,
                          'nerr': nerr,
                          'totN': totN,
                          'kap': kap,
                          'kaperr': kaperr,
                          'rhos': rhos,
                          'rhoserr': rhoserr,
                          'estsize': estsize,
                          'supsize': supsize,
                          'boxdims': boxdims,
                          'com_ave': com_avelst,
                          'com_err': com_errlst,
                         })
    
    return dbase

# -----------------------------------------------------------------------------
usage = '''Basic usage : use 'df_pimc(base_dir, skip, noenergy)' to generate dataframe of your pimc results in your code
    
    base_dir: base directory
    skip: if integer, skip the (skip) number of lines. if 0=<skip<1, skip the ratio. Default: 0
    noenergy: if target estimator file has no energy coloumns (K,V,V_ext,V_int,E,E_mu,K/N,V/N,E/N), then set True. Default: 0
    '''
outputs = """output data columns - 
    'id': pimc id
    'strain': iso-strain
    'mu': chemical potential
    'T': temperature
    'filedic': file_names dictionary ex) 'log', 'estimator', 'super'
    'n': number of particle
    'nerr': standard error of n
    'totN': total number of adsorption sites based on graphene
    'kap': compressibility
    'kaperr': standard error of kap
    'rhos': superfulid fraction
    'rhoserr': maximum error of binning error analysis
    'estsize': number of rows in estimator file
    'supsize': number of rows in superfulid file
    """
# -----------------------------------------------------------------------------
def help():
    print(summary)
    print(usage)
    print(outputs)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("This is module")
    print(summary)
    print(usage)
    print(outputs)
