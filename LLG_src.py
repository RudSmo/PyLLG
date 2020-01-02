#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 16:08:21 2019

@author: rudolf
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
#rc('text', usetex=True)
import time
from scipy.constants import physical_constants
from BoundaryConditions import ob,xpb,ypb,pb

# Auxilliary functions for Heuns predictor-corrector model
global ex
global ey
global ez
ex = np.array([1,0,0])
ey = np.array([0,1,0])
ez = np.array([0,0,1])


def beffective(size,s,magparams):
    """ Auxilliary function to compute the effective magnetic field
        In: size        :: int
            s           :: list of size "size" with elements as numpy nd array
            magparams   :: list of parameters [j,jsd,scd,k,d,hext,bc,mu]
        Out:
            llg         :: list of effective fields at sites i   
    """
    lx,ly = size
    bef = np.zeros((lx*ly,3))
    j,jsd,scd,k,d,a1,a2,hext,bcon,dim,mu,T,dmi,inds = magparams
    if dim == 1:
        if lx == 1 and ly == 1:
            ba = (- jsd * scd[0])/mu - (2*k*s[0][2]*ez)/mu + (2*d*s[0][1]*ey)/mu - hext
            bef[0,0] = ba[0]
            bef[0,1] = ba[1]
            bef[0,2] = ba[2]
        else:
            for i in range(lx):
                 ba = (- jsd * scd[i])/mu - (2*k*s[i][2]*ez)/mu + (2*d*s[i][1]*ey)/mu - hext  
                 if i > 0 and i < (lx-1):
                     bj = -j*(s[i-1]+s[i+1])
                 if i == 0:
                     if bcon == "(1)":
                         bj = -j*(s[i+1]+s[lx-1])
                     elif bcon == "(0)":
                         bj =- j*s[i+1]
                 if i == (lx-1):
                     if bcon == "(1)":
                         bj = -j*(s[i-1]+s[0]) 
                     elif bcon == "(0)":
                         bj = -j*s[i-1]
                 b = -bj/mu - ba
                 bef[i,0] = b[0]
                 bef[i,1] = b[1]
                 bef[i,2] = b[2]
    elif dim == 2:
        l = [(i,np.where(inds[i,:]==True)[0]) for i in range(int(lx*ly))]
        for i in range(lx*ly):
            ba = (-jsd * scd[i])/mu - (2*k*s[i][2]*ez)/mu + (2*d*s[i][1]*ey)/mu - hext + 4*a1*((s[i][0]**3)*ex+(s[i][1]**3)*ey+(s[i][2]**3)*ez)
            snn = [s[ii] for ii in l[i][1]]
            bj = sum(np.array([-j*(snn[ii]) for ii in range(len(snn))]))
            dd = np.array([dmi[0]*ex,dmi[1]*ey,dmi[2]*ez])
            D = np.array([np.cross(np.dot(dd,(s[i]-snn[j])),ez) for j in range(len(snn))]) #/(np.sqrt(np.dot((s[i]-snn[j]),(s[i]-snn[j]))))
            bdmi = sum(np.array([np.cross(D[j],snn[j]) for j in range(len(snn))]))
            bani2 = a2*sum(np.array([s[i][0]*s[ii][0]+s[i][1]*s[ii][1] for ii in l[i][1]]))
            #bth = np.sqrt(2*l/(1+l**2)*kb*T/())*eta
            b = -bj/mu -ba + bdmi/mu + bani2/mu
            bef[i,0] = b[0]
            bef[i,1] = b[1]
            bef[i,2] = b[2]
    else:
        return print("Dimension ",dim," is not implemented! Try dim = 1 or dim = 2 instead.")
    return bef

def LLG(t,s,params):
    g,l,size,magparams = params
    lx,ly = size
    # CONVERT INPUT TO LIST OF INDIVIDUAL SPINS [sk] :: list [s_1,...,s_n] where s_i :: numpy ndarray [s_i^x,s_i^y,s_i^z]
    sk = np.array([[s[0+3*i],s[1+3*i],s[2+3*i]] for i in range(lx*ly)])
    # COMPUTE EFFECTIVE FIELD [beff] :: list [beff_1,...,beff_n] where beff_i :: numpy ndarray [beff_i^x,beff_i^y,beff_i^z]
    beff = beffective(size,sk,magparams)
    # COMPUTE LIST OF RHS OF LLG EQUATION FOR EACH SPIN INDIVIDUALLY <- Vectorize! 
    llg = np.array([-g/(1+l**2)*(np.cross(sk[i],beff[i])+l*np.cross(sk[i],np.cross(sk[i],beff[i]))) for i in range(lx*ly)])
    llgn = llg.flatten()
    return llgn

# CONSTANTS AND INITIALIZATION
def ConstConversion(unit,scd,bc,size,l,j,jsd,hani,hsat,happ):
    if unit == "SI":
        g = physical_constants["electron gyromag. ratio"][0] 
        jn = j*1.60218e-19
        jsdn = jsd*1.60218e-19
        hanin = hani
        hsatn = hsat
        happn = happ
        mus = hsat/(size)
        magparams = [jn,jsdn,scd,hanin,0.0*hsatn,happn,bc,mus]
        params = [g,l,size,magparams]
    
    elif unit == "Hartree":
        g = 1
        jn = j*27.21136
        jsdn = jsd*27.21136
        hanin = hani/(7.99085e-9)   
        hsatn = hsat/(7.99085e-9)
        happn = happ/(7.99085e-9)
        mus = hsat/(size)
        magparams = [jn,jsdn,scd,hanin,0.0*hsatn,happn,bc,mus]
        params = [g,l,size,magparams]
        
    return params

