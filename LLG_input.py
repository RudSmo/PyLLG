#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 16:09:05 2019

@author: rudolf
"""
import numpy as np
import time
from scipy.constants import physical_constants
from scipy.integrate import solve_ivp
from BoundaryConditions import ob,pb,xpb,ypb
from LLG_src import LLG, ConstConversion#, beffective 
import os
import datetime



global ex
global ey
global ez
ex = np.array([1,0,0])
ey = np.array([0,1,0])
ez = np.array([0,0,1])



# Input parameters

si = np.array([[0,1,0],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],
               [0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1]])#/(np.sqrt(0.1**2+0.9**2)))])       #,1,0.1,0,0.5,0.3,0])

#k = np.random.randn(25,3)
#si = np.array([k[i]/(np.sqrt(np.dot(k[i],k[i]))) for i in range(len(k))])
s = si.flatten()

scd = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),
       np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]


bcon = "(10)"             # "open","periodic"

g = physical_constants["electron gyromag. ratio"][0] 

lx = 10
ly = 10

size = (lx,ly)               # (lx,ly)


if bcon == "(00)":
    inds = ob(int(lx*ly),ly,lx)   
elif bcon == "(11)":
    inds = pb(int(lx*ly),ly,lx)
elif bcon == "(10)":
    inds = xpb(int(lx*ly),ly,lx)
elif bcon == "(01)":    
    inds = ypb(int(lx*ly),ly,lx)



l = 0.1                     # float

j = 1                       # in eV 

jsd = 0                     # in eV
 
hani = 0.00                  # in A/m

hsat = 0.0000                    # in A/m

hdemag = 1*hsat                # in A/m

a1 = 0.00

a2 = 0.000                  

hap = 0                  # in A/m

hdir = ez
hdirect = "ez"
    
happ = hap*hdir #10*hdir

dim = 2

mus = 2 #hsat/(size)

dmi = [0.05,0,0]

T = 0 #in Kelvin

# [j,jsd,scd,k,d,hext,bc,mu]
magparams = [j,jsd,scd,hani,hdemag,a1,a2,happ,bcon,dim,mus,T,dmi,inds]

params = [g,l,size,magparams] 

# =============================================================================
# 
#tim = np.linspace(0,20*1e-10,62)
t1 = time.time()
#print(LLG(0,s,params))
sol= solve_ivp(lambda t,s: LLG(t,s,params),[0,20*1e-12], s, method='RK23')#,first_step=500*1e-15,max_step=1000*1e-13)#,max_step=80*1e-15)
t2 = time.time()-t1
print(t2)  
#     
# # Output
t3 = time.time()
pat = "data/TestBC/noBext/zeros/"+str(lx)+"x"+str(ly)+str(bcon)+"1/"
try:
    os.mkdir(pat)
except:
    pass
magn = open(pat+"magnetization.dat","w")
solutions_x = np.zeros((len(sol.t),int(len(sol.y)/3)))
solutions_y = np.zeros((len(sol.t),int(len(sol.y)/3)))
solutions_z = np.zeros((len(sol.t),int(len(sol.y)/3)))
 
for n in range(len(sol.t)):
    for m in range(int(len(sol.y)/3)):
        solutions_x[n,m] = sol.y[0+3*m][n]
        solutions_y[n,m] = sol.y[1+3*m][n]
        solutions_z[n,m] = sol.y[2+3*m][n]
 
for i in range(len(sol.t)):
    magn.write("%s %s %s %s\n"%(sol.t[i],sum(solutions_x[i,:])/(lx*ly),sum(solutions_y[i,:])/(lx*ly),sum(solutions_z[i,:])/(lx*ly)))
magn.close()
 
 
spin_file = "spin"
for i in range(int(len(sol.y)/3)):
    spi = open(pat+spin_file+"_"+str(i+1)+".dat","w")
    for ji in range(len(sol.y[0])):
        spi.write("%s %s %s %s\n"%(sol.t[ji],sol.y[0+3*i][ji],sol.y[1+3*i][ji],sol.y[2+3*i][ji]))
spi.close()
 
t4 = time.time()-t3
print(t4)
# =============================================================================
# path_to_data = "../../data/"
# parent_folder_name =str(lx)+"x"+str(ly)+"/LLG/"
# try:
#     os.mkdir(path_to_data+parent_folder_name)
# except:
#     pass
# folder_name = "j"+str(j)+"_jsd"+str(jsd)+"_l"+str(l)+"_hani"+str(hani)+"_hsat"+str(hsat)+"_happ"+str(hap)+str(hdirect)+"_"+"bc"+str(bcon)+"/"
# direct = path_to_data+parent_folder_name+folder_name
# try:
#     os.mkdir(direct)
# except:
#     pass
# 
# file_prefix = "LLG_"+str(lx)+"x"+str(ly)+"_l"+str(l)+"_"+"_j"+str(j)+"_jsd"+str(jsd)+"_"+bcon
# 
# # 
# log_file = "logfile.dat"
# log = open(direct+log_file,"w")
# log.write("LOGS LLG SOLVER \n")
# log.write("Date/Time                        %s \n"%(datetime.datetime.now()))
# log.write("Initial spin config s_0          %s \n"%s)
# log.write("Charge density spin pol scd      %s \n"%(scd))
# log.write("Boundary conditions              %s \n"%(bcon))
# log.write("Size (lx,ly)                     %s \n"%(size))
# log.write("gilbert damping l                %s \n"%(l))
# log.write("Exchange interaction j           %s \n"%(j))
# log.write("S-d interaction jsd              %s \n"%(jsd))
# log.write("Anisotropy field hani            %s \n"%(hani))
# log.write("Saturation field hsat            %s \n"%(hsat))
# log.write("Applied field happ               %s \n"%(happ))
# log.write("Runtime                          %s s \n"%(t2))
# log.close()
# # 
# spin_file = "spin"
# for i in range(int(len(sol.y)/3)):
#     spi = open(direct+spin_file+"_"+str(i+1)+".dat","w")
#     for ji in range(len(sol.y[0])):
#         spi.write("%s %s %s %s\n"%(sol.t[ji],sol.y[0+3*i][ji],sol.y[1+3*i][ji],sol.y[2+3*i][ji]))
# spi.close()
# # 
# mag_file = "tot_magnetization.dat"
# solutions_x = np.zeros((len(sol.t),int(len(sol.y)/3)))
# solutions_y = np.zeros((len(sol.t),int(len(sol.y)/3)))
# solutions_z = np.zeros((len(sol.t),int(len(sol.y)/3)))
# #print(sol.y)
# for n in range(len(sol.t)):
#     for m in range(int(len(sol.y)/3)):
#         solutions_x[n,m] = sol.y[0+3*m][n]
#         solutions_y[n,m] = sol.y[1+3*m][n]
#         solutions_z[n,m] = sol.y[2+3*m][n]
# 
# magn = open(direct+mag_file,"w")
# for i in range(len(sol.t)):
#     magn.write("%s %s %s %s\n"%(sol.t[i],sum(solutions_x[i,:])/size,sum(solutions_y[i,:])/size,sum(solutions_z[i,:])/size))
# magn.close()
# #=============================================================================
# =============================================================================


#Time-dep Config files    
#time_folder = "time-dependent/"
#try:
#    os.mkdir(direct+time_folder)
#except:
 #   pass
#for i in range(0,len(sol.t),log_tstep):
#    for j in range(int(len(sol.y)/3)):
#        spin_folder = "spin_"+str(j+1)+"/"
#        try:
#            os.mkdir(direct+time_folder+spin_folder)
#        except:
#            pass        
#        time_file = open(direct+time_folder+spin_folder+"spin_config_t"+str(sol.t[i])+"spin"+str(j+1)+".dat","w")
#        time_file.write("%s %s %s"%(sol.y[0+3*j][i],sol.y[1+3*j][i],sol.y[2+3*j][i]))
#        time_file.close()








 





