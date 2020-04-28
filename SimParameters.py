# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 10:32:14 2020

@author: mhdho
"""
n = 200 # number of polymers
N = 200 # Length of polymer chain
SL= 300 # Simulation length

kbo = 1 # bond potential 
l0 = 2 # bond length
kbe=1 # bend Potential
theta0=0 # bend angle
ktor=1 # Torsion potential
rho,epsilon=0.4,1 # LJ potential

CEbo = True #enable bond 
CEbe = True #enable Bend 
CEtr = True #enable Torsion 
CEnb = False #enable nonbond lenord John 
PBC = True # Don't modify
 
ST=1 # Step of jump in each MC step don't change
kb=1.38 #Boltzmann constant don't change
T=1200 # Starting tempreture 
CoolRate=1000000 # cooling rate
FileCoolRate="10000" #should be string only
SaveDataTime=100000000 # time to save data