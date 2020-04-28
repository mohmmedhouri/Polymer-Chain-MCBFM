# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:37:19 2020

@author: mhdho



"""
import matplotlib
import matplotlib.pyplot as plt
import statistics 


import os.path
import math as MT
import numpy as np
X=[]

Y=[]

Z=[]




def RadiousOfGyration(X,Y,Z):
        
        
        SumX = 0
        SumY = 0
        SumZ = 0
        
        N = len(X)

        for i in range (0,len(X)):
            SumX = SumX + X[i]
            SumY = SumY + Y[i]
            SumZ = SumZ + Z[i]
        
        CoM = np.asarray([SumX / N, SumY / N, SumZ / N])
        
        
        SumX = 0
        SumY = 0
        SumZ = 0
        Sum = 0
        for i in range(0, len(X)):
            SumX = SumX + (X[i] - CoM[0]) * (X[i] - CoM[0])
            SumY = SumY + (Y[i] - CoM[1]) * (Y[i] - CoM[1])
            SumZ = SumZ + (Z[i] - CoM[2]) * (Z[i] - CoM[2])

        #   Sum = Sum + MT.pow(self.Chain[i].X - self.CoM[0], 2) + pow(self.Chain[i].Y - self.CoM[1], 2) + pow(
        #     self.Chain[i].Z - self.CoM[2], 2)
        SumX = SumX / len(X)
        SumY = SumY / len(X)
        SumZ = SumZ / len(X)
        Sum = SumX + SumY + SumZ
        RoG = MT.pow(Sum , 0.5)

        return RoG
ROG=[]
for TempR in range (1200,900,-1):
    X=[]
    Y=[]
    Z=[]
    E=[]
    Directory='./Iter1/'+str(TempR)+"/1.txt"
    f=open(Directory, "r")
    
    f1=f.readlines()
    
    for i in f1:
        
        Xt,Yt,Zt,Et = i.split(" ")
        X.append(float(Xt))
        Y.append(float(Yt))
        Z.append(float(Zt))
        E.append(float(Et))
        
    ROG.append(RadiousOfGyration(X,Y,Z))
    plt.plot
    
    #print(RadiousOfGyration(X,Y,Z))
        

fig, ax = plt.subplots()
ax.plot(ROG)

plt.show()
#print(max(ROG))
print(statistics.mean (ROG))
