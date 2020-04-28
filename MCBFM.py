import numpy as np
import random
import math as MT

import os

from Elements import Monomer, Polymer, SimGrid
from SimParameters import kbo, l0, theta0, kbe, ktor, rho, epsilon, CEbo, CEbe, CEnb, CEtr,CoolRate,kb,ST,T

def randsignal ():
    return 1 if random.random() < 0.5 else -1


def Randomizer ():
    Xc = randsignal()*ST*random.uniform(0, 1)
    Yc = randsignal()*ST*random.uniform(0, 1)
    Zc = randsignal()*ST*random.uniform(0, 1)

    return Xc,Yc,Zc

def distanceV(P1, P2):
    dist2 = MT.pow(P1[0] - P2[0], 2) + MT.pow(P1[1] - P2[1], 2) + MT.pow(P1[2] - P2[2], 2)
    return MT.pow(dist2, 0.5)


def Vvdw(R):
    Enb = 4 * epsilon * (MT.pow(rho / R, 12) - MT.pow(rho / R, 6))
    return Enb


class MonteCarlo:
    def __init__(self, PolymerPool, SimulationSpace):
        self.PolymerPool = PolymerPool
        self.SimulationSpace = SimulationSpace

    def BondEnergy(self, PoIn, MoIn, Xc, Yc, Zc):
        Xn = self.PolymerPool[PoIn].Chain[MoIn].X + Xc

        Yn = self.PolymerPool[PoIn].Chain[MoIn].Y + Yc

        Zn = self.PolymerPool[PoIn].Chain[MoIn].Z + Zc

        P0 = np.asarray([Xn, Yn, Zn])

        if (MoIn > 0 and MoIn < len(self.PolymerPool[PoIn].Chain) - 1):

            Xp = self.PolymerPool[PoIn].Chain[MoIn - 1].X

            Yp = self.PolymerPool[PoIn].Chain[MoIn - 1].Y

            Zp = self.PolymerPool[PoIn].Chain[MoIn - 1].Z

            Xa = self.PolymerPool[PoIn].Chain[MoIn + 1].X

            Ya = self.PolymerPool[PoIn].Chain[MoIn + 1].Y

            Za = self.PolymerPool[PoIn].Chain[MoIn + 1].Z

            P1 = np.asarray([Xp, Yp, Zp])
            D1 = distanceV(P0, P1)
            Ebo1 = kbo * (D1 - l0) * (D1 - l0)

            P2 = np.asarray([Xa, Ya, Za])
            D2 = distanceV(P0, P2)
            Ebo2 = kbo * (D2 - l0) * (D2 - l0)

            Ebo = Ebo2 + Ebo1



        elif (MoIn == 0):

            Xa = self.PolymerPool[PoIn].Chain[MoIn + 1].X

            Ya = self.PolymerPool[PoIn].Chain[MoIn + 1].Y

            Za = self.PolymerPool[PoIn].Chain[MoIn + 1].Z

            P2 = np.asarray([Xa, Ya, Za])

            D2 = distanceV(P0, P2)

            Ebo2 = kbo * (D2 - l0) * (D2 - l0)

            Ebo = Ebo2


        elif (MoIn == len(self.PolymerPool[PoIn].Chain) - 1):

            Xp = self.PolymerPool[PoIn].Chain[MoIn - 1].X

            Yp = self.PolymerPool[PoIn].Chain[MoIn - 1].Y

            Zp = self.PolymerPool[PoIn].Chain[MoIn - 1].Z

            P1 = np.asarray([Xp, Yp, Zp])

            D1 = distanceV(P0, P1)

            Ebo1 = kbo * (D1 - l0) * (D1 - l0)

            Ebo = Ebo1

        return Ebo

    def BendEnergy(self, PoIn, MoIn, Xc, Yc, Zc):

        Xn = self.PolymerPool[PoIn].Chain[MoIn].X + Xc

        Yn = self.PolymerPool[PoIn].Chain[MoIn].Y + Yc

        Zn = self.PolymerPool[PoIn].Chain[MoIn].Z + Zc

        P0 = np.asarray([Xn, Yn, Zn])

        Ebe = 0

        if (MoIn > 0 and MoIn < len(self.PolymerPool[PoIn].Chain) - 1):
            Xp = self.PolymerPool[PoIn].Chain[MoIn - 1].X

            Yp = self.PolymerPool[PoIn].Chain[MoIn - 1].Y

            Zp = self.PolymerPool[PoIn].Chain[MoIn - 1].Z

            Xa = self.PolymerPool[PoIn].Chain[MoIn + 1].X

            Ya = self.PolymerPool[PoIn].Chain[MoIn + 1].Y

            Za = self.PolymerPool[PoIn].Chain[MoIn + 1].Z

            P1 = np.asarray([Xp, Yp, Zp])
            D1 = distanceV(P0, P1)

            P2 = np.asarray([Xa, Ya, Za])
            D2 = distanceV(P0, P2)

            D3 = distanceV(P1, P2)

            CosTheta = (D2 * D2 + D1 * D1 - D3 * D3) / (2 * D2 * D1)

            Ebe = kbe * (CosTheta - MT.cos(theta0)) ** 2

        return Ebe

    def TorsionEnergy(self, PoIn, MoIn, Xc, Yc, Zc):
        Etor=0

        if MoIn < len(self.PolymerPool[PoIn].Chain) - 3:
            P0 = [(self.PolymerPool[PoIn].Chain[MoIn].X + Xc), (self.PolymerPool[PoIn].Chain[MoIn].Y + Yc),
                  (self.PolymerPool[PoIn].Chain[MoIn].Z + Zc)]
            
            P1 = [(self.PolymerPool[PoIn].Chain[MoIn + 1].X), (self.PolymerPool[PoIn].Chain[MoIn + 1].Y),
                  (self.PolymerPool[PoIn].Chain[MoIn + 1].Z)]
            
            P2 = [(self.PolymerPool[PoIn].Chain[MoIn + 2].X), (self.PolymerPool[PoIn].Chain[MoIn + 2].Y),
                  (self.PolymerPool[PoIn].Chain[MoIn + 2].Z)]
            
            P3 = [(self.PolymerPool[PoIn].Chain[MoIn + 3].X), (self.PolymerPool[PoIn].Chain[MoIn + 3].Y),
                  (self.PolymerPool[PoIn].Chain[MoIn + 3].Z)]
            

        else:
            P0 = [(self.PolymerPool[PoIn].Chain[MoIn].X + Xc), (self.PolymerPool[PoIn].Chain[MoIn].Y + Yc),
                  (self.PolymerPool[PoIn].Chain[MoIn].Z + Zc)]
            P1 = [(self.PolymerPool[PoIn].Chain[MoIn - 1].X), (self.PolymerPool[PoIn].Chain[MoIn - 1].Y),
                  (self.PolymerPool[PoIn].Chain[MoIn - 1].Z)]
            P2 = [(self.PolymerPool[PoIn].Chain[MoIn - 2].X), (self.PolymerPool[PoIn].Chain[MoIn - 2].Y),
                  (self.PolymerPool[PoIn].Chain[MoIn - 2].Z)]
            P3 = [(self.PolymerPool[PoIn].Chain[MoIn - 3].X), (self.PolymerPool[PoIn].Chain[MoIn - 3].Y),
                  (self.PolymerPool[PoIn].Chain[MoIn - 3].Z)]
        
        
        q1 = np.subtract(P1, P0)  # b -a
        q2 = np.subtract(P2, P1)  # c -b
        q3 = np.subtract(P3, P2)  # d -c
       
        # Calculate cross vectors
        q1_x_q2 = np.cross(q1, q2)
        
        q2_x_q3 = np.cross(q2, q3)
        
        
        if (not (np.sqrt(np.dot(q1_x_q2, q1_x_q2)) == 0 or np.sqrt(np.dot(q2_x_q3, q2_x_q3)) == 0) ):
            n1 = q1_x_q2 / np.sqrt(np.dot(q1_x_q2, q1_x_q2))
            n2 = q2_x_q3 / np.sqrt(np.dot(q2_x_q3, q2_x_q3))

            # Calculate unit vectors
            u1 = n2
            u3 = q2 / (np.sqrt(np.dot(q2, q2)))
            u2 = np.cross(u3, u1)
    
            # Calculate cosine and sine
            cos_theta = np.dot(n1, u1)
    
            # sin_theta = np.dot(n1, u2)
    
            Etor = ktor * cos_theta
            
        else:
            Etor=0
        
        return Etor

    def LJPotential(self, PoIn, MoIn, Xc, Yc, Zc):

        Enb = 0

        P0 = [(self.PolymerPool[PoIn].Chain[MoIn].X + Xc), (self.PolymerPool[PoIn].Chain[MoIn].Y + Yc),
              (self.PolymerPool[PoIn].Chain[MoIn].Z + Zc)]

        for i in range(0, len(self.PolymerPool)):

            for j in range(0, len(self.PolymerPool[i].Chain)):

                if i == PoIn and j == MoIn:
                    continue
                else:
                    P1 = [self.PolymerPool[PoIn].Chain[MoIn].X, self.PolymerPool[PoIn].Chain[MoIn].Y,
                          (self.PolymerPool[PoIn].Chain[MoIn].Z)]
                    Er = distanceV(P1, P0)

                    if ((Er > 0.1 * rho) and (Er < (rho / 6))):
                        Enb = Enb + Vvdw(Er)
        return Enb

    def Hamiltonan(self, PoIn, MoIn, Xc, Yc, Zc):

        Ebo = 0
        Ebe = 0
        ETr = 0
        ENb = 0
        E = 0

        if (CEbo):
            Ebo = self.BondEnergy( PoIn, MoIn, Xc, Yc, Zc)
            #print(Ebo)

        if (CEbe):
            Ebe = self.BendEnergy( PoIn, MoIn, Xc, Yc, Zc)
            #print(Ebe)

        if (CEtr):
            ETr = self.TorsionEnergy( PoIn, MoIn, Xc, Yc, Zc)
            #print(ETr)

        if (CEnb):
            ENb = self.LJPotential( PoIn, MoIn, Xc, Yc, Zc)
            #print(ENb)
        E = Ebo + Ebe + ETr + ENb

        return E

    def MCStep(self,PoIn, MoIn):

        Xc,Yc,Zc = Randomizer()

        Xt,Yt,Zt =self.PolymerPool[PoIn].Chain[MoIn].X + Xc , self.PolymerPool[PoIn].Chain[MoIn].Y + Yc, self.PolymerPool[PoIn].Chain[MoIn].Z + Zc

        if Xt<0 or Xt> self.SimulationSpace.L or Yt<0 or Yt> self.SimulationSpace.L or 0 > Zt > self.SimulationSpace.L:

            return 0

        else:
            E1=self.PolymerPool[PoIn].Chain[MoIn].E

            E2=self.Hamiltonan(PoIn,MoIn,Xc,Yc,Zc)

            r = random.uniform(0, 1)

            if r < MT.exp(-(-E1+E2)/(kb*self.SimulationSpace.T)):
                self.PolymerPool[PoIn].Chain[MoIn].set_monomer(Xt,Yt,Zt,E2)

                return 1
            
    def PrintData(self,CoolR,TimeStep):
        Tfolder=round(self.SimulationSpace.T)
        Folder=CoolR+"/"+str(Tfolder)
        FinalLoc=CoolR+"/"+str(Tfolder)+"/"+str(TimeStep)
        
        if not os.path.exists(CoolR):
            os.mkdir(CoolR)
        if not os.path.exists(Folder):
            os.mkdir(Folder)
        if not os.path.exists(FinalLoc):
            os.mkdir(FinalLoc)
        
        for i in range (0,len(self.PolymerPool)):
            FileName= FinalLoc+"/"+str(i)+".txt"
            f= open(FileName,"w")
            
            for j in range (0,len(self.PolymerPool[i].Chain)):
                f.write(str(self.PolymerPool[i].Chain[j].X)+" ")
                f.write(str(self.PolymerPool[i].Chain[j].Y)+" ")
                f.write(str(self.PolymerPool[i].Chain[j].Z)+" ")
                f.write(str(self.PolymerPool[i].Chain[j].E)+"\n")
                
                    
            f.close()
            
    def TempControll (self,Time):
        self.SimulationSpace.T = T - (Time/CoolRate)*MT.pow(1,-9) # This time is suitable for Kink and slether and Bond fluctuation
            
            
            







