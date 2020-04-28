import numpy as np
import math as MT

BL = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1],
      [1, 1, 1]]  # represents all the blocked Lattice location

class Atom :

    def __init__(self,Type):
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.Charge=0
        self.E = 0
        self.Type = Type
        if (Type == "C"):
            self.Charge = 0
        elif (Type == "O" ):
            self.Charge = -1
        else :
            self.Charge = 1

    def set_Atom (self,Xs,Ys,Zs,Es,Ts):
        self.X = Xs
        self.Y = Ys
        self.Z = Zs
        self.E = Es
        self.Type = Ts
        if (Ts == "C"):
            self.Charge = 0
        elif (Ts == "O" ):
            self.Charge = -1
        else :
            self.Charge = 1

    def updateCoordiate(self, Xs,Ys,Zs,Es):
        self.X = Xs
        self.Y = Ys
        self.Z = Zs
        self.E = Es




class Monomer:

    def __init__(self):
        
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.E = 0
        self.GOA=[Atom("C"),Atom("O"),Atom("H")]
        self.GOA[0].updateCoordiate(self.X,self.Y,self.Z)
        self.GOA[1].updateCoordiate(self.X + 1, self.Y + 1, self.Z + 1)
        self.GOA[2].updateCoordiate(self.X - 1, self.Y - 1, self.Z - 1)


    def set_monomer(self, Xs, Ys, Zs,Es,GOAs):
        # self.MID = MIDs
        self.X = Xs
        self.Y = Ys
        self.Z = Zs
        self.E = Es
        self.GOA = GOAs

    def CenterLocation (self):
        SumX = 0
        SumY = 0
        SumZ = 0
        N = len(self.GOA)

        for i in self.GOA:
            SumX = SumX + i.X
            SumY = SumY + i.Y
            SumZ = SumZ + i.Z

        self.X=  SumX / N
        self.Y = SumY / N
        self.Z = SumZ / N



class Polymer:

    def __init__(self, N):
        self.Chain = []
        for i in range(0, N):
            m_temp = Monomer()
            self.Chain.append(m_temp)

        self.CoM = np.asarray([0, 0, 0])  # Center of Mass

        self.RoG = 0  # Radius of gyration

    def __set_Polymer__(self, Polytemp):
        self.Chain = Polytemp

    def UpdatePolymer(self, Chainlocation):
        """
        This function is used to update the coordinates of the polymer

        :rtype: none
        """
        self.Chain = Chainlocation

        self.CenterOfMass()
        self.RadiousOfGyration()





    def UpdateLocation(self, MonomerIndex, X, Y, Z):
        self.Chain[MonomerIndex].SetCoordinate(X, Y, Z)


    def CenterOfMass(self):  #
        """
        This function is used to calculate the center of mass
        :rtype: none
        """
        SumX = 0
        SumY = 0
        SumZ = 0
        N = len(self.Chain)

        for i in self.Chain:
            i.CenterLocation()
            SumX = SumX + i.X
            SumY = SumY + i.Y
            SumZ = SumZ + i.Z

        self.CoM = np.asarray([SumX / N, SumY / N, SumZ / N])

        return self.CoM

    def RadiousOfGyration(self):
        SumX = 0
        SumY = 0
        SumZ = 0
        Sum = 0

        for i in range(0, len(self.Chain)):
            SumX = SumX + (self.Chain[i].X - self.CoM[0]) * (self.Chain[i].X - self.CoM[0])
            SumY = SumY + (self.Chain[i].Y - self.CoM[1]) * (self.Chain[i].Y - self.CoM[1])
            SumZ = SumZ + (self.Chain[i].Z - self.CoM[1]) * (self.Chain[i].Z - self.CoM[2])

        #   Sum = Sum + MT.pow(self.Chain[i].X - self.CoM[0], 2) + pow(self.Chain[i].Y - self.CoM[1], 2) + pow(
        #     self.Chain[i].Z - self.CoM[2], 2)
        SumX = SumX / len(self.Chain)
        SumY = SumY / len(self.Chain)
        SumZ = SumZ / len(self.Chain)
        sum = SumX + SumY + SumZ
        self.RoG = MT.pow(Sum / len(self.Chain), 0.5)

        return self.RoG

class SimGrid : # this is the simulation space


    def __init__(self, L, T):
        self.L=L
        self.T=T
    def set_T (self,T0 ):
        self.T=T0
