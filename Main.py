import  random 
from Elements import Monomer, Polymer, SimGrid
#from Elements import T
from MCBFM import MonteCarlo
from SimParameters import CoolRate, FileCoolRate, n,N,SaveDataTime,SL,T


MonomerListTemp =[]

polymerTemp =[]
polymers = []

Ht=round(N/10)
Pt=round(N/20)

for i in range (0,n):
    MonomerListTemp = []
    polymerTemp = Polymer(N)
    for j in range (0,Ht):

        for k in range (0,Pt):
            MonomerListTemp.append(Monomer())
            MonomerListTemp[-1].set_monomer(i+1,j*2+1,k*2+1,0)



    polymerTemp.__set_Polymer__(MonomerListTemp)
    polymers.append(polymerTemp)


SimSpace = SimGrid(SL,T)

DataSim=MonteCarlo (polymers,SimSpace)

print("Calculating the energy Level")

for i in range (0,len(DataSim.PolymerPool)):
    
    for j in range (0,len(DataSim.PolymerPool[i].Chain)):
        
        DataSim.PolymerPool[i].Chain[j].E= DataSim.Hamiltonan(i,j,0,0,0)


#DataSim.PrintData("100",100)
print("Relax the created polymers")
for i in range (0,100000000000):
    
    for j in range (0,len(DataSim.PolymerPool)):
        
        PoIn=random.randint(0,len(DataSim.PolymerPool)-1)
        MoIn=random.randint(0,len(DataSim.PolymerPool[PoIn].Chain)-1)
        
        DataSim.MCStep(PoIn,MoIn)
        
DataSim.PrintData(FileCoolRate,0)
print("Start the Simulation")
for i in range (0,100000000000): # Simulation time
    DataSim.TempControll(i)
    if (i%SaveDataTime==0):
        DataSim.PrintData(FileCoolRate,i)
    
    for j in range (0,len(DataSim.PolymerPool)):
        
        PoIn=random.randint(0,len(DataSim.PolymerPool)-1)
        MoIn=random.randint(0,len(DataSim.PolymerPool[PoIn].Chain)-1)
        
        DataSim.MCStep(PoIn,MoIn)


    
        








