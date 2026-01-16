import numpy as np
#from James_v37 import getAtomicCovalentRadius as getAcr

step = 1   #distance between gridpoints in angstrem
angleVectorChange = [90, 90, 90]   #changes of angles around X, Y and Z axis - in degrees
fileName = "benzene.xyz"   #file with a molecule - here benzene

def getAcr(atomSymb):   #function that returns atomic covalent radius - shortened version of what can be found in James_v37
 if atomSymb == 'H':
  number = 0.31
 if atomSymb == 'C':
  number = 0.76
 return 1.16*number

def getVdw(atomSymb):   #function that returns van der Waals radius - shortened version of what can be found in James_v37
 if atomSymb == 'H':
  number = 1.1
 if atomSymb == 'C':
  number = 1.7
 return 1.2*number

def ReadXYZ(xyzFile):   #function reading an xyz file of molecule
 f = open(xyzFile, "r")
 allLines = f.readlines()
 allLines = [x.strip().split() for x in allLines]
 coords = np.array([np.array([float(x[1]), float(x[2]), float(x[3])]) for x in allLines[2:]])   #saving coords as numpy 2D array after deleting first 2 lines (number of atoms and comment line)
 symbolList = [x[0] for x in allLines[2:]]   #saving symbols of each atom - first column in xyz file
 return symbolList, coords

def CheckIfOverlap(sl, al1, al2):   #function that checks if two molecules overlap - takes symbol list, and both coord lists
 for i,a1 in enumerate(al1):
  for j,a2 in enumerate(al2):
   radSum = getAcr(sl[i])+getAcr(sl[j])
   v = np.subtract(a2,a1)
   dist = np.sqrt(np.dot(v,v))
   if dist < radSum:
    return True   #if at least one atom of first molecule is within covalent radius of at least one atom of second molecule: molecules overlap
 return False

def CheckIfInteract(sl, al1, al2):   #function similar to CheckIfOverlap but using van der waals radii
 for i,a1 in enumerate(al1):
  for j,a2 in enumerate(al2):
   radSum = getVdw(sl[i])+getVdw(sl[j])
   v = np.subtract(a2,a1)
   dist = np.sqrt(np.dot(v,v))
   if dist < radSum:
    return True
 return False

def TranslateMol(al, gridpoint, step):   #translating molecule by a gridpoint*stepsize
 gridX = gridpoint[0]*step
 gridY = gridpoint[1]*step
 gridZ = gridpoint[2]*step
 al_translated = []
 for each in al:
  al_translated.append([each[0]+gridX, each[1]+gridY, each[2]+gridZ])
 return al_translated

def CutInHalf(sl, al):   #function that cuts in half the pairs and only returns the second half (omits the original molecule at 0,0,0)
 sl = sl[12:]
 al = al[12:]
 return sl, al

def Check2Mols(sl, al1, al2):   #function that checks if two molecules are overlaping, interacting or non-interacting - used to calculate energies
 overlap = CheckIfOverlap(sl, al1, al2)
 if overlap == True:
  return "Not valid"
 else:
  interact = CheckIfInteract(sl, al1, al2)
  if interact == False:
   return "Non-interacting"
  else:
   return "Interacting"

def ReadResults(resultsFile):   #function that reads result file containing pairs names and energies of interaction (defined as Epair - E1 - E2 and given in kcal/mol)
 f = open(resultsFile, "r")
 allLines = f.readlines()
 allLines = [x.strip().split() for x in allLines]
 nameList = [x[0] for x in allLines]
 energyList = [x[1] for x in allLines]
 al_all = []   #this will contain all coordinates of second molecule in the same order as names
 for name in nameList:
  sl, al = ReadXYZ(name+'.xyz')
  sl, al = CutInHalf(sl, al)
  al_all.append(al)
 sl, al = ReadXYZ(nameList[0]+'.xyz')
 sl = sl[:12]
 al_original = al[:12]
 return sl, nameList, energyList, al_all, al_original   #function returns names, energies and coords as well as list of symbols and original coords of a molecule in 0,0,0


#list contating which rotation results in which position i.e. after rotation 4 and then rotation 1, resulted rotation is combRot[4][1] so rotation 5 - to be automated in the future
CombRot = [[0,1,2,3,4,5], [1,0,5,4,3,2], [2,3,0,1,5,4], [3,2,4,5,1,0], [4,5,3,2,0,1], [5,4,1,0,2,3]]

#list containing ordinals of those pairs that have rotation 0 to simplify the problem for now
onlyFlat = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,36,39,40,41,42,43,46,47,48,49,52,53,62,63,64,65,70,71,72,73,74,75,76,77,104,105,110,111,112,113,114,115,118,119,120,121,122,123,130,131,132,133,134,135,136,141,152,153,164,165,168,173,178,179,181,182,186,193,204,205,212,222,227,228,229,230,241,242,245,246,247,252,253,269,274,275,276,277,290,311,318,319,323,330,331,334,339,340,341,342,343,344,347,350,351,352,353,354,363,364,365,366,367,371,376,377,378,389,390,394,399,400,403,404,405,406]


def GetPossibleCrystal(howMany, sl, nameList, al_all):   #function that returns a recipe for a crystal given by 3 pair numbers and sum of energies of their interactions
 noOfPairs = howMany
 energies = []
 recipe = []
 for i in range(noOfPairs-2):
  if i in onlyFlat:
   for j in range(i+1, noOfPairs-1):
    if j in onlyFlat:
     for k in range(j+1, noOfPairs):
      if k in onlyFlat:
       infoMol = [Check2Mols(sl, al_all[i], al_all[j]), Check2Mols(sl, al_all[i], al_all[k]), Check2Mols(sl, al_all[j], al_all[k])]   #information if molecules are valid
#     if "Not valid" not in infoMol and 'Non-interacting' not in infoMol:
       if "Not valid" not in infoMol:   #if all pairs are valid with each other it is a valid crystal structure
        x = nameList[i].split('-')
        y = nameList[j].split('-')
        z = nameList[k].split('-')
        if int(x[4])==0 and int(y[4])==0 and int(z[4])==0:   #simplifying problem for now to only rotation 1, so the number of structures is in 1000s not 1000000s
#        if (int(x[4])==0 or int(x[4])==4) and (int(y[4])==0 or int(y[4])==4) and (int(z[4])==0 or int(z[4])==4):   #only rotation 0 and 4
         energy1 = float(energyList[i])   #getting energies
         energy2 = float(energyList[j])
         energy3 = float(energyList[k])
         name1 = [int(x[0]), int(x[1]), int(x[2]), int(x[4])]   #getting names
         name2 = [int(y[0]), int(y[1]), int(y[2]), int(y[4])]
         name3 = [int(z[0]), int(z[1]), int(z[2]), int(z[4])]      
         rotNew1 = CombRot[name1[3]][name2[3]]    #calculating relative rotation
         rotNew2 = CombRot[name1[3]][name3[3]]
         rotNew3 = CombRot[name2[3]][name3[3]]
         nameNew1 = str(abs(name2[0]-name1[0]))+'-'+str(abs(name2[1]-name1[1]))+'-'+str(abs(name2[2]-name1[2]))+'--'+str(rotNew1)   #figuring out how each of 2 pairs see each other
         nameNew2 = str(abs(name3[0]-name1[0]))+'-'+str(abs(name3[1]-name1[1]))+'-'+str(abs(name3[2]-name1[2]))+'--'+str(rotNew2)
         nameNew3 = str(abs(name3[0]-name2[0]))+'-'+str(abs(name3[1]-name2[1]))+'-'+str(abs(name3[2]-name2[2]))+'--'+str(rotNew3)
         energy4 = 0   #setting energy of pairs as 0, to be overwritten if they interact
         energy5 = 0
         energy6 = 0
         for m, name in enumerate(nameList):
          if nameNew1 == name:
            energy4 = float(energyList[m])
          if nameNew2 == name:
            energy5 = float(energyList[m])
          if nameNew3 == name:
            energy6 = float(energyList[m])
         energies.append([energy1, energy2, energy3, energy4, energy5, energy6])   #energies 1-3 are from molecule in the center interacting with second molecule; energies 4-6 are second molecules interacting with each other
         recipe.append([i, j, k])
 return recipe, energies

def OneOr2(letter):   #necessary to generate .gen file for DFTB+ calculation
 if letter == 'C':
  return '1'
 if letter == 'H':
  return '2'

def WriteCrystal(sL, recipe, energies, al_all, al_original, noA):   #takes recipe and creates a .gen file containing crystal parameters and coords
 counterOfAllCrystals = 0
 for i, energy in enumerate(energies):
  if sum(energy) < 0:   #considering only those crystals where interactions of 1-3 stabilize it more than 4-6 potentially destabilize it
   counterOfAllCrystals += 1
   [x, y, z] = recipe[i]
   fN = "c"+str(x)+"-"+str(y)+"-"+str(z)+".gen"
   savingFile = open(fN,"w")
   firstLine = str(noA)+"    S\n"
   secondLine = "C   H\n"
   savingFile.write(firstLine)   #writing the beginning of .gen file and then coords of all pairs in recipe
   savingFile.write(secondLine)
   for j, atom in enumerate(al_original):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+1)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   for j, atom in enumerate(al_all[x]):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+13)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   for j, atom in enumerate(al_all[y]):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+25)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   for j, atom in enumerate(al_all[z]):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+37)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   nameX = nameList[x].split('-')
   nameY = nameList[y].split('-')
   nameZ = nameList[z].split('-')
   nameX = [int(nameX[0]), int(nameX[1]), int(nameX[2]), int(nameX[4])]
   nameY = [int(nameY[0]), int(nameY[1]), int(nameY[2]), int(nameY[4])]
   nameZ = [int(nameZ[0]), int(nameZ[1]), int(nameZ[2]), int(nameZ[4])]
   al1 = TranslateMol(al_all[x], nameY[:3], step)   #calculating 4 more molecules that will be created from translational combinations of pairs
   al2 = TranslateMol(al_all[x], nameZ[:3], step)
   al3 = TranslateMol(al_all[y], nameZ[:3], step)
   al4 = TranslateMol(al1, nameZ[:3], step)
   vectorX = nameX[:3]   #writing half of the vector in each direction
   vectorY = nameY[:3]
   vectorZ = nameZ[:3]
   for j, atom in enumerate(al1):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+49)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   for j, atom in enumerate(al2):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+61)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   for j, atom in enumerate(al3):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+73)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   for j, atom in enumerate(al4):
#    thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"   #for xyz
    thisLine = "   "+str(j+85)+"   "+OneOr2(sL[j])+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
    savingFile.write(thisLine)
   thisLine = "   "+'{:.6f}'.format(0)+"   "+'{:.6f}'.format(0)+"   "+'{:.6f}'.format(0)+"\n"   # 0 0 0 point required for .gen file
   savingFile.write(thisLine)
   thisLine = "   "+'{:.6f}'.format(2*step*vectorX[0])+"   "+'{:.6f}'.format(2*step*vectorX[1])+"   "+'{:.6f}'.format(2*step*vectorX[2])+"\n"   #vectors of crystal lattice
   savingFile.write(thisLine)
   thisLine = "   "+'{:.6f}'.format(2*step*vectorY[0])+"   "+'{:.6f}'.format(2*step*vectorY[1])+"   "+'{:.6f}'.format(2*step*vectorY[2])+"\n"
   savingFile.write(thisLine)
   thisLine = "   "+'{:.6f}'.format(2*step*vectorZ[0])+"   "+'{:.6f}'.format(2*step*vectorZ[1])+"   "+'{:.6f}'.format(2*step*vectorZ[2])+"\n"
   savingFile.write(thisLine)
   savingFile.close()
 print(counterOfAllCrystals)




sl, nameList, energyList, al_all, al_original = ReadResults('results.dat')   #reading results and preparing coords for crystal generation
recipe, energies = GetPossibleCrystal(407, sl, nameList, al_all)   #obtaining valid recipies for crystals
WriteCrystal(sl, recipe, energies, al_all, al_original, 96)   #creating and writing crystal files
