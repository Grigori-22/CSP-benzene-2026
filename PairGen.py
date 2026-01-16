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

def cos0(ang):
 return np.round(np.cos(ang), 6)

def sin0(ang):
 return np.round(np.sin(ang), 6)

'''
# matrices of symmetry elements of isotropic grid - could be used to reduce number of gridpoint using symmetry
cubicMatrices = np.array([[[ 1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0], [ 1, 0, 0], [0, 0, 1]], [[ 1, 0, 0], [0, 0, 1], [0, 1, 0]], [[0, 1, 0], [0, 0, 1], [ 1, 0, 0]], [[0, 0, 1], [ 1, 0, 0], [0, 1, 0]], [[0, 0, 1], [0, 1, 0], [ 1, 0, 0]],
[[-1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0,-1, 0], [ 1, 0, 0], [0, 0, 1]], [[-1, 0, 0], [0, 0, 1], [0, 1, 0]], [[0,-1, 0], [0, 0, 1], [ 1, 0, 0]], [[0, 0,-1], [ 1, 0, 0], [0, 1, 0]], [[0, 0,-1], [0, 1, 0], [ 1, 0, 0]],
[[ 1, 0, 0], [0,-1, 0], [0, 0, 1]], [[0, 1, 0], [-1, 0, 0], [0, 0, 1]], [[ 1, 0, 0], [0, 0,-1], [0, 1, 0]], [[0, 1, 0], [0, 0,-1], [ 1, 0, 0]], [[0, 0, 1], [-1, 0, 0], [0, 1, 0]], [[0, 0, 1], [0,-1, 0], [ 1, 0, 0]],
[[-1, 0, 0], [0,-1, 0], [0, 0, 1]], [[0,-1, 0], [-1, 0, 0], [0, 0, 1]], [[-1, 0, 0], [0, 0,-1], [0, 1, 0]], [[0,-1, 0], [0, 0,-1], [ 1, 0, 0]], [[0, 0,-1], [-1, 0, 0], [0, 1, 0]], [[0, 0,-1], [0,-1, 0], [ 1, 0, 0]],
[[ 1, 0, 0], [0, 1, 0], [0, 0,-1]], [[0, 1, 0], [ 1, 0, 0], [0, 0,-1]], [[ 1, 0, 0], [0, 0, 1], [0,-1, 0]], [[0, 1, 0], [0, 0, 1], [-1, 0, 0]], [[0, 0, 1], [ 1, 0, 0], [0,-1, 0]], [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
[[-1, 0, 0], [0, 1, 0], [0, 0,-1]], [[0,-1, 0], [ 1, 0, 0], [0, 0,-1]], [[-1, 0, 0], [0, 0, 1], [0,-1, 0]], [[0,-1, 0], [0, 0, 1], [-1, 0, 0]], [[0, 0,-1], [ 1, 0, 0], [0,-1, 0]], [[0, 0,-1], [0, 1, 0], [-1, 0, 0]],
[[ 1, 0, 0], [0,-1, 0], [0, 0,-1]], [[0, 1, 0], [-1, 0, 0], [0, 0,-1]], [[ 1, 0, 0], [0, 0,-1], [0,-1, 0]], [[0, 1, 0], [0, 0,-1], [-1, 0, 0]], [[0, 0, 1], [-1, 0, 0], [0,-1, 0]], [[0, 0, 1], [0,-1, 0], [-1, 0, 0]],
[[-1, 0, 0], [0,-1, 0], [0, 0,-1]], [[0,-1, 0], [-1, 0, 0], [0, 0,-1]], [[-1, 0, 0], [0, 0,-1], [0,-1, 0]], [[0,-1, 0], [0, 0,-1], [-1, 0, 0]], [[0, 0,-1], [-1, 0, 0], [0,-1, 0]], [[0, 0,-1], [0,-1, 0], [-1, 0, 0]]])
'''

def GenRotations(angleVectorChange):   #generated rotation matrices for every possible rotation given by angleVectorChange
 rotMaxList = []
 allMatCount = 0
 redMatCount = 0
 angle1 = angleVectorChange[0]*np.pi/180
 angle2 = angleVectorChange[1]*np.pi/180
 angle3 = angleVectorChange[2]*np.pi/180
 angles1 = np.arange(0,2*np.pi,angle1)
 angles2 = np.arange(0,2*np.pi,angle2)
 angles3 = np.arange(0,2*np.pi,angle3)
 for thetaX in angles1:
  for thetaY in angles2:
   for thetaZ in angles3:
    Rx = np.array([[1, 0, 0], [0, cos0(thetaX), -1*sin0(thetaX)], [0, sin0(thetaX), cos0(thetaX)]])
    Ry = np.array([[cos0(thetaY), 0, sin0(thetaY)], [0, 1, 0], [-1*sin0(thetaY), 0, cos0(thetaY)]])
    Rz = np.array([[cos0(thetaZ), -1*sin0(thetaZ), 0], [sin0(thetaZ), cos0(thetaZ), 0], [0, 0, 1]])
    rotMatrix = Rx@Ry@Rz   #matrix is created for every possible combination of angles
    allMatCount += 1
    isInRotMaxList = []
    for matrix in rotMaxList:   #removing duplicate matrices
     isInRotMaxList.append((rotMatrix == matrix).all())
    if True not in isInRotMaxList:
     rotMaxList.append(rotMatrix)
     redMatCount += 1
 print("Generated "+str(redMatCount)+"/"+str(allMatCount)+" rotation matrices. The rest were redundant.")
 return rotMaxList

def MySort(sl, al):   #sorting function that sorts atoms by symbol, then X, then Y, then Z and round all coords to 6 decimal digits
 tempArray = []
 for i, atom in enumerate(al):
  tempArray.append([sl[i], atom[0], atom[1], atom[2]])
 tempArray.sort(key=lambda seq: (seq[0], seq[1], seq[2], seq[3]))
 sortedArray = np.array([np.array([np.round(float(x[1]),6), np.round(float(x[2]),6), np.round(float(x[3]),6)]) for x in tempArray])
 return sortedArray

def RotateMol(sl, al, rotMaxList):   #generating all possible rotations of molecule using rotation matrices and removing duplicates
 al_all = []
 al = [np.array(x) for x in al]
 for rotMatrix in rotMaxList:
  al_rotated = np.transpose(rotMatrix@np.transpose(al))
  al_rotated = np.array([np.array([np.round(float(x[0]),6), np.round(float(x[1]),6), np.round(float(x[2]),6)]) for x in al_rotated])
  al_all.append(al_rotated)
 al = MySort(sl, al)
 toCompare = [al]
 al_red = [al]
 for al_rot in al_all:
  al_sort = MySort(sl, al_rot)
  areTheSame = []
  for each in toCompare:
   isTheSame = (each == al_sort).all()
   areTheSame.append(isTheSame)
  if True not in areTheSame:
   toCompare.append(al_sort)
   al_red.append(al_rot)
 print(str(len(al_red))+" out of "+str(len(al_all))+" structures are unique due to molecule's symmetry.")
 return al_red

def TranslateMol(al, gridpoint, step):   #translating molecule by a gridpoint*stepsize
 gridX = gridpoint[0]*step
 gridY = gridpoint[1]*step
 gridZ = gridpoint[2]*step
 al_translated = []
 for each in al:
  al_translated.append([each[0]+gridX, each[1]+gridY, each[2]+gridZ])
 return al_translated

def GetFurthest(sl, al):   #calculates what is the furthest ever necessary gridpoint using most distant atom and it's vdW radius
 dist = 0
 for i, each in enumerate(al):
  d = np.sqrt(each[0]**2+each[1]**2+each[2]**2)+getVdw(sl[i])
  if d > dist:
   dist = d
 return dist

def MakeGrid(sl, al, step):   #generating grid around a molecule using stepsize and GetFurthest
 maxDist = GetFurthest(sl, al)
 maxPoint = int(np.ceil(2*maxDist/step))
 pointList = [x for x in range(-maxPoint, maxPoint+1)]
 gridList = []
 for each1 in pointList:
  for each2 in pointList:
   for each3 in pointList:
    gridList.append([each1,each2,each3])
 return gridList

def GenValidPoint(sl, al, al_rot, gridpoint, step):   #generating valid molecule positions for rotated molecule
 al_translated = TranslateMol(al_rot, gridpoint, step)
 overlap = CheckIfOverlap(sl, al, al_translated)
 if overlap == True:
  return "Not valid"   #returns appropiate string if rotated molecule in a given gridpoint will overlap with a molecule at 0,0,0
 else:
  interact = CheckIfInteract(sl, al, al_translated)
  if interact == False:
   return "Not valid"   #returns appropiate string if rotated molecule in a given gridpoint will not interact with a molecule at 0,0,0
  else:
   return al_translated   #returns coords if molecule interacts but not overlaps with molecule at 0,0,0

def WriteAFile(fN, noA, sL, aL1, aL2, gridpoint, angle):   #fileName, noOfAtoms, symbolList, atomList 1 and 2, which gridpoint and which angle (for name)
 savingFile = open(fN,"w")                                 #this function generates xyz file of a pair for a valid pairing
 firstLine = str(noA)+"\n"
 secondLine = "Gridpoint "+str(gridpoint)+"; position No. "+str(angle)+" \n"
 savingFile.write(firstLine)
 savingFile.write(secondLine)
 for i, atom in enumerate(aL1):
  thisLine = sL[i]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
  savingFile.write(thisLine)
 for j, atom in enumerate(aL2):
  thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
  savingFile.write(thisLine)
 savingFile.close()

def WriteRotations(sL, al_red):   #function that creates a file containing all rotations of a molecule (not necessary but a good practice to check if they look right)
 fN = "All_Rotations.xyz"
 savingFile = open(fN,"w")
 noA = len(symbolList)
 for i, al_rot in enumerate(al_red):
  firstLine = str(noA)+"\n"
  secondLine = "Position No. "+str(i)+" \n"
  savingFile.write(firstLine)
  savingFile.write(secondLine)
  for j, atom in enumerate(al_rot):
   thisLine = sL[j]+"   "+'{:.6f}'.format(atom[0])+"   "+'{:.6f}'.format(atom[1])+"   "+'{:.6f}'.format(atom[2])+"\n"
   savingFile.write(thisLine)
 savingFile.close()

symbolList, atomList = ReadXYZ(fileName)               #reading benzene
gridList = MakeGrid(symbolList, atomList, step)        #creating grid aroung it
rotMaxList = GenRotations(angleVectorChange)           #generating rotation matrices
al_red = RotateMol(symbolList, atomList, rotMaxList)   #generating all rotated molecules
noOfStruc = 0
WriteRotations(symbolList, al_red)                     #making an xyz file with all rotations
#xx = [[],[],[],[],[],[]]   #these are for visualization purposes
#yy = [[],[],[],[],[],[]]
#zz = [[],[],[],[],[],[]]
for gridPoint in gridList:                                               #for each gridpoint
 for i, al_rot in enumerate(al_red):                                     #for each rotation
  point = GenValidPoint(symbolList, atomList, al_rot, gridPoint, step)   #checks if it's a valid position for a second molecule
  if point == "Not valid":
   continue
  else:
   if gridPoint[0] >= 0 and gridPoint[1] >=0 and gridPoint[2] >= 0:      #for symmetry reasons we generate only 1/8th of a grid with only posistive x, y and z
    noOfStruc += 1
    fN = str(gridPoint[0])+"-"+str(gridPoint[1])+"-"+str(gridPoint[2])+"--"+str(i)+".xyz"
    noA = 2*len(symbolList)
    WriteAFile(fN, noA, symbolList, atomList, point, gridPoint, i)       #writing all pairs into files
#    xx[i].append(gridPoint[0])   #these are for visualization purposes
#    yy[i].append(gridPoint[1])
#    zz[i].append(gridPoint[2])
print("Number of pairs is", noOfStruc)


'''
#This is for creating a picture of gridpoints and rotations that are valid in them
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure()
gs = GridSpec(12, 6, figure=fig)
colours = ['red','blue','green','orange','purple','grey']
ax1 = fig.add_subplot(gs[:6, :2], projection='3d')
ax2 = fig.add_subplot(gs[6:, :2], projection='3d')
ax3 = fig.add_subplot(gs[:6, 2:4], projection='3d')
ax4 = fig.add_subplot(gs[6:, 2:4], projection='3d')
ax5 = fig.add_subplot(gs[:6, 4:], projection='3d')
ax6 = fig.add_subplot(gs[6:, 4:], projection='3d')
weirdSol = [ax1, ax2, ax3, ax4, ax5, ax6]
for i in range(6):
 ax = weirdSol[i]
 ax.scatter(xx[i], yy[i], zz[i], marker='.', c=colours[i])
 ax.set_xticks([])
 ax.set_yticks([])
 ax.set_zticks([])
plt.show()
'''