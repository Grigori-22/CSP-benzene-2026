def readposcar(poscar):

       import sys
       #import numpy
       #from numpy import *
       
       poscarfile = open(poscar)
       allLines=poscarfile.readlines()
       poscarfile.close()
       numIons=0
       cell = []
       coords = []
       lc=float(allLines[1].strip())
       atomtypes=[]       

       atomlabels=allLines[0].split()
       j=0
       v5=0

       for i in range(len(allLines)):
              if j>1 and j<5:
                     line = allLines[i]
                     fline = line.split()
                     fline = [float(fline[0])*lc, float(fline[1])*lc, float(fline[2])*lc]
                     cell.append(fline)

              if j==5:
                      c=0
                      try:
                           test=int(allLines[i].split()[0])
                      except:
                           atomlabels=allLines[i].split()
                           j=j-1
                           v5=1
                           c=1

                      if c==0:       
                           line = allLines[i]
                           numIonsList=line.strip().split()
                           numspecies=len(numIonsList)
                           for z in range(numspecies):
                               numIons=numIons+float(numIonsList[z])
                               for c in range(int(numIonsList[z])):
                                   atomtypes.append(atomlabels[z])

                            
              if j>5: 

                 if allLines[i].strip()=='' or c>numIons:
                      break
                 else:
                      try:
                           test=float(allLines[i].split()[0])
                      except:
                           c=0
                      else:
                            c+=1
                            line = allLines[i]
                            fline = line.strip().split()
                            fline = [float(fline[0]), float(fline[1]), float(fline[2])]
                            coords.append(fline)

                
              j=j+1 

       #coords2 = array(coords)
       #cell2 = array(cell)
       return [coords, cell, atomtypes]

def pytester():
       print('hello')

def dotproduct(v1, v2):
  import math

  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  import math
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  import math
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def add(v1, v2):
  import math
  return [sum(x) for x in zip(v1,v2)]

def subtract(v1, v2):
  import math
  return [a-b for (a,b) in zip(v1, v2)]


def addArrays(a, v):
  import math
  import numpy

  x=numpy.array(a)+numpy.array(v)
  return arrayToList(x)

def subtractArrays(a, v):
  import math
  import numpy

  x=numpy.array(a)-numpy.array(v)
  return arrayToList(x)

def matrixmultiply(a,b):
       c=a@b
       return c


def compareMolecules(atype,btype):
       return sorted(atype)==sorted(btype) 

def sortAtomtypes(coords,atype):
       newcoords=[]
       newatypes=[]

       sortcoords=[]
       sorta=[]

       i=0
       for i1 in atype:
          toadd=1
          for i2 in range(len(sorta)):
               if i1==sorta[i2]:
                    sortcoords[i2].append(coords[i])
                    toadd=0
          if toadd==1:
               sorta.append(i1) 
               sortcoords.append([])
               a1=len(sortcoords)
               sortcoords[a1-1].append(coords[i])
          i=i+1

       b=0
       for b1 in sorta:
           for b2 in sortcoords[b]:
                newcoords.append(b2)
                newatypes.append(b1)
           b=b+1
                
       return [newcoords, newatypes]

def groupMoleculeTypes(coords, atypes, ncomp):
       # the new set of coordinates and atomtypes
       newcoords=[]
       newatypes=[]
       newncomp=[]

       # a list of the distinct molecular compositions (mlist), and the full sorted list of atomic coordinates from each molecular type
       mlistcoords=[]
       matypes=[]
       mlist=[]

       # loop over the elements of ncomp, and use it to fill up mlist and mlistcoords
       # counters - i is iterating ncomp, and j is iterating the atomic coordinates
       i=0
       j=0
       at=makeCopyList1D(atypes)
       for i1 in ncomp:

         currenttype=atypes[i:i+i1]
         currentcoord=coords[i:i+i1]
         print('currently considered mol: ',currenttype)

         toadd=1
         for i2 in range(len(mlist)):
               if compareMolecules(mlist[i2],currenttype):
                    for j1 in range(len(currenttype)):
                        mlistcoords[i2].append(currentcoord[j1])
                        matypes[i2].append(currenttype[j1])
                    toadd=0
         if toadd==1:
               mlist.append(currenttype) 
               mlistcoords.append([])
               matypes.append([])
               a1=len(mlistcoords)
               for j1 in range(len(currenttype)):
                    mlistcoords[a1-1].append(currentcoord[j1])
                    matypes[a1-1].append(currenttype[j1])
         i=i+i1
                        
       # make newncomp
       print('unique molecular compositions :',mlist)
       x=0
       for x1 in mlist:            
            tcoords,tatypes=sortAtomtypes(mlistcoords[x],matypes[x])
            for x2 in tcoords:
                  newcoords.append(x2)
            for x3 in tatypes:
                  newatypes.append(x3)
            newncomp.append(len(tcoords))
            x=x+1
       
       return newcoords, newatypes, newncomp          
              

def getCenter(atomset):

       xcenter=0.0
       ycenter=0.0
       zcenter=0.0

       numatoms=len(atomset)

       for i,bi in enumerate(atomset):
              xcenter=xcenter+bi[0]
              ycenter=ycenter+bi[1]
              zcenter=zcenter+bi[2]

       xcenter=xcenter/numatoms
       ycenter=ycenter/numatoms
       zcenter=zcenter/numatoms

       return [xcenter, ycenter, zcenter]

def getVDWRadius(atom):

# the Average radii from the Wikipedia entry for "Van der Waals radius"

       if atom == 'H':
              number = 1.10
       elif atom == 'Li':
              number = 1.82
       elif atom == 'Be':
              number = 1.53
       elif atom == 'B':
              number = 1.92
       elif atom == 'C':
              number = 1.70
       elif atom == 'N':
              number = 1.55
       elif atom == 'O':
              number = 1.52
       elif atom == 'F':
              number = 1.47
       elif atom == 'Na':
              number = 2.27
       elif atom == 'Mg':
              number = 1.73
       elif atom == 'Al':
              number = 1.84
       elif atom == 'Si':
              number = 2.1
       elif atom == 'P':
              number = 1.8
       elif atom == 'S':
              number = 1.8
       elif atom == 'Cl':
              #number = 1.75
              number = 1.95
       elif atom == 'Ar':
              number = 1.88
       elif atom == 'K':
              number = 2.75
       elif atom == 'Ca':
              number = 2.31
       elif atom == 'Ti':
              number = 1.7
       elif atom == 'Fe':
              number = 1.94
       elif atom == 'Cu':
              number = 1.40
       elif atom == 'Zn':
              number = 1.39
       elif atom == 'Ga':
              number = 1.87
       elif atom == 'Ge':
              number = 2.11
       elif atom == 'As':
              number = 1.85
       elif atom == 'Rb':
              number = 3.0
       elif atom == 'Rh':
              number = 1.72
       elif atom == 'Ag':
              number = 1.72
       elif atom == 'Sn':
              number = 2.17
       elif atom == 'Sb':
              number = 2.06
       elif atom == 'I':
              number = 2.06
       elif atom == 'Cs':
              number = 3.43
       elif atom == 'Pt':
              number = 1.75
       elif atom == 'Au':
              number = 1.66
       elif atom == 'Hg':
              number = 1.65
       elif atom == 'Re':  # look this up
              number = 1.25
       elif atom == 'Pb':
              number = 2.02

       return 1.2*number

def getAtomicCovalentRadius(atom):

# the Average radii from the Wikipedia entry for "Covalent radius"

       if atom == 'H':
              number = 0.31
       elif atom == 'Li':
              number = 1.28
       elif atom == 'Be':
              number = 0.96
       elif atom == 'B':
              number = 0.84
       elif atom == 'C':
              number = 0.76
       elif atom == 'N':
              number = 0.71
       elif atom == 'O':
              number = 0.66
       elif atom == 'F':
              number = 0.57
       elif atom == 'Na':
              number = 1.67
       elif atom == 'Mg':
              number = 1.41
       elif atom == 'Al':
              number = 1.21
       elif atom == 'Si':
              number = 1.11
       elif atom == 'P':
              number = 1.07
       elif atom == 'S':
              number = 1.05
       elif atom == 'Cl':
       #       number = 0.4*1.32
              number = 0.8*1.32
       elif atom == 'Ar':
              number = 1.06
       elif atom == 'K':
              number = 2.03
       elif atom == 'Ca':
              number = 1.76
       elif atom == 'Mn':
              number = 1.61
       #       number = 1.00
       elif atom == 'Ti':
              number = 1.40
       elif atom == 'Fe':
              number = 1.32
       elif atom == 'Cu':
              number = 1.32
       elif atom == 'Zn':
              number = 1.22
       elif atom == 'Ga':
              number = 1.22
       elif atom == 'Ge':
              number = 1.20
       elif atom == 'As':
              number = 1.19
       elif atom == 'Br':
              number = 1.14
       elif atom == 'Rb':
              number = 2.20
       elif atom == 'Rh':
              number = 1.42
       elif atom == 'Ag':
              number = 1.45
       elif atom == 'Sn':
              number = 1.39
       elif atom == 'Sb':
              number = 1.39
       elif atom == 'I':
              number = 1.39
       elif atom == 'Cs':
              number = 2.44
       elif atom == 'Ce':
              number = 2.04
       elif atom == 'Sm':
              number = 1.98
       elif atom == 'Ir':
              number = 1.41
       elif atom == 'Pt':
              number = 1.38
       elif atom == 'Au':
              number = 1.36
       elif atom == 'Hg':
              number = 1.72
       elif atom == 'Re':    # look up
              number = 1.25
       elif atom == 'Pb':
              number = 1.46

       return number*1.16

def getAtomicNumber(atom):

       if atom == 'H':
              number = 1
       elif atom == 'Li':
              number = 3
       elif atom == 'Be':
              number = 4
       elif atom == 'B':
              number = 5
       elif atom == 'C':
              number = 6
       elif atom == 'N':
              number = 7
       elif atom == 'O':
              number = 8
       elif atom == 'F':
              number = 9
       elif atom == 'Na':
              number = 11
       elif atom == 'Mg':
              number = 12
       elif atom == 'Al':
              number = 13
       elif atom == 'Si':
              number = 14
       elif atom == 'P':
              number = 15
       elif atom == 'Cl':
              number = 17
       elif atom == 'Ar':
              number = 18
       elif atom == 'K':
              number = 19
       elif atom == 'Ca':
              number = 20
       elif atom == 'Mn':
              number = 25
       elif atom == 'Fe':
              number = 26
       elif atom == 'Cu':
              number = 29
       elif atom == 'Zn':
              number = 30
       elif atom == 'Ga':
              number = 31
       elif atom == 'Ge':
              number = 32
       elif atom == 'As':
              number = 33
       elif atom == 'Rb':
              number = 37
       elif atom == 'Ag':
              number = 47
       elif atom == 'Sn':
              number = 50
       elif atom == 'Sb':
              number = 51
       elif atom == 'I':
              number = 53
       elif atom == 'Cs':
              number = 55
       elif atom == 'Ce':
              number = 58
       elif atom == 'Sm':
              number = 62
       elif atom == 'Re':
              number = 75
       elif atom == 'Ir':
              number = 77
       elif atom == 'Au':
              number = 79
       elif atom == 'Hg':
              number = 80
       elif atom == 'Pb':
              number = 82
       elif atom == 'Cf':
              number = 98

       return number

def getAtomicSymbol(atom):

       #print atom

       table = []
       table.append('zip')
       table.append('H')
       table.append('He')
       table.append('Li')
       table.append('Be')
       table.append('B')
       table.append('C')
       table.append('N')
       table.append('O')
       table.append('F')
       table.append('Ne')
       table.append('Na')
       table.append('Mg')
       table.append('Al')
       table.append('Si')
       table.append('P')
       table.append('S')
       table.append('Cl')
       table.append('Ar')
       table.append('K')
       table.append('Ca')
       table.append('Sc')
       table.append('Ti')
       table.append('Y')
       table.append('Cr')
       table.append('Mn')
       table.append('Fe')
       table.append('Co')
       table.append('Ni')
       table.append('Cu')
       table.append('Zn')
       table.append('Ga')
       table.append('Ge')
       table.append('As')
       table.append('Se')
       table.append('Br')
       table.append('Kr')
       table.append('Rb')
       table.append('Sr')
       table.append('Y')
       table.append('Zr')
       table.append('Nb')
       table.append('Mo')
       table.append('Tc')
       table.append('Ru')
       table.append('Rh')
       table.append('Pd')
       table.append('Ag')
       table.append('Cd')
       table.append('In')
       table.append('Sn')
       table.append('Sb')
       table.append('Te')
       table.append('I')
       table.append('Xe')
       table.append('Cs')
       table.append('Ba')
       table.append('La')
       table.append('Ce')
       table.append('Pr')
       table.append('Nd')
       table.append('Pm')
       table.append('Sm')
       table.append('Eu')
       table.append('Gd')
       table.append('Tb')
       table.append('Dy')
       table.append('Ho')
       table.append('Er')
       table.append('Tm')
       table.append('Yb')
       table.append('Lu')
       table.append('Hf')
       table.append('Ta')
       table.append('W') 
       table.append('Re')
       table.append('Os')
       table.append('Ir')
       table.append('Pt')
       table.append('Au')
       table.append('Hg')
       table.append('Tl')
       table.append('Pb')
       table.append('Bi')
       table.append('Po')
       table.append('At')
       table.append('Rn')
       table.append('Fr')
       table.append('Ra')
       table.append('Ac')
       table.append('Th')
       table.append('Pa')
       table.append('U')
       table.append('Np')
       table.append('Pu')
       table.append('Am')
       table.append('Cm')
       table.append('Bk')
       table.append('Cf')


       return table[int(atom)]


def transposeList(m):
       c=list(zip(*m))
       return c

def rotateAboutX(matrix, angle):

       import numpy
       import math
       radangle=float(angle)*math.pi/180

       transMatrix=numpy.zeros((3,3), float)
       transMatrix[0][0]=1.0
       transMatrix[1][1]=math.cos(angle)
       transMatrix[1][2]=math.sin(angle)
       transMatrix[2][1]=-1.0*math.sin(angle)
       transMatrix[2][2]=math.cos(angle)

       matrix=numpy.array(matrix)
       newmatrix=matrixmultiply(matrix, transMatrix)
       newmatrix=arrayToList(newmatrix)

       return newmatrix

def rotateAboutY(matrix, angle):

       import numpy
       import math
       radangle=float(angle)*math.pi/180

       transMatrix=numpy.zeros((3,3), float)
       transMatrix[1][1]=1.0
       transMatrix[0][0]=math.cos(angle)
       transMatrix[2][0]=math.sin(angle)
       transMatrix[0][2]=-1.0*math.sin(angle)
       transMatrix[2][2]=math.cos(angle)

       matrix=numpy.array(matrix)
       newmatrix=matrixmultiply(matrix, transMatrix)
       newmatrix=arrayToList(newmatrix)

       return newmatrix       

def rotateAboutZ(matrix, angle):

       import numpy
       import math

       radangle=float(angle)*math.pi/180

       transMatrix=numpy.zeros((3,3), float)
       transMatrix[2][2]=1.0
       transMatrix[0][0]=math.cos(angle)
       transMatrix[0][1]=math.sin(angle)
       transMatrix[1][0]=-1.0*math.sin(angle)
       transMatrix[1][1]=math.cos(angle)

       matrix=numpy.array(matrix)
       newmatrix=matrixmultiply(matrix, transMatrix)
       newmatrix=arrayToList(newmatrix)

       return newmatrix       

def arrayToList(thearray):

        import numpy

        shape=numpy.shape(thearray)
        numatoms=shape[0]

        list=[]
        for i in range(numatoms):
                list.append([])
                for j in range(shape[1]):
                        list[i].append(thearray[i][j])

        return list

def writeXYZ(coords, name, atomtypes):

       filename = name
       outfile = open(filename,'w')
       numIons=len(coords)
       
       outfile.write(str(numIons)+'\n')
       outfile.write(name+'\n')
       for j in range(numIons):
              outfile.write('    '+str(atomtypes[j])+'     '+'%.8f'%(coords[j][0])+'     '+str('%.8f'%(coords[j][1]))+'     '+str('%.8f'%(coords[j][2]))+'  '+'\n')


def readXYZMovie(inputfilename):

        import os

        inputfile=open(inputfilename, 'r')
        numcoords=0
        coords=[]
        atomtypes=[]
        numatoms=int(inputfile.readline().strip())
        numStructs = int(int(os.popen('wc -l '+ inputfilename + " | awk '{print $1}'").read().strip())/(numatoms+2))
        toprint = 'There are '+str(numStructs)+' structures being read from '+inputfilename
        print('')
        inputfile.readline()
        
        coords.append([])
        for z in range(numatoms):
              line=inputfile.readline()
              toadd=[float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])]
              coords[0].append(toadd)
              atomtypes.append(line.strip().split()[0])
        
        for z2 in range(numStructs-1):
              inputfile.readline()
              inputfile.readline()
        
              coords.append([])
              for z1 in range(numatoms):
                 line=inputfile.readline()
                 toadd=[float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])]
                 coords[z2+1].append(toadd)
 
        inputfile.close()

        return [coords,atomtypes] 

def getAngle(u,v):

        import math

        udotv = u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

        umod = math.sqrt((u[0]*u[0])+(u[1]*u[1])+(u[2]*u[2]))
        vmod = math.sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]))
        
        oper=udotv*(1/umod)*(1/vmod)

        thetaR=math.acos(oper) 

        thetaD=thetaR*180.0/math.pi

        return thetaR,thetaD

def getAtomCoord(coords,n):

        atomc=[]
        for a in range(coords[n-1]):
            atomc.append(coords[n-1][a])

        return atomc

def getAtomSetCoords(coords,nlist):

        aset=[]
        for aT in range(len(nlist)):
            aset.append(coords[nlist[aT]-1][:])
 
        return aset

def allFracToCart(coords, cell):

        numatoms=len(coords)
        cartcoords=[]

        for i in range(numatoms):
                atom=coords[i]
                cartatom=fracToCart(atom, cell)
                cartcoords.append(cartatom)

        return cartcoords


def allCartToFrac(coords, cell):

        numatoms=len(coords)
        fraccoords=[]

        for i in range(numatoms):
                atom=coords[i]
                fracatom=cartToFrac(atom, cell)
                fraccoords.append(fracatom)

        return fraccoords

def fracToCart(vector, cell):

        import numpy

        cell_T = transposeList(cell)
        cartcoords = matrixmultiply(numpy.array(cell_T),numpy.array(vector))

        return [cartcoords[0],cartcoords[1],cartcoords[2]]


def cartToFrac(vector, cell):
   
        import numpy

        cell_T = transposeList(cell)
        cell_I = numpy.linalg.inv(numpy.array(cell_T))
        newfrac = matrixmultiply(cell_I,numpy.array(vector))

        return [newfrac[0],newfrac[1],newfrac[2]]

def findMolecules(coords, cell, atomtypes, dists):

        import math
        import numpy

        #print('initial atomtypes :', atomtypes)
        natoms=len(coords)
        molList=[]
        pruneList=[]

        for i in range(natoms):
                  assign=0
                  bounddict=[]
                  # build a list, bounddict, of atoms bound to atom i
                  for z in range(i):
                       covR=getAtomicCovalentRadius(atomtypes[i])+getAtomicCovalentRadius(atomtypes[z])
                       if dists[i][z]<covR:
                            if dists[i][z]<0.1:
                                print('atoms are too close: ',i,' and ',z) 
                                print('                       removing ',i) 
                                pruneList.append(i)
                                assign=-1
                            else:
                                bounddict.append(z)
                                assign=1 


                  if assign==1:
                      # make a list, tojoin, of the elements of molList to which atom i should belong
                      counter=0
                      tojoin=[]
                      for a in molList:
                            for q in bounddict:
                                if q in a:
                                      tojoin.append(counter)
                                      break
                            
                            counter+=1 
                                    
                      # merge molList entries together
                      tojoin.sort()
                      basem=tojoin.pop(0)
                      molList[basem].append(i)
                      tojoin.reverse()
                      for a in tojoin:
                           cut=molList.pop(a)
                           for a1 in cut:
                              if a1 not in molList[basem]:
                                      molList[basem].append(a1)

                  # create a new molList entry for i               
                  elif assign==0:
                      molList.append([i])

        # reorder the coords and atomtypes files
        newcoords=[]
        newatomtypes=[]
        numpermol=[]
        for a in molList:
            a.sort()
            numpermol.append(len(a))
            for a1 in a:
               newcoords.append(coords[a1][:])
               newatomtypes.append(atomtypes[a1][:])

        #print('newatomtypes :', newatomtypes)
        return newcoords, cell, newatomtypes, numpermol


def calcRadialPOSCAR(coords, cell):

        import math
        import numpy

        natoms=len(coords)
        dists=[]
        for b in range(natoms):
             dists.append(natoms*[0.0]) 

        distcut=10.0
        cella=cell[0][:]
        la=int(distcut/length(cella))
        z1r=[float(i-1-la) for i in range(2*la+3)]
        cellb=cell[1][:]
        lb=int(distcut/length(cellb))
        z2r=[float(i-1-lb) for i in range(2*lb+3)]
        cellc=cell[2][:]
        lc=int(distcut/length(cellc))
        z3r=[float(i-1-lc) for i in range(2*lc+3)]

        for i in range(natoms):
                        candidate = [coords[i][0], coords[i][1], coords[i][2]]
                        for z in range(i):
                                dist=1000.0
                                for z1 in z1r:
                                        for z2 in z2r:
                                                for z3 in z3r:
                                                        image = add(coords[z][:],[z1,z2,z3])
                                                        testdiff = subtract(candidate,image)
                                                        testdiff2 = fracToCart(testdiff, cell)
                                                        testdist = length(testdiff2) 
                                                        if testdist<dist:
                                                               dist=testdist
                                dists[i][z]=dist
                                dists[z][i]=dist

        return dists

def getNNmolList(coords,atomtypes,cell,dists,ncomp):

        import math
        import numpy

        nmols=len(ncomp)
        natoms=len(coords)
        nnlist=[]
        centermass=[]
        mollookup=[]
        a=0
        m1=0
        for j1 in ncomp:
               molave=[0.0,0.0,0.0]
               for a1 in range(j1):
                   molave[0]=molave[0]+coords[a1+a][0]
                   molave[1]=molave[1]+coords[a1+a][1]
                   molave[2]=molave[2]+coords[a1+a][2]
                   mollookup.append(m1)

               molave[0]=molave[0]/float(j1)
               molave[1]=molave[1]/float(j1)
               molave[2]=molave[2]/float(j1)
               a=a+j1
               m1+=1

               centermass.append(molave)
               nnlist.append([])

        shift=ncomp[0]
        for a1 in range(nmols-1):
              nmatoms=ncomp[a1+1]
              alist=list(range(shift,nmatoms+shift))
              olist=list(range(0,shift))
              for a2 in alist:
                  for a3 in olist:
                         abR=getVDWRadius(atomtypes[a2])+getVDWRadius(atomtypes[a3])
                         if dists[a2][a3]<abR: 
                                 imap=getShortestTranslate(coords[a3][:],coords[a2][:],cell)
                                 if [mollookup[a2],imap] not in nnlist[mollookup[a3]]:
                                      nnlist[mollookup[a3]].append([mollookup[a2],imap])
                  for b1 in alist:
                         abR=getVDWRadius(atomtypes[a2])+getVDWRadius(atomtypes[b1])
                         imaplist=getVDWTranslations(coords[a2][:],coords[b1][:],cell,abR)
                         for b2 in imaplist:
                                 cand=[mollookup[b1],b2]
                                 if cand not in nnlist[mollookup[a2]]:
                                      nnlist[mollookup[a2]].append(cand)

              shift+=nmatoms

        return nnlist 
                  
def printNNNList(coords,atomtypes,cell,ncomp,nnlist,mol):

        import math
        import numpy
        #print(nnlist[mol])

        nmols=len(nnlist)
        nnflist=[[] for i in range(nmols)]
        for a in range(nmols):
             for b in nnlist[a]:
                 nnflist[a].append(b)
                 if b[0]!=a:
                    nnflist[b[0]].append([a,[float(-1.0*i) for i in b[1]]])

        #print(nnflist[mol])

        # center molecule - called cm
        numcm=ncomp[mol]
        shift=sum([i for i in ncomp[:mol]])
        cmlist=list(range(shift,numcm+shift))
        cmcoords=makeCopyList2D(coords[shift:numcm+shift][:])
        cmxyz=allFracToCart(cmcoords,cell)
        cmatypes=atomtypes[shift:numcm+shift][:]

        n=len(nnflist[mol])
        for n1 in range(n):
             # take paired m as first end-point - is a[0]
             a=nnflist[mol][n1]
             numa=ncomp[a[0]]
             alist=list(range(sum(ncomp[:a[0]]),sum(ncomp[:a[0]+1])))

             rest=[i+(n1+1) for i in list(range(n-(n1+1)))]
             for n2 in rest:
                # take n2 as second end-point - is b[0]
                b=nnflist[mol][n2]
                numb=ncomp[b[0]]
                blist=list(range(sum(ncomp[:b[0]]),sum(ncomp[:b[0]+1])))

                # 'p' is the triplet list of molecules - starts from center
                # 'e' is the doublet list of end molecules 
                pcoords=makeCopyList2D(cmxyz)
                patypes=cmatypes[:]
                ecoords=[]
                eatypes=[]

                for a1 in alist:
                   coord=add(coords[a1][:], a[1])
                   coordxyz=fracToCart(coord,cell)
                   pcoords.append(coordxyz)
                   patypes.append(atomtypes[a1])
                   ecoords.append(coordxyz)
                   eatypes.append(atomtypes[a1])

                for b1 in blist:
                   coord=add(coords[b1][:], b[1])
                   coordxyz=fracToCart(coord,cell)
                   pcoords.append(coordxyz)
                   patypes.append(atomtypes[b1])
                   ecoords.append(coordxyz)
                   eatypes.append(atomtypes[b1])


                tagp='mol'+str(a[0])+'-'+str(mol)+'-'+str(b[0])+'_a'+str(int(a[1][0]))+str(int(a[1][1]))+str(int(a[1][2]))+'_b'+str(int(b[1][0]))+str(int(b[1][1]))+str(int(b[1][2]))+'.xyz'
                tage='mol'+str(a[0])+'-'+str(mol)+'-'+str(b[0])+'_a'+str(int(a[1][0]))+str(int(a[1][1]))+str(int(a[1][2]))+'_b'+str(int(b[1][0]))+str(int(b[1][1]))+str(int(b[1][2]))+'ENDS.xyz'
                writeXYZ(pcoords,tagp,patypes)
                writeXYZ(ecoords,tage,eatypes)

def printMolList(coords,atomtypes,cell,ncomp,nnlist):

        import math
        import numpy

        shift=0
        for a1 in range(len(ncomp)):
             numm1=ncomp[a1]
             m1list=list(range(shift,numm1+shift))
             m1coords=makeCopyList2D(coords[shift:numm1+shift][:])
             m1coordsxyz=allFracToCart(m1coords,cell)
             m1atypes=atomtypes[shift:numm1+shift][:]
             tag='mol'+str(a1)+'.xyz'
             writeXYZ(m1coordsxyz,tag,m1atypes)
             shift=shift+numm1



def printNNList(coords,atomtypes,cell,ncomp,nnlist):

        import math
        import numpy

        shift=0
        for a1 in range(len(nnlist)):
             numm1=ncomp[a1]
             m1list=list(range(shift,numm1+shift))
             m1coords=makeCopyList2D(coords[shift:numm1+shift][:])
             m1coordsxyz=allFracToCart(m1coords,cell)
             m1atypes=atomtypes[shift:numm1+shift][:]
             for a2 in nnlist[a1]:
                 tag='mol'+str(a1)+'-'+str(a2[0])+'_'+str(int(a2[1][0]))+str(int(a2[1][1]))+str(int(a2[1][2]))+'.xyz'
                 numm2=ncomp[a2[0]]
                 m2list=list(range(sum(ncomp[:a2[0]]),sum(ncomp[:a2[0]+1])))
                 pcoords=makeCopyList2D(m1coordsxyz)
                 patypes=m1atypes[:]
                 for a3 in m2list:
                       coord=add(coords[a3][:], a2[1])
                       coordxyz=fracToCart(coord,cell)
                       pcoords.append(coordxyz)
                       patypes.append(atomtypes[a3])

                 writeXYZ(pcoords,tag,patypes)
             shift=shift+numm1


def makeCohesive(coords,cell,atomtypes,dists,ncomp):

        import math
        import numpy

        nmols=len(ncomp)
        newcoords=makeCopyList2D(coords)
        shift=0
        for i in range(nmols):
               natoms=ncomp[i]
               alist=list(range(shift,natoms+shift))
               active=[alist[0]]
               while len(active)>0:
                  for i1 in active:
                     alist.remove(i1)
                  newactive=[]
                  for i2 in active:
                     for z in alist:
                          covR=getAtomicCovalentRadius(atomtypes[i2])+getAtomicCovalentRadius(atomtypes[z])
                          if dists[i2][z]<covR and z not in newactive:
                               newactive.append(z)
                               newcoords[z]=getShortest(newcoords[i2],newcoords[z],cell)
                  active=newactive

               shift+=natoms

        return newcoords 
                  

def getVDWTranslations(v1,v2,cell,abR):

        import math
        import numpy

        cella=cell[0][:]
        la=int(10.0/length(cella))
        z1r=[float(i-1-la) for i in range(2*la+3)]
        cellb=cell[1][:]
        lb=int(10.0/length(cellb))
        z2r=[float(i-1-lb) for i in range(2*lb+3)]
        cellc=cell[2][:]
        lc=int(10.0/length(cellc))
        z3r=[float(i-1-lc) for i in range(2*lc+3)]

        dist=1000.0
        newv2=[0.,0.,0.]
        imaplib=[]
        for z1 in z1r:
               for z2 in z2r:
                        for z3 in z3r:
                                 image = add(v2,[z1,z2,z3])
                                 testdiff = subtract(v1,image)
                                 testdiff2 = fracToCart(testdiff, cell)
                                 testdist = length(testdiff2) 
                                 if testdist<abR:
                                          idist=length([z1,z2,z3])
                                          if idist>0.02:
                                              imap=[z1,z2,z3]
                                              imaplib.append(imap)

        return imaplib


def getShortestTranslate(v1,v2,cell):

        import math
        import numpy

        cella=cell[0][:]
        la=int(10.0/length(cella))
        z1r=[float(i-1-la) for i in range(2*la+3)]
        cellb=cell[1][:]
        lb=int(10.0/length(cellb))
        z2r=[float(i-1-lb) for i in range(2*lb+3)]
        cellc=cell[2][:]
        lc=int(10.0/length(cellc))
        z3r=[float(i-1-lc) for i in range(2*lc+3)]

        dist=1000.0
        newv2=[0.,0.,0.]
        for z1 in z1r:
               for z2 in z2r:
                        for z3 in z3r:
                                 image = add(v2,[z1,z2,z3])
                                 testdiff = subtract(v1,image)
                                 testdiff2 = fracToCart(testdiff, cell)
                                 testdist = length(testdiff2) 
                                 if testdist<dist:
                                          dist=testdist
                                          imap=[z1,z2,z3]

        return imap

def getShortest(v1,v2,cell):

        import math
        import numpy

        cella=cell[0][:]
        la=int(15.0/length(cella))
        z1r=[float(i-1-la) for i in range(2*la+3)]
        cellb=cell[1][:]
        lb=int(15.0/length(cellb))
        z2r=[float(i-1-lb) for i in range(2*lb+3)]
        cellc=cell[2][:]
        lc=int(15.0/length(cellc))
        z3r=[float(i-1-lc) for i in range(2*lc+3)]

        dist=1000.0
        newv2=[0.,0.,0.]
        for z1 in z1r:
               for z2 in z2r:
                        for z3 in z3r:
                                 image = add(v2,[z1,z2,z3])
                                 testdiff = subtract(v1,image)
                                 testdiff2 = fracToCart(testdiff, cell)
                                 testdist = length(testdiff2) 
                                 if testdist<dist:
                                          dist=testdist
                                          newv2=[image[0],image[1],image[2]]

        return newv2 

def makeCopyList2D(l):

        import math
        import numpy

        newl=[]
        nd1=len(l)
        for z1 in range(nd1):
                newl.append([])
                nd2=len(l[z1])
                for z2 in range(nd2):
                      newl[z1].append(l[z1][z2])

        return newl

def makeCopyList1D(l):

        import math
        import numpy

        newl=[]
        nd1=len(l)
        for z1 in range(nd1):
               newl.append(l[z1])

        return newl

def reordermols(coords, atype, ncomp, order):

       # the new set of coordinates and atomtypes
       newcoords=[]
       newatypes=[]
       newncomp=[]
       print('hello')
       print(ncomp)

       for i in order:
           a=sum(ncomp[:i])
           b=a+ncomp[i]         
           print(b)
           newncomp.append(ncomp[i])
           print([a+m for m in range(b-a)])
           for j in [a+m for m in range(b-a)]:
              newcoords.append(makeCopyList1D(coords[j]))
              newatypes.append(atype[j])
           print(newatypes)

       return newcoords, newatypes, newncomp

def writePoscar(coords, cell, name, atomtype):

        atype=makeCopyList1D(atomtype)

        filename = name
        outfile = open(filename,'w')
        numIons=len(coords)

        asymbols='  '

        count = 1
        lasts=atype.pop(0)
        numspecies=[]

        for i in atype:
              if i==lasts:
                  count+=1
              else:
                  asymbols=asymbols+lasts+' '
                  numspecies.append(count)
                  lasts=i
                  count=1
                  
        asymbols=asymbols+lasts
        numspecies.append(count)
            

        outfile.write(asymbols+'  \n')
        outfile.write('1.000000'+'\n')
        for i1 in range(3):
                outfile.write('    '+'%13.10f'%(cell[i1][0])+'   '+'%13.10f'%(cell[i1][1])+'   '+'%13.10f'%(cell[i1][2])+'\n')
        outfile.write(asymbols+'  \n')
        for i2 in range(len(numspecies)):
                outfile.write(str(numspecies[i2])+' ')
        outfile.write('\n')
        outfile.write('Direct'+'\n')
        for j in range(numIons):
                outfile.write('   '+'%13.10f'%(coords[j][0])+'   '+'%13.10f'%(coords[j][1])+'   '+'%13.10f'%(coords[j][2])+'\n')

