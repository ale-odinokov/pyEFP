# this file contains definitions of classes and functions aimed at the 
# reading, writing and modification of the data from the EFP files obtained
# as a result of the GAMESS US calculation

import openbabel
import pybel
import numpy as np

Bohr = 0.52917721092 

NumValElec = {"H":1, "C":4, "N":5, "O":6, "S":6}
Symb = ("H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar")

# Converts string from the EFP file to OBAtom object
def EFP2atom(s):
    at=openbabel.OBAtom()
    st=s.split()
    at.SetAtomicNum(int(float(st[-1])))
    at.SetVector(float(st[1])*Bohr, float(st[2])*Bohr, float(st[3])*Bohr)
    at.SetTitle(st[0][3:len(st[0])])
    return at

# Reads EFP file and returns OBMol object
def EFP2mol(fname):
    mol0=openbabel.OBMol()
    fp=open(fname,"r")
    s=fp.readline()
    while s:
        if (s.strip()=="$FRAGNAME"):
            break
        s=fp.readline()
    if not s:
        return None
    s=fp.readline()
    s=fp.readline()
    s=fp.readline()
    while s and s[0:1]=="A":
#        print s,
        at=EFP2atom(s)
        mol0.AddAtom(at)
        s=fp.readline()
    fp.close()

    mol0.ConnectTheDots()
    mol0.PerceiveBondOrders()

    return mol0

# Writes block of coordinates in GAMESS US format from OBMol object
def WriteGAMESS_DATA(fp, mol):
    for at in openbabel.OBMolAtomIter(mol):
        fp.write(" {0:3s}{1:3.1f}{2:18.8f}{3:18.8f}{4:18.8f}\n".format(Symb[at.GetAtomicNum()-1], float(at.GetAtomicNum()), at.x(), at.y(), at.z()))

# Writes block of coordinates in XYZ format from OBMol object
def WriteOBMol2xyz(fp, mol):
    for at in openbabel.OBMolAtomIter(mol):
        fp.write(" {0:3s}{1:18.8f}{2:18.8f}{3:18.8f}\n".format(Symb[at.GetAtomicNum()-1], at.x(), at.y(), at.z()))

# Writes block of coordinates in PDB format from OBMol object
def WriteOBMol2PDB(fp, mol, start, resid):
    i = 0
    for at in openbabel.OBMolAtomIter(mol):
            fp.write("ATOM  {0:<5s} {1:4s} LIG {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f}{6:6.2f}{7:6.2f}          {8:>2s}  \n".format(str(i+start)[-5:], \
                    str(Symb[at.GetAtomicNum()-1]).ljust(4)[0:4], resid, at.x(), at.y(), at.z(), 0.0, 0.0, Symb[at.GetAtomicNum()-1]))
            i += 1

# Best fit of two groups of atoms defined via OBMol objects and lists of indices
# Procedure follows Kabsch algorithm as described in
# W. Kabsch, Acta Cryst. (1978). A34, 827-828
# It returns COM vectors for the given and reference groups and rotation (matrix)

def SolveKabsch(molFit, molRef, maskFit, maskRef):
    # Extracting coordinates and masses of the reference group
    massRef=[]
    coordsRef=[]
    for i in maskRef:
        at=molRef.GetAtom(i)
        massRef.append(at.GetAtomicMass())
        coordsRef.append([at.GetX(), at.GetY(), at.GetZ()])

    # Extracting coordinates and masses of the fitted group
    massFit=[]
    coordsFit=[]
    for i in maskFit:
        at=molFit.GetAtom(i)
        massFit.append(at.GetAtomicMass())
        coordsFit.append([at.GetX(), at.GetY(), at.GetZ()])

    # Defining Center of Mass (COM) position of the reference group
    mRef=0.0
    COMRef=[0.0, 0.0, 0.0]
    for i in xrange(len(massRef)):
        mRef+=massRef[i]
        COMRef[0]+=massRef[i]*coordsRef[i][0]
        COMRef[1]+=massRef[i]*coordsRef[i][1]
        COMRef[2]+=massRef[i]*coordsRef[i][2]
    COMRef=np.array(COMRef)/mRef

    # Defining Center of Mass (COM) position of the fitted group
    mFit=0.0
    COMFit=[0.0, 0.0, 0.0]
    for i in xrange(len(massFit)):
        mFit+=massFit[i]
        COMFit[0]+=massFit[i]*coordsFit[i][0]
        COMFit[1]+=massFit[i]*coordsFit[i][1]
        COMFit[2]+=massFit[i]*coordsFit[i][2]
    COMFit=np.array(COMFit)/mFit

    # Shifting reference group
    crdRef=[]
    for v in coordsRef:
        crdRef.append([v[0]-COMRef[0], v[1]-COMRef[1], v[2]-COMRef[2]])

    # Shifting fitted group
    crdFit=[]
    for v in coordsFit:
        crdFit.append([v[0]-COMFit[0], v[1]-COMFit[1], v[2]-COMFit[2]])

    # Calculating R matrix
    Rmat=[]
    for i in xrange(3):
        Rmat.append([])
        for j in xrange(3):
            temp=0
            for k in xrange(len(massFit)):
                if (massFit[k]!=massRef[k]):
                    print "WARNING: masses of atoms from the reference and fitted groups are not pairwise identical"
                temp += massFit[k]*crdRef[k][i]*crdFit[k][j]
            Rmat[i].append(temp)
    Rmat=np.array(Rmat)

    # Calculating RR matrix
    RRmat = np.dot(Rmat.T,Rmat)

    # Calculating vectors a
    mu, a = np.linalg.eig(RRmat)
    if(np.linalg.det(a)<0):
        for i in xrange(3):
            a[i][2] = -1*a[i][2]
    a = a.T
    sol = sorted(zip(mu, a, [0,1,2]), key = lambda x: x[0], reverse=True)
    mu = np.array(zip(*sol)[0])
    a = np.array(zip(*sol)[1])
    order = np.array(zip(*sol)[2])

    # Calculating vectors b
    b=[]
    for i in xrange(2):
        b.append([])
        for j in xrange(3):
            temp=0
            for k in xrange(3):
                temp += Rmat[j][k]*a[i][k]/np.sqrt(mu[i])
            b[i].append(temp)
    b.append([])
    b[2].append(b[0][1]*b[1][2]-b[0][2]*b[1][1])
    b[2].append(b[0][2]*b[1][0]-b[0][0]*b[1][2])
    b[2].append(b[0][0]*b[1][1]-b[0][1]*b[1][0])
    b = np.array(b)

    newa = [1, 2, 3]
    for i in xrange(len(a)):
        newa[order[i]] = a[i]
    newa = np.array(newa)
    newb = [1, 2, 3]
    for i in xrange(len(b)):
        newb[order[i]] = b[i]
    newb = np.array(newb)

    # Calculating rotation matrix
    rot = np.dot(newb.T,newa)

    return rot, COMRef, COMFit


# Applying least squares fit to the OBMol object 
def ApplyKabsch(molFit, COMFit, COMRef, rot):
    crd=[]
    for at in openbabel.OBMolAtomIter(molFit):
        crd.append([at.GetX()-COMFit[0], at.GetY()-COMFit[1], at.GetZ()-COMFit[2]])
    crd=np.array(crd)
    for i in xrange(len(crd)):
        crd[i]=np.dot(rot,crd[i])
    i = 0
    for at in openbabel.OBMolAtomIter(molFit):
        at.SetVector(crd[i][0]+COMRef[0], crd[i][1]+COMRef[1], crd[i][2]+COMRef[2])
        i+=1
    return True

# Applying least squares fit to the EFPdata object 
# multipole definitions are erased
def ApplyKabschEFP(molFit, COMFit, COMRef, rot):
    crd=[]
    for at in molFit.atoms:
        crd.append([Bohr*at.x-COMFit[0], Bohr*at.y-COMFit[1], Bohr*at.z-COMFit[2]])
    for b in molFit.bonds:
        crd.append([Bohr*b.x-COMFit[0], Bohr*b.y-COMFit[1], Bohr*b.z-COMFit[2]])
    for c in molFit.centroids:
        crd.append([Bohr*c.x-COMFit[0], Bohr*c.y-COMFit[1], Bohr*c.z-COMFit[2]])
    crd=np.array(crd)
    for i in xrange(len(crd)):
        crd[i]=np.dot(rot,crd[i])
    for i in xrange(len(crd)):
        crd[i]=[crd[i][0]+COMRef[0], crd[i][1]+COMRef[1], crd[i][2]+COMRef[2]]
    i = 0
    newMol = EFPdata()
    newMol.name = molFit.name
    for i in xrange(molFit.natoms):
        newMol.natoms += 1
        newMol.atoms.append(EFPatom())
        newMol.atoms[i].name = molFit.atoms[i].name
        newMol.atoms[i].x = crd[i][0]/Bohr
        newMol.atoms[i].y = crd[i][1]/Bohr
        newMol.atoms[i].z = crd[i][2]/Bohr
        newMol.atoms[i].mass = molFit.atoms[i].mass
        newMol.atoms[i].num = molFit.atoms[i].num
        newMol.atoms[i].charge = molFit.atoms[i].charge
    for i in xrange(molFit.nbonds):   
        newMol.nbonds += 1
        newMol.bonds.append(EFPbond())
        newMol.bonds[i].idx1 = molFit.bonds[i].idx1
        newMol.bonds[i].idx2 = molFit.bonds[i].idx2
        newMol.bonds[i].x = crd[i+molFit.natoms][0]/Bohr
        newMol.bonds[i].y = crd[i+molFit.natoms][1]/Bohr
        newMol.bonds[i].z = crd[i+molFit.natoms][2]/Bohr
        newMol.bonds[i].charge = molFit.bonds[i].charge 
    for i in xrange(molFit.ncentroids): 
        newMol.ncentroids += 1
        newMol.centroids.append(EFPcentroid())
        newMol.centroids[i].x = crd[i+molFit.natoms+molFit.nbonds][0]/Bohr
        newMol.centroids[i].y = crd[i+molFit.natoms+molFit.nbonds][1]/Bohr
        newMol.centroids[i].z = crd[i+molFit.natoms+molFit.nbonds][2]/Bohr
    return newMol

# -------------------------------- Definition of classes -------------------------------------------------

# Class desribing Atom entry in the EFP
class EFPatom:
    def __init__(self):
        self.name = "NONAME"
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.mass = 0.0
        self.num = 0.0
        self.charge = 0.0
        self.dipole = (0.0, 0.0, 0.0)
        self.quadrupole = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.octupole = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.screen2 = (0.0, 0.0)
        self.screen3 = (0.0, 0.0)

    def FillFromString(self, s):
        st = s.split()
        self.name = st[0][3:len(st[0])]
        self.x = float(st[1])
        self.y = float(st[2])
        self.z = float(st[3])
        self.mass = float(st[4])
        self.num = float(st[5])

    def AddCharge(self, s):
        st = s.split()
        self.charge = float(st[1])

    def AddDipole(self, s):
        st = s.split()
        self.dipole = ( float(st[1]), float(st[2]), float(st[3]) )

    def AddQuadrupole(self, s1, s2):
        st1 = s1.split()
        st2 = s2.split()
        self.quadrupole = (float(st1[1]), float(st1[2]), float(st1[3]), float(st1[4]), float(st2[0]), float(st2[1]))

    def AddOctupole(self, s1, s2, s3):
        st1 = s1.split()
        st2 = s2.split()
        st3 = s3.split()
        self.octupole = ( float(st1[1]), float(st1[2]), float(st1[3]), float(st1[4]), float(st2[0]),
                          float(st2[1]), float(st2[2]), float(st2[3]), float(st3[0]), float(st3[1]) )

    def AddScreen2(self, s):
        st = s.split()
        self.screen2 = ( float(st[1]), float(st[2]) )

    def AddScreen3(self, s):
        st = s.split()
        self.screen3 = ( float(st[1]), float(st[2]) )

    # Print to the screen in the human-readable format
    def Print(self):
        print "Name: {0:10s}".format(self.name)
        print "Number: {0:3.1f}".format(self.num)
        print "Mass: {0:10.7f}".format(self.mass)
        print "x: {0:15.10f}".format(self.x)
        print "y: {0:15.10f}".format(self.y)
        print "z: {0:15.10f}".format(self.z)
        print "Charge: {0:15.10f}".format(self.charge)
        print "Dipole: {0:16.10f}{1:16.10f}{2:16.10f}".format(self.dipole[0], self.dipole[1], self.dipole[2])
        print "Quadrupole: {0:16.10f}{1:16.10f}{2:16.10f}{3:16.10f} >"\
              .format(self.quadrupole[0], self.quadrupole[1], self.quadrupole[2], self.quadrupole[3])
        print "            {0:16.10f}{1:16.10f}".format(self.quadrupole[4], self.quadrupole[5])
        print "Octupole: {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >"\
              .format(self.octupole[0], self.octupole[1], self.octupole[2], self.octupole[3])
        print "          {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >"\
              .format(self.octupole[4], self.octupole[5], self.octupole[6], self.octupole[7])
        print "          {0:17.9f}{1:17.9f}".format(self.octupole[8], self.octupole[9])
        print "Screen2:{0:14.9f}{1:14.9f}".format(self.screen2[0], self.screen2[1])
        print "Screen3:{0:14.9f}{1:14.9f}".format(self.screen3[0], self.screen3[1])

# -------------------------------------------------------------------------------------------------------
# Class desribing Bond entry in the EFP
class EFPbond:
    def __init__(self):
        self.idx1 = 0
        self.idx2 = 0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.charge = 0.0
        self.dipole = (0.0, 0.0, 0.0)
        self.quadrupole = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.octupole = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.screen2 = (0.0, 0.0)
        self.screen3 = (0.0, 0.0)

    def FillFromString(self, s):
        st = s.split()
        ss = st[0]
        ss = ss[2:len(ss)]
        if (len(ss)==2):
            self.idx1 = int(ss[0])
            self.idx2 = int(ss[1])
        if (len(ss)==3):
            self.idx1 = int(ss[0:2])
            self.idx2 = int(ss[2])
        if (len(ss)==4):
            self.idx1 = int(ss[0:2])
            self.idx2 = int(ss[2:4])
        self.x = float(st[1])
        self.y = float(st[2])
        self.z = float(st[3])

    def AddCharge(self, s):
        st = s.split()
        self.charge = float(st[1])

    def AddDipole(self, s):
        st = s.split()
        self.dipole = ( float(st[1]), float(st[2]), float(st[3]) )

    def AddQuadrupole(self, s1, s2):
        st1 = s1.split()
        st2 = s2.split()
        self.quadrupole = (float(st1[1]), float(st1[2]), float(st1[3]), float(st1[4]), float(st2[0]), float(st2[1]))

    def AddOctupole(self, s1, s2, s3):
        st1 = s1.split()
        st2 = s2.split()
        st3 = s3.split()
        self.octupole = ( float(st1[1]), float(st1[2]), float(st1[3]), float(st1[4]), float(st2[0]),
                          float(st2[1]), float(st2[2]), float(st2[3]), float(st3[0]), float(st3[1]) )

    def AddScreen2(self, s):
        st = s.split()
        self.screen2 = ( float(st[1]), float(st[2]) )

    def AddScreen3(self, s):
        st = s.split()
        self.screen3 = ( float(st[1]), float(st[2]) )

    # Print to the screen in the human-readable format
    def Print(self):
        print "index 1: {0:10d}".format(self.idx1)
        print "index 2: {0:10d}".format(self.idx2)
        print "x: {0:15.10f}".format(self.x)
        print "y: {0:15.10f}".format(self.y)
        print "z: {0:15.10f}".format(self.z)
        print "Charge: {0:15.10f}".format(self.charge)
        print "Dipole: {0:15.10f}  {1:15.10f}  {2:15.10f}".format(self.dipole[0], self.dipole[1], self.dipole[2])
        print "Quadrupole: {0:16.10f}{1:16.10f}{2:16.10f}{3:16.10f} >"\
              .format(self.quadrupole[0], self.quadrupole[1], self.quadrupole[2], self.quadrupole[3])
        print "            {0:16.10f}{1:16.10f}".format(self.quadrupole[4], self.quadrupole[5])
        print "Octupole: {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >"\
              .format(self.octupole[0], self.octupole[1], self.octupole[2], self.octupole[3])
        print "          {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >"\
              .format(self.octupole[4], self.octupole[5], self.octupole[6], self.octupole[7])
        print "          {0:17.9f}{1:17.9f}".format(self.octupole[8], self.octupole[9])
        print "Screen2:{0:14.9f}{1:14.9f}".format(self.screen2[0], self.screen2[1])
        print "Screen3:{0:14.9f}{1:14.9f}".format(self.screen3[0], self.screen3[1])

# -------------------------------------------------------------------------------------------------------
# Class desribing Centroid entry in the EFP
class EFPcentroid:
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.polar = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    def FillFromString(self, s):
        st = s.split()
        self.x = float(st[1])
        self.y = float(st[2])
        self.z = float(st[3])

    def AddPolar(self, s1, s2, s3):
        st1 = s1.split()
        st2 = s2.split()
        st3 = s3.split()
        self.polar = ( float(st1[0]), float(st1[1]), float(st1[2]), float(st1[3]), float(st2[0]),
                          float(st2[1]), float(st2[2]), float(st2[3]), float(st3[0]) )

    # Print to the screen in the human-readable format
    def Print(self):
        print "x: {0:15.10f}".format(self.x)
        print "y: {0:15.10f}".format(self.y)
        print "z: {0:15.10f}".format(self.z)
        print "Polarizability: {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >"\
              .format(self.polar[0], self.polar[1], self.polar[2], self.polar[3])
        print "                {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >"\
              .format(self.polar[4], self.polar[5], self.polar[6], self.polar[7])
        print "                {0:17.9f}".format(self.polar[8])
        

# description of capping: indices of atoms and fragments, excess charge 
class Capping:
    def __init__(self):
        self.oldindex = 0       # old index of the bonding atom
        self.bondindex = (0, 0) # new indices of the bond between fragments       
        self.bond = False       # whether a middlepoint of the bond between fragments is present
        self.charge = 0.0       # excess charge

# -------------------------------------------------------------------------------------------------------
# Class containing all data describing EFP fragment and methods
# used to read it from the file and remove some atoms
class EFPdata:
    def __init__(self):
        self.name = "None"
        self.natoms = 0
        self.nbonds = 0
        self.ncentroids = 0
        self.atoms = []
        self.bonds = []
        self.centroids = []
        self.valel = []
        self.cappings = []

    def FillName(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>1):
                if (s.strip()[0:1]=="$"):
                    break
            s = fp.readline()
        st = s.strip()
        self.name = st[1:len(st)]

    def FillCoords(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>18):
                if (s.strip()[0:18]=="COORDINATES (BOHR)"):
                    break
            s = fp.readline()
        s = fp.readline()
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            if (s[0:1]=="A"):
                self.natoms += 1
                self.atoms.append(EFPatom())
                self.atoms[-1].FillFromString(s)
            if (s[0:1]=="B"):
                self.nbonds += 1
                self.bonds.append(EFPbond())
                self.bonds[-1].FillFromString(s)
            s = fp.readline()  

    def FillCharges(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>9):
                if (s.strip()[0:9]=="MONOPOLES"):
                    break
            s = fp.readline()
        s = fp.readline()
        idxa = 0
        idxb = 0
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            if (s[0:1]=="A"):
                self.atoms[idxa].AddCharge(s)
                idxa += 1
            if (s[0:1]=="B"):
                self.bonds[idxb].AddCharge(s)
                idxb += 1
            s = fp.readline()  

    def FillDipoles(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>7):
                if (s.strip()[0:7]=="DIPOLES"):
                    break
            s = fp.readline()
        s = fp.readline()
        idxa = 0
        idxb = 0
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            if (s[0:1]=="A"):
                self.atoms[idxa].AddDipole(s)
                idxa += 1
            if (s[0:1]=="B"):
                self.bonds[idxb].AddDipole(s)
                idxb += 1
            s = fp.readline()  

    def FillQuadrupoles(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>11):
                if (s.strip()[0:11]=="QUADRUPOLES"):
                    break
            s = fp.readline()
        s = fp.readline()
        idxa = 0
        idxb = 0
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            s2 = fp.readline()  
            if (s[0:1]=="A"):
                self.atoms[idxa].AddQuadrupole(s,s2)
                idxa += 1
            if (s[0:1]=="B"):
                self.bonds[idxb].AddQuadrupole(s,s2)
                idxb += 1
            s = fp.readline()  

    def FillOctupoles(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>9):
                if (s.strip()[0:9]=="OCTUPOLES"):
                    break
            s = fp.readline()
        s = fp.readline()
        idxa = 0
        idxb = 0
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            s2 = fp.readline()  
            s3 = fp.readline()  
            if (s[0:1]=="A"):
                self.atoms[idxa].AddOctupole(s,s2,s3)
                idxa += 1
            if (s[0:1]=="B"):
                self.bonds[idxb].AddOctupole(s,s2,s3)
                idxb += 1
            s = fp.readline()  

    def FillScreen2(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>7):
                if (s.strip()[0:7]=="SCREEN2"):
                    break
            s = fp.readline()
        s = fp.readline()
        idxa = 0
        idxb = 0
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            s = s.strip()
            if (s[0:1]=="A"):
                self.atoms[idxa].AddScreen2(s)
                idxa += 1
            if (s[0:1]=="B"):
                self.bonds[idxb].AddScreen2(s)
                idxb += 1
            s = fp.readline()  

    def FillScreen3(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>7):
                if (s.strip()[0:7]=="SCREEN3"):
                    break
            s = fp.readline()
        s = fp.readline()
        idxa = 0
        idxb = 0
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
            s = s.strip()
            if (s[0:1]=="A"):
                self.atoms[idxa].AddScreen3(s)
                idxa += 1
            if (s[0:1]=="B"):
                self.bonds[idxb].AddScreen3(s)
                idxb += 1
            s = fp.readline() 

    def FillCentroids(self, fp):
        s = fp.readline()
        while s:
            if (len(s)>18):
                if (s.strip()[0:18]=="POLARIZABLE POINTS"):
                    break
            s = fp.readline()
        s = fp.readline()
        while s:
            if (len(s)>4):
                if (s.strip()[0:4]=="STOP"):
                    break
                self.ncentroids += 1
                self.centroids.append(EFPcentroid())
                self.centroids[-1].FillFromString(s)
                s1 = fp.readline()
                s2 = fp.readline()
                s3 = fp.readline()
                self.centroids[-1].AddPolar(s1,s2,s3)
            s = fp.readline()  


    def SetName(self, NewName):
        self.name = NewName

    def ReadFromFile(self, fname):
        fp = open(fname, "r")
        self.FillName(fp)
        self.FillCoords(fp)
        self.FillCharges(fp)
        self.FillDipoles(fp)
        self.FillQuadrupoles(fp)
        self.FillOctupoles(fp)
        self.FillCentroids(fp)
        self.FillScreen3(fp)
        self.FillScreen2(fp)
        fp.close()

    # Data handlers
    def FindClosestAtom(self,x,y,z):
        d0 = 1.0E+15
        ndx = -1
        i = 0
        for a in self.atoms:
            dx = (a.x-x)
            dy = (a.y-y)
            dz = (a.z-z)
            d = dx*dx + dy*dy + dz*dz
            if (d<d0):
                ndx = i   
                d0 = d
            i += 1
        return ndx, d0

    def FindClosestBond(self,x,y,z):
        d0 = 1.0E+15
        ndx = -1
        i = 0
        for b in self.bonds:
            dx = (b.x-x)
            dy = (b.y-y)
            dz = (b.z-z)
            d = dx*dx + dy*dy + dz*dz
            if (d<d0):
                ndx = i   
                d0 = d
            i += 1
        return ndx, d0

    def Bonded(self, i):
        res = []
        for b in self.bonds:
            if b.idx1 == i:
                res.append(b.idx2)
            if b.idx2 == i:
                res.append(b.idx1)
        return res

    def AssignValEl(self):
        self.valel = []
        for c in self.centroids:
            ai, ad = self.FindClosestAtom(c.x, c.y, c.z)
            bi, bd = self.FindClosestBond(c.x, c.y, c.z)
            if (bd <= ad):
                self.valel.append( ( self.bonds[bi].idx1, self.bonds[bi].idx2 ) )
            else:
                self.valel.append( ( ai+1, ai+1 ) )

    def GetTotalCharge(self):
        q = 0.0
        for a in self.atoms:
            q += a.charge + a.num
        for b in self.bonds:
            q += b.charge
        return q

    def GetNumberOfPoints(self):
        return self.natoms + self.nbonds + self.ncentroids

    # adding excess charge to the capping with bonded atom i
    def AddCapCharge(self, i, q):
        for c in self.cappings:
            if c.oldindex == i:
                c.charge += q
                break

    # Deleting atoms with all additional information
    # atoms from the list1 are deleted with all centroids and bonds
    # the bonds and centroids connecting atom from the list2 and atom not from the list1/list2 are not deleted
    # indices in list1 and list2 start from 0
    def DeleteAtoms(self, list1, list2):
        dellist = list1[:]
        dellist.extend(list2)
        dellist.sort(reverse=True)
        newlist = []
        j = 0
        for i in xrange(self.natoms):
            if (i not in dellist):
                newlist.append(j)
                j += 1
            else:
                newlist.append(-1)
        neighbors = {}
        # --------- assigning cappings ----------
        for i in xrange(self.natoms):
            if (i in dellist) and (self.atoms[i].num>1):
                self.cappings.append(Capping())
                neighbors[i+1] = [i+1]
                self.cappings[-1].oldindex = i+1
                for b in self.bonds:
                    if (b.idx1 == i+1):
                        neighbors[i+1].append(b.idx2)
                    elif (b.idx2 == i+1):
                        neighbors[i+1].append(b.idx1)           
        # ---------------------------------------
        for i in dellist:
            if (self.atoms[i].num>1):
#                print "+", self.atoms[i].charge
                self.AddCapCharge(i+1, self.atoms[i].charge+self.atoms[i].num)
            else:
                for j in neighbors:
                    if i+1 in neighbors[j]:
#                        print "+", self.atoms[i].charge
                        self.AddCapCharge(j, self.atoms[i].charge+self.atoms[i].num)
                        break
            self.natoms += -1
 #           print "-", self.atoms[i].charge
            self.atoms.pop(i)
        for b in self.bonds[:]:
            if (b.idx1-1 not in dellist) and (b.idx2-1 not in dellist):
                b.idx1 = newlist[b.idx1-1] + 1
                b.idx2 = newlist[b.idx2-1] + 1
            elif (b.idx1-1 in list1) or (b.idx2-1 in list1):
                if (b.idx1 in neighbors.keys()):
#                   print "+", b.charge
                    self.AddCapCharge(b.idx1, b.charge)
                elif (b.idx2 in neighbors.keys()):
#                    print "+", b.charge
                    self.AddCapCharge(b.idx2, b.charge)
 #               print "-", b.charge
                self.bonds.remove(b)
                self.nbonds += -1
            elif (b.idx1-1 in dellist) and (b.idx2-1 in dellist):
                if (b.idx1 in neighbors.keys()):
#                    print "+", b.charge
                    self.AddCapCharge(b.idx1, b.charge)
                elif (b.idx2 in neighbors.keys()):
#                    print "+", b.charge
                    self.AddCapCharge(b.idx2, b.charge)
#                print "-", b.charge
                self.bonds.remove(b)
                self.nbonds += -1
            else:
                for c in self.cappings:
                    if (c.oldindex == b.idx1) or (c.oldindex == b.idx2):
                        c.bond = True
                        c.bondindex = (newlist[b.idx1-1] + 1, newlist[b.idx2-1] + 1)   
                        break                 
                b.idx1 = newlist[b.idx1-1] + 1
                b.idx2 = newlist[b.idx2-1] + 1
        i = 0
        valelDel = []
        # Delete all centroids between different fragments!!
        for c in self.centroids[:]:
            if (self.valel[i][0]-1 in dellist):
                 valelDel.append(i)
#            if (self.valel[i][0]-1 not in dellist) and (self.valel[i][1]-1 not in dellist):
#                self.valel[i] = (newlist[self.valel[i][0]-1]+1, newlist[self.valel[i][1]-1]+1)
#            elif (self.valel[i][0]-1 in list1) or (self.valel[i][1]-1 in list1):
#                valelDel.append(i)
#            elif (self.valel[i][0]-1 in dellist) and (self.valel[i][1]-1 in dellist):
#                valelDel.append(i)
#            else:
#                self.valel[i] = (newlist[self.valel[i][0]-1]+1, newlist[self.valel[i][1]-1]+1)
            i += 1
        valelDel.sort(reverse=True)
        for i in valelDel:
            self.centroids.pop(i)
            self.ncentroids += -1

    # Output methods
    def Print(self):
        print "\nFragment name: ", self.name
        print "\n ============== ATOMS ==============="
        i = 0
        for a in self.atoms:
            i += 1
            print "\n    ------ atom ", i, "-------"
            a.Print() 
        print "\n ============== BONDS ==============="
        i = 0
        for b in self.bonds:
            i += 1
            print "\n    ------ bond ", i, "-------"
            print "\n",
            b.Print() 
        print "\n ============== CENTROIDS ==============="
        i = 0
        for c in self.centroids:
            i += 1
            print "\n    ------ centroid ", i, "-------"
            print "\n",
            c.Print() 
        print self.valel
        print "\n ============ CAPPINGS ============="
        for c in self.cappings:
            print c.bondindex, c.charge

    def WriteAll2EFP(self, fp):
        fp.write(" $"+self.name+"\n")
        fp.write(" Automatically generated electrostatics and polarization\n")
        fp.write(" COORDINATES (BOHR)\n")
        i = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:15.10f}{3:15.10f}{4:15.10f}{5:12.7f}{6:5.1f}\n".format(\
                     i+1, a.name.ljust(4)[0:4], a.x, a.y, a.z, a.mass, a.num))
            i += 1
        for b in self.bonds:
            fp.write(" BO{0:d}{1:<s}{2:15.10f}{3:15.10f}{4:15.10f}{5:12.7f}{6:5.1f}\n".format(\
                     b.idx1, str(b.idx2).ljust(5-len(str(b.idx1))), b.x, b.y, b.z, 0.0, 0.0))
            
        fp.write(" STOP\n")
        fp.write(" MONOPOLES\n")
        i = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:15.10f}{3:10.5f}\n".format(\
                     i+1, a.name.ljust(4)[0:4], a.charge, a.num))
            i += 1
        for b in self.bonds:
            fp.write(" BO{0:d}{1:<s}{2:15.10f}{3:10.5f}\n".format(\
                     b.idx1, str(b.idx2).ljust(5-len(str(b.idx1))), b.charge, 0.0))
        fp.write(" STOP\n")
        fp.write(" DIPOLES\n")
        i = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:16.10f}{3:16.10f}{4:16.10f}\n".format(\
                     i+1, a.name.ljust(4)[0:4], a.dipole[0], a.dipole[1], a.dipole[2]))
            i += 1
        for b in self.bonds:
            fp.write(" BO{0:d}{1:<s}{2:16.10f}{3:16.10f}{4:16.10f}\n".format(\
                     b.idx1, str(b.idx2).ljust(5-len(str(b.idx1))), b.dipole[0], b.dipole[1], b.dipole[2]))
        fp.write(" STOP\n")
        fp.write(" QUADRUPOLES\n")
        i = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:16.10f}{3:16.10f}{4:16.10f}{5:16.10f} >\n".format(\
                     i+1, a.name.ljust(4)[0:4], a.quadrupole[0], a.quadrupole[1], a.quadrupole[2], a.quadrupole[3]))
            fp.write("        {0:16.10f}{1:16.10f}\n".format(a.quadrupole[4], a.quadrupole[5]))
            i += 1
        for b in self.bonds:
            fp.write(" BO{0:d}{1:<s}{2:16.10f}{3:16.10f}{4:16.10f}{5:16.10f} >\n".format(\
                     b.idx1, str(b.idx2).ljust(5-len(str(b.idx1))), b.quadrupole[0], b.quadrupole[1], b.quadrupole[2], b.quadrupole[3]))
            fp.write("        {0:16.10f}{1:16.10f}\n".format(b.quadrupole[4], b.quadrupole[5]))
        fp.write(" STOP\n")
        fp.write(" OCTUPOLES\n")
        i = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:17.9f}{3:17.9f}{4:17.9f}{5:17.9f} >\n".format(\
                     i+1, a.name.ljust(4)[0:4], a.octupole[0], a.octupole[1], a.octupole[2], a.octupole[3]))
            fp.write("        {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >\n".format(\
                     a.octupole[4], a.octupole[5], a.octupole[6], a.octupole[7]))
            fp.write("        {0:17.9f}{1:17.9f}\n".format(a.octupole[8], a.octupole[9]))
            i += 1
        for b in self.bonds:
            fp.write(" BO{0:d}{1:<s}{2:17.9f}{3:17.9f}{4:17.9f}{5:17.9f} >\n".format(\
                     b.idx1, str(b.idx2).ljust(5-len(str(b.idx1))), b.octupole[0], b.octupole[1], b.octupole[2], b.octupole[3]))
            fp.write("        {0:17.9f}{1:17.9f}{2:17.9f}{3:17.9f} >\n".format(\
                     b.octupole[4], b.octupole[5], b.octupole[6], b.octupole[7]))
            fp.write("        {0:17.9f}{1:17.9f}\n".format(b.octupole[8], b.octupole[9]))
        fp.write(" STOP\n")
        fp.write(" POLARIZABLE POINTS\n")
        i = 0
        for c in self.centroids:
            fp.write(" CT{0:s}{1:16.10f}{2:16.10f}{3:16.10f}\n".format(str(i+1).ljust(5), c.x, c.y, c.z))
            fp.write("        {0:16.10f}{1:16.10f}{2:16.10f}{3:16.10f} >\n".format(\
                    c.polar[0], c.polar[1], c.polar[2], c.polar[3]))
            fp.write("        {0:16.10f}{1:16.10f}{2:16.10f}{3:16.10f} >\n".format(\
                    c.polar[4], c.polar[5], c.polar[6], c.polar[7]))
            fp.write("        {0:16.10f}\n".format(c.polar[8]))
            i += 1
        fp.write(" STOP\n")
        fp.write(" SCREEN2\n")
        i = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:16.10f}{3:14.10f}\n".format(\
                     i+1, a.name.ljust(4)[0:4], a.screen2[0], a.screen2[1]))
            i += 1
        for b in self.bonds:
            fp.write(" BO{0:d}{1:<s}{2:16.10f}{3:14.10f}\n".format(\
                     b.idx1, str(b.idx2).ljust(5-len(str(b.idx1))), b.screen2[0], b.screen2[1]))
        fp.write(" STOP\n")
        fp.write(" $END\n")  


    def WriteCoords2EFP(self, fp):
        fp.write(" FRAGNAME="+self.name+"\n")
        nat = 0
        for a in self.atoms:
            fp.write(" A{0:02d}{1:4s}{2:18.10f}{3:18.10f}{4:18.10f}\n".format(nat+1, a.name.ljust(4)[0:4], a.x*Bohr, a.y*Bohr, a.z*Bohr))
            nat += 1
            if (nat>=3):
               break

    def WriteCoords2EFP_psi4(self, fp):
        fp.write("--\n")
        fp.write(" scf "+self.name+"\n")
        nat = 0
        for a in self.atoms:
            fp.write(" {0:18.10f}{1:18.10f}{2:18.10f}\n".format(a.x*Bohr, a.y*Bohr, a.z*Bohr))
            nat += 1
            if (nat>=3):
               break

    def WriteAtoms2PDB(self, fp, start, resid):
        i = 0
        for a in self.atoms:
            fp.write(" {0:<5d} {1:4s} LIG {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f}{6:6.2f}{7:6.2f}          {8:>2s}  \n".format(\
                     i+start, str(Symb[int(a.num)-1]).ljust(4)[0:4], resid, a.x*Bohr, a.y*Bohr, a.z*Bohr, 0.0, 0.0, Symb[int(a.num)-1]))
            i += 1

    def WriteAll2PDB(self, fp, start, resid):
        i = 0
        for a in self.atoms:
            fp.write("ATOM  {0:<5d} {1:4s} LIG {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f}{6:6.2f}{7:6.2f}          {8:>2s}  \n".format(\
                     i+start, str(Symb[int(a.num)-1]).ljust(4)[0:4], resid, a.x*Bohr, a.y*Bohr, a.z*Bohr, 0.0, 0.0, Symb[int(a.num)-1]))
            i += 1
        for b in self.bonds:
            fp.write("ATOM  {0:<5d} {1:4s} LIG {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f}{6:6.2f}{7:6.2f}          {8:>2s}  \n".format(\
                     i+start, "He".ljust(4)[0:4], resid, b.x*Bohr, b.y*Bohr, b.z*Bohr, 0.0, 0.0, "He"))
            i += 1
        for c in self.centroids:
            fp.write("ATOM  {0:<5d} {1:4s} LIG {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f}{6:6.2f}{7:6.2f}          {8:>2s}  \n".format(\
                     i+start, "Xe".ljust(4)[0:4], resid, c.x*Bohr, c.y*Bohr, c.z*Bohr, 0.0, 0.0, "Xe"))
            i += 1

