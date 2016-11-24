# this file contains classes and functions created to read fragment database
# and automatically subdivide molecule into known fragments

import sys
import openbabel
import pybel
import pyEFP
import numpy as np

# Checks if elements of the list1 are contained in the list2
def IsIncluded(list1, list2):
    for l in list1:
        if (l not in list2):
            return False
    return True

# Checks if two lists share at least one element
def ListOverlap(list1, list2):
    for l in list1:
        if l in list2:
            return True
    return False

# contains description of the single entry in the fragment database
class DBentry:
    def __init__(self):
        self.name = ""
        self.smarts = ""
        self.bondpatterns = {}

    def FillFromFile(self,fp):
        s=fp.readline()
        while s:
            s.strip()
            if (len(s)>0):
                if (s[0:5]=="$frag"):
                    s=fp.readline()
                    break
            s=fp.readline()
        self.name = s.strip()
        s=fp.readline()
        self.smarts = s.strip()
        s=fp.readline()
        while s:
            s.strip()
            if (len(s)>0):
                if (len(s)>8):
                    if (s[0:8]=="$endfrag"):
                        return True
                st=s.split()
                num = int(st.pop(0))
                self.bondpatterns[num] = st
            s=fp.readline()
        return False

# contains description of the fragment database
class fragDB:
    def __init__(self):
        self.nentries = 0
        self.entries = []
        self.dirname = ""

    def ReadFromFile(self, dirname):
        self.dirname = dirname
        fp = open(dirname+"/fragments", "r")
        res = True
        while (res == True):
            self.entries.append(DBentry())
            self.nentries += 1
            res = self.entries[-1].FillFromFile(fp)
        fp.close()
        self.nentries += -1
        self.entries.pop(-1)

    def Smarts(self, name):
        for e in self.entries:
            if (e.name == name):
                return e.smarts
        return ""

    def Bondpatterns(self, name, nbonds):
        for e in self.entries:
            if (e.name == name):
                if nbonds in e.bondpatterns.keys():
                    return e.bondpatterns[nbonds]
                else:
                    return []

    def Print(self):
        i = 0
        for f in self.entries:
            i += 1
            print "========== ENTRY ", i, " =========="
            print "name:    ", f.name
            print "smarts:  ", f.smarts
            print "bond patterns:  "
            for b in f.bondpatterns:
                print b, f.bondpatterns[b]  

# instance of the fragment (part of the molecule)
# all indices start from 1
class frag:
    def __init__(self):
        self.ename = ""  # name of the entry in the database
        self.fname = ""  # name of the file with EFP parameters 
        self.molcore = []   # list of indices of the fragment core in the molecule
        self.moltotal = []  # list of indices of the total fragment (with connections) in the molecule
        self.efpcore = []   # list of indices of the fragment core in the EFP
        self.efptotal = []  # list of indices of the total fragment (with connections) in the EFP
        self.subMol = openbabel.OBMol()    # corresponding submolecule as a part of total molecule
        self.links = []     # list of connections with other residues
        self.rotat = []     # parameters of the best fit
        self.refndx = []    # list of indices of all atoms in the fragment according to the initial method of numbering

# molecule as a union of fragments
class molunion:
    def __init__(self):
        self.frags = []

    def construct(self, db, mol):
        cores = []
        support = []
        nfrags = 0
        for e in db.entries:
            fr = pybel.Smarts(e.smarts)
            found = fr.findall(mol)
            if len(found)>0:
                cores.append([e.name, found])
                for f in found:
                    support.append([e.name, f])
                    nfrags += 1

        print "\nTotal number of database matches: ", nfrags, "\n"

        # constructing overlap matrix of the fragments
        print "Constructing overlap matrix"
        overlap = []
        for i in xrange(len(support)):
            overlap.append([])
            for j in xrange(len(support)):
                overlap[-1].append(ListOverlap(support[i][1], support[j][1]))

        # subdivide molecule into non-overlapping fragments without any atoms left
        print "Creating a set of unique fragments"
        numat = mol.OBMol.NumHvyAtoms()
        complete = False
        result = None
        nums = []
        for i in xrange(len(support)):
            nums.append([i])
        n0 = 0
        n = len(nums)
        ntry = 0
        for i in xrange(len(support)-1):
            for nn in nums[n0:n]:
                for j in xrange(len(support)):
                    if j not in nn:
                        good = 1
                        for k in nn:
                            if overlap[j][k]:
                                good = 0
                                break
                        if good:
                            newl = nn[:]
                            newl.append(j)
                            sstop = 0 
                            for nn2 in nums[n:]:
                                if IsIncluded(nn2,newl):
                                    sstop = 1
                                    break
                            if (sstop == 0):
                                ntry += 1
                                nat = 0
                                for k in newl:
                                    nat += len(support[k][1])
                                if (nat == numat):
                                    complete = True
                                    result = newl
                                    break
                                nums.append(newl)
#                                print i, j, newl
                if (complete):
                    break
            if (complete):
                break
            n0 = n
            n = len(nums)

        if (result == None):
            for nn in nums:
                nat = 0
                for k in nn:
                    nat += len(support[k][1])
                if (nat == numat):
                    result = nn
                    break
#            result = [0]

        # convert data to the appropriate format
        coredict = {}
        nfrag = 0
        for ss in result:
            if support[ss][0] not in coredict:
                coredict[support[ss][0]] = [support[ss][1]]
                nfrag += 1
            else:
                coredict[support[ss][0]].append(support[ss][1])
                nfrag += 1
        cores = []
        print "The molecule was subdivided into ", nfrag, " non-overlapping fragments after ", ntry+1, " attempts\n" 
        for c in coredict:
            cores.append([c, coredict[c]])

        for c in cores:
#            print c
            for cc in c[1]:
                self.frags.append(frag())
                self.frags[-1].ename = c[0]
                self.frags[-1].refndx = list(cc[:])
                WithBonds0 = list(cc[:])
                for i in cc:
                    for at in openbabel.OBAtomAtomIter(mol.OBMol.GetAtom(i)):
                        if (at.GetIndex()+1 not in cc):
                            if (at.GetAtomicNum() > 1):
                                WithBonds0.append(at.GetIndex()+1) 
                            else:
                                self.frags[-1].refndx.append(at.GetIndex()+1)
                # Creating submolecule from the fragment
#                subMol = openbabel.OBMol()
                for at in openbabel.OBMolAtomIter(mol.OBMol):
                    if (at.GetIndex()+1) in WithBonds0:
                        at1 = openbabel.OBAtom()
                        at1.Duplicate(at)
                        self.frags[-1].subMol.AddAtom(at1)
                self.frags[-1].subMol.ConnectTheDots()
                self.frags[-1].subMol.PerceiveBondOrders()
                self.frags[-1].subMol.AddHydrogens()    
                # creating new lists for submolecule
                frB = pybel.Smarts(pybel.Molecule(self.frags[-1].subMol).write("smi"))
#                print "========================  fragment", self.frags[-1].ename, " ==========================="
#                print self.frags[-1].ename, pybel.Molecule(self.frags[-1].subMol).write("smi")
                fr = pybel.Smarts(db.Smarts(c[0]))
                self.frags[-1].moltotal = frB.findall(pybel.Molecule(self.frags[-1].subMol))[0]
                newcore = fr.findall(pybel.Molecule(self.frags[-1].subMol))[0]  # indices in the submolecule
                self.frags[-1].molcore = []
                for f in self.frags[-1].moltotal:
                    if f in newcore:
                        self.frags[-1].molcore.append(f)
                # creating new lists for EFP fragment
                nbd = len(self.frags[-1].moltotal) - len(self.frags[-1].molcore)  # number of connections
                bpat = db.Bondpatterns(c[0], nbd)
                found = 0
                for fn in bpat:
                    efpmol = pyEFP.EFP2mol("./"+db.dirname+"/"+str(c[0])+"/"+fn+".efp")
#                    print self.frags[-1].fname, pybel.Molecule(efpmol).write("smi")
                    newlB = frB.findall(pybel.Molecule(efpmol))
#                    print c[0], fn, efpmol.NumAtoms(), len(self.frags[-1].moltotal) - len(self.frags[-1].molcore), newlB
#                    print pybel.Molecule(efpmol).write("smi")
#                    print "frB = ", pybel.Molecule(self.frags[-1].subMol).write("smi")
#                    print frB.findall(mol)
                    if len(newlB) > 0:
                        self.frags[-1].efptotal = newlB[0]
                        found = 1
                        break
                if (found==0):
                    print "Cannot find bond pattern for fragment ", self.frags[-1].ename, " corresponding to the structure ", pybel.Molecule(self.frags[-1].subMol).write("smi").strip(), "( ", nbd, " connections )"
                    return 1
                self.frags[-1].fname = fn
                newl = fr.findall(pybel.Molecule(efpmol))[0]
                self.frags[-1].efpcore = []
                for f in self.frags[-1].efptotal:
                    if f in newl:
                        self.frags[-1].efpcore.append(f)
        # calculating links
        # first loop to identify index in local reference
        for f in self.frags:
            for i in xrange(len(f.efptotal)):
                if (f.efptotal[i] not in f.efpcore):
                    f.links.append([f.efptotal[i], f.moltotal[i]])
        # second loop to identify linked residue and index in the corresponding reference
        ifr0 = 0
        for f in self.frags:
            for xx in xrange(len(f.links)):
                i = f.links[xx][1]
                x0 = f.subMol.GetAtom(i).x()
                y0 = f.subMol.GetAtom(i).y()
                z0 = f.subMol.GetAtom(i).z()
                found = 0
                ifr = 0
                for f2 in self.frags:
#                    if (f2.molcore[0] != f.molcore[0]):
                    if (ifr0 != ifr):
                        for j in f2.molcore:
                            x1 = f2.subMol.GetAtom(j).x()
                            y1 = f2.subMol.GetAtom(j).y()
                            z1 = f2.subMol.GetAtom(j).z()
                            dx = x1-x0
                            dy = y1-y0
                            dz = z1-z0
                            dr2 = dx*dx + dy*dy + dz*dz
                            if dr2 < 0.0001:
                                found = 1
                                f.links[xx].append(ifr)
#                                print xx, len(f.links), f2.moltotal.index(j), len(f2.efptotal), f2.fname, f2.ename
#                                sys.stdout.flush()
                                f.links[xx].append(f2.efptotal[f2.moltotal.index(j)])
                                f.links[xx].append(j)
                                break
                        if (found == 1):
                            break
                    ifr += 1
            ifr0 += 1
        return 0


    def FitCoords(self, db):
        Fitmol = openbabel.OBMol()
        for f in self.frags:
            # constructing submolecule from the fragment
            efpmol = pyEFP.EFP2mol("./"+db.dirname+"/"+f.ename+"/"+f.fname+".efp")
            # fit atoms from the core lists
            sys.stdout.flush()
            rot, COMRef, COMFit = pyEFP.SolveKabsch(efpmol, f.subMol, f.efpcore, f.molcore)
            f.rotat = [rot, COMRef, COMFit]
            # deleting extra atoms
            complete = 0
            while (complete == 0):
                complete = 1
                for at in openbabel.OBMolAtomIter(efpmol):
                    extra = 0
                    if at.GetIndex()+1 not in f.efpcore:
                        extra = 1
                        for at1 in openbabel.OBAtomAtomIter(at):
                            if (at1.GetIndex()+1 in f.efpcore) and (at.GetAtomicNum() == 1):
                                extra = 0
                                break
                    if (extra == 1):
                        efpmol.DeleteAtom(at)
                        complete = 0
            pyEFP.ApplyKabsch(efpmol, COMFit, COMRef, rot)
            Fitmol += efpmol
        return Fitmol

    def ComputeBondResidual(self, db):
        for f in self.frags:
            # constructing submolecule from the fragment
            efpmol = pyEFP.EFP2mol("./"+db.dirname+"/"+f.ename+"/"+f.fname+".efp")
            # fit atoms from the core lists
            rot, COMRef, COMFit = pyEFP.SolveKabsch(efpmol, f.subMol, f.efpcore, f.molcore)
            # deleting extra atoms
            complete = 0
            while (complete == 0):
                complete = 1
                for at in openbabel.OBMolAtomIter(efpmol):
                    extra = 0
                    if at.GetIndex()+1 not in f.efptotal:
                        extra = 1
                        for at1 in openbabel.OBAtomAtomIter(at):
                            if (at1.GetIndex()+1 in f.efptotal) and (at.GetAtomicNum() == 1):
                                extra = 0
                                break
                    if (extra == 1):
                        efpmol.DeleteAtom(at)
                        complete = 0
            efpmol.AddHydrogens()
            pyEFP.ApplyKabsch(efpmol, COMFit, COMRef, rot)
            for i in xrange(len(f.links)):
                dx = efpmol.GetAtom(f.links[i][0]).x() - self.frags[f.links[i][2]].subMol.GetAtom(f.links[i][4]).x() 
                dy = efpmol.GetAtom(f.links[i][0]).y() - self.frags[f.links[i][2]].subMol.GetAtom(f.links[i][4]).y() 
                dz = efpmol.GetAtom(f.links[i][0]).z() - self.frags[f.links[i][2]].subMol.GetAtom(f.links[i][4]).z()
                f.links[i].append(np.sqrt(dx*dx + dy*dy + dz*dz))
           

    def Print(self):
        i = 0
        for f in self.frags:
            i += 1
            print "===========  FRAGMENT ", i, " ============"
            print "database name:    ", f.ename
            print "file fname:    ", f.fname
            print "molcore:  ", f.molcore
            print "moltotal: ", f.moltotal
            print "efpcore:  ", f.efpcore
            print "efptotal: ", f.efptotal
            print "number of atoms in the submolecule:  ", f.subMol.NumAtoms()
            print "connections: "
            for l in f.links:
                print l
            print "reference:   ", f.refndx

