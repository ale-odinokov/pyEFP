import sys
import argparse
import math
import openbabel
import pybel
import pyEFP
import fragrepr as fr

parser = argparse.ArgumentParser(description="Constructs cluster for the advanced EFP electrostatics calculation")
parser.add_argument("infile", help="file with input coordinates in pdb format")
parser.add_argument("--cuttype", choices=["mol", "frag", "none"], default="none", help="method of box stacking and cutting, mol is for keeping molecules whole, frag is for keeping fragments whole, none is for keeping system as is, just defying a new center")
parser.add_argument("--rcut", type=float, default=10.0, help="radius for cutting, in Angstroms")
parser.add_argument("--qm", nargs="*", type=int, default=[], help="a list of molecules in the QM subsystem")
parser.add_argument("--dbdir", default="efp_db", help="name of the directory with fragment database")

args = parser.parse_args()

if (args.cuttype=="frag"):
    dbase = fr.fragDB()
    dbase.ReadFromFile(args.dbdir)
    print "\n Reading fragment database as described in the file "+args.dbdir+"/fragments\n"

#print args.infile
#print args.cuttype
#print args.rcut
#print args.qm

openbabel.obErrorLog.SetOutputLevel(0);
mols = pybel.readfile("pdb", args.infile)

mols_sep = []
for mol in mols:
    newm = mol.OBMol.Separate()
    for m in newm:
        mols_sep.append(pybel.Molecule(m))
print "\n", len(mols_sep), "molecules found\n"
openbabel.obErrorLog.SetOutputLevel(1);

unitcell = openbabel.toUnitCell(mol.OBMol.GetData(openbabel.UnitCell))
#print unitcell.GetA()
#print unitcell.GetB()
#print unitcell.GetC()
#print unitcell.GetAlpha()
#print unitcell.GetBeta()
#print unitcell.GetGamma()

a = unitcell.GetA()
b = unitcell.GetB()
c = unitcell.GetC()
cos_a = math.cos(unitcell.GetAlpha()*math.pi/180)
cos_b = math.cos(unitcell.GetBeta()*math.pi/180)
sin_b = math.sin(unitcell.GetBeta()*math.pi/180)
cos_g = math.cos(unitcell.GetGamma()*math.pi/180)
sin_g = math.sin(unitcell.GetGamma()*math.pi/180)
cos_A = (cos_a-cos_b*cos_g)/(sin_b*sin_g)
sin_A = math.sqrt(1-cos_A*cos_A)
cos_th = (cos_a - cos_g*cos_b)/sin_g

# definig center-of-mass of the QM sybsystem,
# if no QM subsystem is used, then define COM of the total system
if (len(args.qm)>0):
    ndxlist = args.qm[:]
else:
    ndxlist = range(1,len(mols_sep)+1)
x = 0.0
y = 0.0
z = 0.0
m = 0.0
for i in ndxlist:
    mol = mols_sep[i-1]
    for at in openbabel.OBMolAtomIter(mol.OBMol):
        x += at.x()*at.GetAtomicMass()
        y += at.y()*at.GetAtomicMass()
        z += at.z()*at.GetAtomicMass()
        m += at.GetAtomicMass()
if (m>0.0):
    x = x/m
    y = y/m
    z = z/m
com = [x, y, z]

# calculates position of the periodic image
def ImageShift(na, nb, nc):
    x = na*a + nb*b*cos_g + nc*c*cos_b
    y = nb*b*sin_g + nc*c*cos_th
    z = nc*c*sin_b*sin_A
    return [x, y, z]

# constructs new box around the origin
for mol in mols_sep:
    x = 0.0
    y = 0.0
    z = 0.0
    m = 0.0
    for at in openbabel.OBMolAtomIter(mol.OBMol):
        x += at.x()*at.GetAtomicMass()
        y += at.y()*at.GetAtomicMass()
        z += at.z()*at.GetAtomicMass()
        m += at.GetAtomicMass()
    if (m>0.0):
        x = x/m
        y = y/m
        z = z/m
    z_c = (z-com[2])/sin_b/sin_A
    y_b = (y-com[1]-z_c*cos_th)/sin_g
    x_a = x - com[0] - y_b*cos_g - z_c*cos_b
    na = x_a/a
    nb = y_b/b
    nc = z_c/c
    da = 0
    if (na >= 0.0):
        while (na-da*1.0>=0.5):
            da += 1
    else:
        while (na-da*1.0<-0.5):
            da -= 1
    db = 0
    if (nb >= 0.0):
        while (nb-db*1.0>=0.5):
            db += 1
    else:
        while (nb-db*1.0<-0.5):
            db -= 1
    dc = 0
    if (nc >= 0.0):
        while (nc-dc*1.0>=0.5):
            dc += 1
    else:
        while (nc-dc*1.0<-0.5):
            dc -= 1

    shift = ImageShift(da, db, dc)
    for at in openbabel.OBMolAtomIter(mol.OBMol):
        x = at.x() - com[0] - shift[0]
        y = at.y() - com[1] - shift[1]
        z = at.z() - com[2] - shift[2]
        at.SetVector(x, y, z)

if (args.cuttype=="none"):
    fpdb = open("cluster.pdb", "w")
    fpdb.write("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f} P 1           1\n".format(a, b, c, unitcell.GetAlpha(), unitcell.GetBeta(), unitcell.GetGamma()))
    nat = 0
    nmol = 0
    for mol in mols_sep:
        pyEFP.WriteOBMol2PDB(fpdb, mol.OBMol, nat+1, nmol+1)
        nat += mol.OBMol.NumAtoms()
        nmol += 1
    fpdb.close()

elif (args.cuttype=="frag"):
    clust = []
    QMcoords = []
    for i in args.qm:
        QMcoords.extend( [atom.coords for atom in mols_sep[i-1]] )
        clust.append(mols_sep[i-1].OBMol)
#    print len(QMcoords), len(QMcoords[0])
    nmax = int(2+args.rcut/max([a, b, c]))
    rmax = args.rcut*args.rcut
    for i in xrange(1,len(mols_sep)+1):
        mol = mols_sep[i-1]
        if i in args.qm: 
            for nx in xrange(-nmax+1,nmax):
                for ny in xrange(-nmax+1,nmax):
                    for nz in xrange(-nmax+1,nmax):
                        shift = ImageShift(nx, ny, nz)
                        inc = 0
                        crd = [atom.coords for atom in mol]
                        for r1 in crd:
                            inner = 0
                            for r2 in QMcoords:
                                dx = r1[0] + shift[0] - r2[0] 
                                dy = r1[1] + shift[1] - r2[1] 
                                dz = r1[2] + shift[2] - r2[2] 
                                dd = dx*dx + dy*dy + dz*dz
                                if (dd <= rmax):
                                    inner = 1
                                    break
                            if (inner):
                                inc = 1
                                break
                        if (inc):
                            if (nx==0) and (ny==0) and (nz==0):
                                pass
                            else: 
                                newmol = openbabel.OBMol(mol.OBMol)
                                newmol.Translate(openbabel.vector3(*shift))
                                clust.append(newmol)
        else: 
            un1 = fr.molunion()
            un1.construct(dbase,mol)
            fcom = []
            fcrd = []
            for f in un1.frags:
                x = 0.0
                y = 0.0
                z = 0.0
                m = 0.0
                fcrd.append([])
                for j in f.refndx:
                    at = mol.OBMol.GetAtom(j)
                    dm = at.GetAtomicMass()
                    x += at.x()*dm
                    y += at.y()*dm
                    z += at.z()*dm
                    fcrd[-1].append(j)
                    m += dm
                fcom.append([x/m, y/m, z/m])
            for nx in xrange(-nmax+1,nmax):
                for ny in xrange(-nmax+1,nmax):
                    for nz in xrange(-nmax+1,nmax):
                        shift = ImageShift(nx, ny, nz)
                        newmol = openbabel.OBMol()
                        for j in xrange(len(fcrd)):
                            inc = 0
                            for k in fcrd[j]:
                                at = mol.OBMol.GetAtom(k)
                                inc2 = 0
                                ax = at.x()
                                ay = at.y()
                                az = at.z()
                                for r2 in QMcoords:
                                    dx = ax + shift[0] - r2[0] 
                                    dy = ay + shift[1] - r2[1] 
                                    dz = az + shift[2] - r2[2] 
                                    dd = dx*dx + dy*dy + dz*dz
                                    if (dd <= rmax):
                                        inc2 = 1
                                        break
                                if (inc2==1):
                                    inc = 1
                                    break
                            if (inc):
                                for k in fcrd[j]:
                                    newmol.AddAtom(mol.OBMol.GetAtom(k))
#                        print "================================================================"
#                        print "                  molecule ", i, " has ", newmol.NumAtoms(), " atoms"
#                        print "================================================================"
                        if (newmol.NumAtoms()>0):
                            newmol.ConnectTheDots()
                            newmol.PerceiveBondOrders()
#                            newmol.AddHydrogens() 
#                            newmol.Translate(openbabel.vector3(*shift))
                            for at in openbabel.OBMolAtomIter(newmol):  # for some reason, Translate does not work
                                at.SetVector(at.x()+shift[0], at.y()+shift[1], at.z()+shift[2])
                            clust.append(newmol)
    fpdb = open("cluster.pdb", "w")
#    fpdb.write("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f} P 1           1\n".format(a, b, c, unitcell.GetAlpha(), unitcell.GetBeta(), unitcell.GetGamma()))
    nat = 0
    nmol = 0
    for mol in clust:
#        print "printing mol ", nmol+1, " with ", mol.NumAtoms(), " atoms"
#        sys.stdout.flush()
#        for at in openbabel.OBMolAtomIter(mol):
#            print at.GetAtomicNum(), at.x(), at.y(), at.z()
#            print at.GetAtomicNum()
#            sys.stdout.flush()
        pyEFP.WriteOBMol2PDB(fpdb, mol, nat+1, nmol+1)
        nat += mol.NumAtoms()
        nmol += 1
    fpdb.close()

elif (args.cuttype=="mol"):
    clust = []
    QMcoords = []
    for i in args.qm:
        QMcoords.extend( [atom.coords for atom in mols_sep[i-1]] )
        clust.append(mols_sep[i-1].OBMol)
    print len(QMcoords), len(QMcoords[0])
    nmax = int(2+args.rcut/max([a, b, c]))
    rmax = args.rcut*args.rcut
    for nx in xrange(-nmax+1,nmax):
        for ny in xrange(-nmax+1,nmax):
            for nz in xrange(-nmax+1,nmax):
                shift = ImageShift(nx, ny, nz)
                for i in xrange(1,len(mols_sep)+1):
                    mol = mols_sep[i-1]
                    inc = 0
                    crd = [atom.coords for atom in mol]
                    for r1 in crd:
                        inner = 0
                        for r2 in QMcoords:
                            dx = r1[0] + shift[0] - r2[0] 
                            dy = r1[1] + shift[1] - r2[1] 
                            dz = r1[2] + shift[2] - r2[2] 
                            dd = dx*dx + dy*dy + dz*dz
                            if (dd <= rmax):
                                inner = 1
                                break
                        if (inner):
                            inc = 1
                            break
                    if (inc):
                        if (i in args.qm) and (nx==0) and (ny==0) and (nz==0):
                            pass
                        else: 
                            newmol = openbabel.OBMol(mol.OBMol)
                            newmol.Translate(openbabel.vector3(*shift))
                            clust.append(newmol)

    fpdb = open("cluster.pdb", "w")
#    fpdb.write("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f} P 1           1\n".format(a, b, c, unitcell.GetAlpha(), unitcell.GetBeta(), unitcell.GetGamma()))
    nat = 0
    nmol = 0
    for mol in clust:
        pyEFP.WriteOBMol2PDB(fpdb, mol, nat+1, nmol+1)
        nat += mol.NumAtoms()
        nmol += 1
    fpdb.close()
                    
      
    


