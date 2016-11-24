import sys
import os
import argparse
import openbabel
import pybel
import fragrepr as fr
import pyEFP


parser = argparse.ArgumentParser(description="Reads coordinates from the .pdb file and prepares GAMESS US input file for pure EFP or QM/EFP calculation")
parser.add_argument("infile", help="file with coordinates in the .pdb format")
parser.add_argument("--dbdir", default="efp_db", help="name of the directory with fragment database")
parser.add_argument("--qm", nargs="*", type=int, default=[], help="a list of molecules in the QM subsystem")
parser.add_argument("--format", choices=["gamess", "psi4"], default="gamess", help="input file format")

args = parser.parse_args()

db = fr.fragDB()
db.ReadFromFile("./"+args.dbdir)
print "\n Reading fragment database as described in the file ./"+args.dbdir+"/fragments\n"
#db.Print()

openbabel.obErrorLog.SetOutputLevel(0);
print  "\n Reading coordinates from the file ", args.infile
mols = list(pybel.readfile("pdb", args.infile))
mols_sep = []
for mol in mols:
    newm = mol.OBMol.Separate()
    for m in newm:
        mols_sep.append(pybel.Molecule(m))
#        print mols_sep[-1].write("smi")
print len(mols_sep), " molecules found\n"
openbabel.obErrorLog.SetOutputLevel(1);

#mols_sep[0].write("pdb", "test_out.pdb")

FR = []
newMol = []
QMmols = []

ndx = 0
for mol in mols_sep:

    ndx += 1
    if (ndx not in args.qm):

        nfr = 0
        FR.append([])
        newMol.append([])

        un1 = fr.molunion()
        if (un1.construct(db,mol)!=0):
            sys.exit("An error occured, aborting run")

#        un1.Print()

        Fitmol = un1.ComputeBondResidual(db)
        Fitmol = un1.FitCoords(db)

        for i in xrange(len(un1.frags)):
            f = un1.frags[i]
            FR[-1].append(pyEFP.EFPdata())
            FR[-1][nfr+i].ReadFromFile(args.dbdir+"/"+f.ename+"/"+f.fname+".efp")
            FR[-1][nfr+i].name = str(ndx)+"_"+str(nfr+i+1)
            FR[-1][nfr+i].AssignValEl()
            list1 = []
            list2 = []
            for j in xrange(len(f.links)):
                for ln in un1.frags[f.links[j][2]].links:
                    if ln[2] == i:
                        res = ln[-1]
                        break
                if f.links[j][-1] < res:
                    list2.append(f.links[j][0]-1)
            for k in xrange(FR[-1][nfr+i].natoms):
                if (k+1 not in f.efpcore) and (k not in list2):
                    if FR[-1][nfr+i].atoms[k].num > 1:
                        list1.append(k)
                    else:
                        if fr.IsIncluded(FR[-1][nfr+i].Bonded(k+1), f.efpcore) == False:
                            list1.append(k)
            FR[-1][nfr+i].DeleteAtoms(list1,list2)

            newMol[-1].append(pyEFP.ApplyKabschEFP(FR[-1][nfr+i], f.rotat[2], f.rotat[1], f.rotat[0]))

        for i in xrange(len(FR[-1])):
            for c in FR[-1][nfr+i].cappings:
                if c.bond:
                    for b in FR[-1][nfr+i].bonds:
                        if (b.idx1 in c.bondindex) and (b.idx2 in c.bondindex):
                            b.charge += c.charge
                            break
                else:
                    for ln in un1.frags[i].links:
                        if ln[0] == c.oldindex:
                            for ln2 in un1.frags[ln[2]].links:
                                if ln2[2] == i:
                                    for c2 in FR[-1][nfr+ln[2]].cappings:
                                        if (c2.oldindex == ln2[0]):
                                            for b in FR[-1][nfr+ln[2]].bonds:
                                                if (b.idx1 in c2.bondindex) and (b.idx2 in c2.bondindex):
                                                    b.charge += c.charge
                                                    break
                                    break
                            break
    else:
        QMmols.append(mol.OBMol)

if (args.format == "gamess"):
    header_efp = ""
    fopt = open("templates/header_QMefp_GS", "r")
    s = fopt.readline()
    while s:
        header_efp += s
        s = fopt.readline()
    fopt.close()

    fp = open("efp_input.inp", "w")

    fp.write(header_efp.rstrip('\r\n'))
    fp.write("\n\n $DATA\n\n")
    fp.write("C1\n")
    for m in QMmols:
        pyEFP.WriteGAMESS_DATA(fp, m)
    fp.write(" $END\n")

    fp.write(" $EFRAG\n")
    fp.write("  iscrpol=1 iscrelec=0 NOEXREP NOCHTR NODISP \n")
    for i in xrange(len(FR)): 
        for j in xrange(len(FR[i])): 
           newMol[i][j].WriteCoords2EFP(fp)
    fp.write(" $END\n")
    for i in xrange(len(FR)): 
        for j in xrange(len(FR[i])): 
            FR[i][j].WriteAll2EFP(fp)
#            FR[i][j].Print()
#    fp.write(" $END\n")

    fp.close()

elif (args.format == "psi4"):
    os.mkdir("psi4_input")
    fp = open("./psi4_input/psi4_input.dat", "w")
    fp.write("molecule {\n")
    for m in QMmols:
        pyEFP.WriteOBMol2xyz(fp, m)
    for i in xrange(len(FR)): 
        for j in xrange(len(FR[i])): 
           newMol[i][j].WriteCoords2EFP_psi4(fp)
    fp.write("}\n")
    fp.write("\nset {\n    basis 6-31g\n}\n")
    fp.write("\nenergy('PBE0')\n")
    fp.close()

    for i in xrange(len(FR)): 
        for j in xrange(len(FR[i])): 
            fp = open("./psi4_input/"+FR[i][j].name+".efp", "w")
            FR[i][j].WriteAll2EFP(fp)
            fp.write("\n")
            fp.close()

fpdb = open("efp.pdb", "w")

nat = 0
nmol = 0
for i in xrange(len(QMmols)):
    pyEFP.WriteOBMol2PDB(fpdb, QMmols[i], nat+1, nmol+1)
    nat += QMmols[i].NumAtoms()
    nmol += 1
for i in xrange(len(FR)): 
    for j in xrange(len(FR[i])): 
        newMol[i][j].WriteAll2PDB(fpdb, nat+1, nmol+1)
        nat += newMol[i][j].GetNumberOfPoints()
    nmol += 1

fpdb.close()





