# This script helps to create input files for the QC calculations of the new fragment
# Only one fragment is treated, but multiple connection patterns are considered 

import sys
import copy
import os
import shutil
import argparse
import openbabel
import pybel
import fragrepr as fr
import pyEFP

parser = argparse.ArgumentParser(description="Prepares input files for the EFP calculation and adds new fragment to the database")
parser.add_argument("infile", help="file describing new fragment using SMARTS")
parser.add_argument("--task", choices=["opt", "efp", "merge"], default="opt", help="defines what to do: prepare files for geometry optimization (opt), efp calculation (efp) or merge data with the existing database")
parser.add_argument("--dbdir", default="efp_db", help="name of the directory with fragment database")

args = parser.parse_args()

# parses QC output and returns geometry of the molecule
def ExtractGeometry(fname):
    try:
        fp = open(fname, "r")
        s = fp.readline()
        found = 0
        while s:
            if (s.strip() == "***** EQUILIBRIUM GEOMETRY LOCATED *****"):
                found = 1
                break
            s = fp.readline()
        if (found == 0):
            print "ERROR: cannot find equilibrium geometry in the file ", fname
            return 1
        s = fp.readline()
        s = fp.readline()
        s = fp.readline()
        s = fp.readline()
        geom = ""
        while s:
            if (len(s.strip())==0):
                break
            else:   
                geom  += s
            s = fp.readline()
        fp.close()
        return geom
    except IOError:
        print "ERROR: cannot open file ", fname
        return 1

# reading input file with the desription of the desired new fragment
fname = args.infile
print "Reading description of the new fragment from the file ", fname
fp = open(fname, "r")
fragname = fp.readline().strip()
corestring = fp.readline().strip()
fragstring = []
s = fp.readline()
while s:
    if ( len(s.strip())>0 ): 
        fragstring.append(s.strip())
    s = fp.readline()
fp.close()

core = pybel.readstring("smi",corestring)
coreSMARTS = pybel.Smarts(corestring)

# checking for self-consistency
for f in fragstring:
    fmol = pybel.readstring("smi", f)
    res = coreSMARTS.findall(fmol)
    if ( len(res) == 0 ):
        sys.exit("ERROR: cannot find core "+corestring+" in the structure "+f)
    for i in xrange(1,fmol.OBMol.NumAtoms()+1):
        atm = fmol.OBMol.GetAtom(i)
        if (i not in res[0]) and ( atm.IsHydrogen() == False):
            bonded = 0
            for at in openbabel.OBAtomAtomIter(fmol.OBMol.GetAtom(i)):
                if at.GetIndex()+1 in res[0]:
                    bonded = 1
                    break
            if (bonded == 0):
                sys.exit("ERROR: an atom in the structure "+f+" is not directly connected to the core of the fragment")

# removing redundant bond patterns
redundant = []
for i in xrange(len(fragstring)):
    f1 = fragstring[i]
    f1mol = pybel.readstring("smi", f1)
    f1SMARTS = pybel.Smarts(f1)
    for j in xrange(i+1, len(fragstring)):
        f2 = fragstring[j]
        f2mol = pybel.readstring("smi", f2)
        f2SMARTS = pybel.Smarts(f2)
        res1 = f1SMARTS.findall(f2mol)
        res2 = f2SMARTS.findall(f1mol)
        if (len(res1)>0) and (len(res2)>0):
            redundant.append(f1)
            break

ndel = 0
for rfr in redundant:
    fragstring.remove(rfr)
    ndel += 1
print "Found", ndel, "redundant bond patterns, they will be removed"

# generating all needed bond patterns, including those obtained after the bond breaking
finalstring = [corestring]
mols = []
for f in fragstring[:]:
    fmol = pybel.readstring("smi", f)
    rescore = coreSMARTS.findall(fmol)
    batoms = []
    bonds = []
    for i in xrange(1,fmol.OBMol.NumAtoms()+1):
        atm = fmol.OBMol.GetAtom(i)
        if (i not in rescore[0]) and (atm.IsHydrogen() == False):
            batoms.append(i)
            bonds.append([i, []])
            for ba in openbabel.OBAtomAtomIter(fmol.OBMol.GetAtom(i)):
                bonds[-1][1].append(ba)
    ndx = [[], []]
    for i in batoms:
        ndx[-1].append( set([i]) )
    for i in xrange(len(batoms)-1):     # creating all possible indices for the given bond pattern 
        ndx.append([])
        for s in ndx[-2]:
            for j in batoms:
                newset = copy.deepcopy(s)
                if (j not in newset):
                    newset.add(j)
                    unique = 1
                    for ss in ndx[-1]:
                        if (ss == newset):
                            unique = 0
                            break
                    if (unique==1):
                        ndx[-1].append(newset)
    for n in ndx:          # saving only unique bond patterns
        for nn in n:
            newmol = openbabel.OBMol(fmol.OBMol)
            delatoms = []
            for i in xrange(1,newmol.NumAtoms()+1):
#                print rescore, nn, i
                if (i not in rescore[0]) and (i not in nn):
                    delatoms.append(newmol.GetAtom(i))
            for da in delatoms:
                    newmol.DeleteAtom(da)

            newmol.PerceiveBondOrders()
            newmol.AddHydrogens()   
            newstring = pybel.Molecule(newmol).write("smi")
            newSMARTS = pybel.Smarts(newstring)
            unique = 1
            for ff in finalstring:
                ffSMARTS = pybel.Smarts(ff)
                ffmol = pybel.readstring("smi", ff)
                res1 = ffSMARTS.findall( pybel.Molecule(newmol) )
                res2 = newSMARTS.findall(ffmol)
                if (len(res1)>0) and (len(res2)>0):
                    unique = 0
                    break 
            if (unique==1):
                finalstring.append(newstring.strip())

print "A complete set of ", len(finalstring), "non-redundant bond patterns has been generated"

# reading current state of the fragment database
db = fr.fragDB()
print "\nReading fragment database as described in the file "+args.dbdir+"/fragments\n"
try:
    db.ReadFromFile(args.dbdir+"/fragments")
except Exception:
    pass

ndel = 0
db_data = fr.DBentry()
if (len(db.entries)>0):
    # looking for the new core in the database
    for e in db.entries:
        eSMARTS = pybel.Smarts(e.smarts)
        emol = pybel.readstring("smi",e.smarts)
        res1 = eSMARTS.findall(core)
        res2 = coreSMARTS.findall(emol)
        if (len(res1)>0) and (len(res2)>0):
            print "Found the same core in the database, changing name ", fragname, " to ", e.name
            fragname = e.name
            db_data = e
            break

    for f in finalstring[:]:      #  deleting bond patterns that are already in the database 
        ffSMARTS = pybel.Smarts(f)
        ffmol = pybel.readstring("smi", f)
        for bp in db_data.bondpatterns.values():
            for bps in bp:
                dbmol = pyEFP.EFP2mol(args.dbdir+"/"+fragname+"/"+bps+".efp")
                dbstring = pybel.Molecule(dbmol).write("smi")
                dbSMARTS = pybel.Smarts(dbstring)
                res1 = dbSMARTS.findall(ffmol)
                res2 = ffSMARTS.findall(pybel.Molecule(dbmol))
                print f, dbstring
                if (len(res1)>0) and (len(res2)>0):
                    finalstring.remove(f)
                    print "DELETED!"
                    ndel += 1

print "Removing ", ndel, " bond patterns for the fragment ", fragname, " that are already in the database"

if (args.task == "opt"):
    # Creating input files for quantum-chemical calculations
    header_opt = ""
    fopt = open("templates/header_opt", "r")
    s = fopt.readline()
    while s:
        header_opt += s
        s = fopt.readline()
    fopt.close()

    dirname = "QC_data_"+fragname
    if os.path.exists(dirname):
        sys.exit("ERROR: The directory "+dirname+" already exists")
    else:
        print "Input files for quantum-chemical calculations will be stored in the directory ", dirname
        os.makedirs(dirname)

    print "Relationship between structures and file names:"
    n = 1
    for ff in finalstring:
        fname = dirname+"/"+fragname+"_"+str(n)+"_opt.inp"
        print ff, "\t", fragname+"_"+str(n) 
        ffmol = pybel.readstring("smi", ff)
        ffmol.make3D(forcefield="mmff94", steps=50)
#        ffmol.write(format="pdb", filename=dirname+"/"+db_data.name+"_"+str(n)+".pdb")
        ffmol.write(format="inp", filename=fname)

        content = ""
        fp = open(fname, "r")
        s  = fp.readline()
        s  = fp.readline()
        while s:
            content += s
            s = fp.readline()
        fp.close()

        fp = open(fname, "w")
        fp.write(header_opt.rstrip('\r\n') + '\n' + content)
        fp.close()

        n += 1

elif (args.task == "efp"):
    # Creating input files for quantum-chemical calculations
    header_efp = ""
    fefp = open("templates/header_efp", "r")
    s = fefp.readline()
    while s:
        header_efp += s
        s = fefp.readline()
    fefp.close()

    dirname = "QC_data_"+fragname
    if os.path.exists(dirname):
        print "Input files for quantum-chemical calculations will be stored in the directory ", dirname
    else:
        sys.exit("ERROR: The directory "+dirname+" does not exist")

    print "Relationship between structures and file names:"
    n = 1
    for ff in finalstring:
        fname0 = dirname+"/"+fragname+"_"+str(n)+"_opt.out"
        print ff, "\t", fragname+"_"+str(n) 
        geom = ExtractGeometry(fname0)
        if (geom == 1):
            sys.exit("Error while parsing file "+fname0)
        else:
            fname = dirname+"/"+fragname+"_"+str(n)+"_efp.inp"
            fp = open(fname, "w")
            fp.write(header_efp.rstrip('\r\n') + '\n' + "\n $DATA\n\nC1\n" + geom + " $END\n")
            fp.close()
        n += 1

elif (args.task == "merge"):
    # Merging new data with the existing database
    dirname = "QC_data_"+fragname
    if os.path.isfile(args.dbdir+"/fragments"):
        print "saving file ",  args.dbdir+"/fragments", " as ", args.dbdir+"/#fragments#"
        os.rename(args.dbdir+"/fragments", args.dbdir+"/#fragments#")

    newMetaData = {}   # preparing new metadata
    nat_core = core.OBMol.NumHvyAtoms()
    print "Relationship between structures and file names:"
    n = 1
    for ff in finalstring:
        fname0 = dirname+"/"+fragname+"_"+str(n)+"_opt.out"
        print ff, "\t", fragname+"_"+str(n) 
        ffmol = pybel.readstring("smi", ff)
        nat = ffmol.OBMol.NumHvyAtoms() - nat_core
        if nat not in newMetaData.keys():
            newMetaData[nat] = []
        newMetaData[nat].append(fragname+"_"+str(n))
        n += 1

    fp = open(args.dbdir+"/fragments", "w")  # creating new metafile
    done = 0
    for e in db.entries:
        if (e.name != fragname):
            fp.write("$frag\n")
            fp.write(e.name+"\n")
            fp.write(e.smarts+"\n")
            for b in e.bondpatterns.keys():
                fp.write(str(b))
                for s in e.bondpatterns[b]:
                    fp.write(" "+s)
                fp.write("\n")
            fp.write("$endfrag\n\n")
        else:
            done = 1
            fp.write("$frag\n")
            fp.write(e.name+"\n")
            fp.write(e.smarts+"\n")
            for b in e.bondpatterns.keys():
                fp.write(str(b))
                if (b in newMetaData.keys()):
                    for s in newMetaData[b]:
                        name = s.split("_")[0]
                        n = int(s.split("_")[1])
                        NameList = []
                        for ss in e.bondpatterns.values():
                            NameList.extend(ss)
                        while (name+"_"+str(n) in NameList):    
                            n += 1
#                            print "try: ", name+"_"+str(n)
                        e.bondpatterns[b].append(name+"_"+str(n))
                        shutil.copyfile("QC_data_"+fragname+"/"+s+"_efp.efp", args.dbdir+"/"+fragname+"/"+name+"_"+str(n)+".efp")
#                        print s, name+"_"+str(n)
                for s in e.bondpatterns[b]:
                    fp.write(" "+s)
                fp.write("\n")
            fp.write("$endfrag\n\n")
            break
    if (done == 0):
        os.makedirs(args.dbdir+"/"+fragname)
        fp.write("$frag\n")
        fp.write(fragname+"\n")
        fp.write(corestring+"\n")
        for b in newMetaData.keys():
            fp.write(str(b))
            for s in newMetaData[b]:
                shutil.copyfile("QC_data_"+fragname+"/"+s+"_efp.efp", args.dbdir+"/"+fragname+"/"+s+".efp")
                fp.write(" "+s)
            fp.write("\n")
        fp.write("$endfrag\n\n")
    fp.close()

"""    
    for d in newMetaData.keys():
        print d, 
        for dd in newMetaData[d]:
            print dd,
        print ""
"""



