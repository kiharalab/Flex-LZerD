#!/usr/bin/env python2
# -*- coding: UTF-8 -*-


import argparse
import sys
import os
import itertools
import tempfile
import shutil
import math
import subprocess
import functools
from Bio.PDB import *
from Bio.Data.SCOPData import protein_letters_3to1
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
from numpy import *
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
import warnings
import numpy
from scipy.linalg import sqrtm
import scipy
#from Bio.PDB.Vector import rotaxis
#from Bio.PDB.vectors import rotaxis
#from Bio.PDB.vectors import calc_dihedral
from Bio.PDB import rotaxis
from Bio.PDB import calc_dihedral
import Bio

import time
import tarfile

import prody

import prody.dynamics.rtb
from prody.dynamics.rtb import RTB
from prody.dynamics.anm import ANM
from prody.dynamics.gamma import Gamma
from prody.dynamics.gamma import GammaStructureBased
from prody.dynamics.gamma import GammaVariableCutoff

from scipy.optimize import minimize

import shutil

prody.confProDy(auto_secondary=False)

def doPhenixGeometryMinimization(struct, keepLog=False, receptor=None):
	struct = struct.copy()

	lstruct = struct
	lchains = set(lstruct.getChids())

	if receptor is None:
		pass
	else:
		struct = struct + receptor

	t = tempfile.mkdtemp(dir="/dev/shm")
	f = t+"/struct.pdb"
	
	if ("phenix_extra_input" in os.environ) and (os.environ["phenix_extra_input"] != ""):
		extra_f = os.environ["phenix_extra_input"]
		shutil.copy(extra_f, t+"/"+extra_f)
		
	
	prody.writePDB(f, struct)

	try:
		ret = subprocess.check_output('bash ~/hingebin/dophenixgeomincwd_v2 '+os.getcwd(), shell=True, cwd=t)
	except subprocess.CalledProcessError as e:
		print(e.returncode)
		print(e.output)
		exit(1)
		

	modstruct = prody.parsePDB(t+"/struct_minimized.pdb")
	modstruct = modstruct.select("chain %s" % " ".join(lchains))

	shutil.rmtree(t)

	return modstruct

# ARGS

bound_file = sys.argv[1]
unbound_file = sys.argv[2]
#eigenvalue_file = sys.argv[3]
#eigenvector_file = sys.argv[4]
outpref = sys.argv[6]
receptor_file = sys.argv[7]

shouldRenumber = (sys.argv[3] == "renumber")

# use BioPython to set up nice PDB files with matching residues
parser = PDBParser()
bound_struct = parser.get_structure(bound_file, bound_file)
unbound_struct = parser.get_structure(unbound_file, unbound_file)
receptor_struct = parser.get_structure(receptor_file, receptor_file)
print("[info] loaded PDB files")

def get_res_crid(r):
	return (r.get_parent().get_id(), r.get_id())

present = set([get_res_crid(r) for r in bound_struct.get_residues()]) & set([get_res_crid(r) for r in unbound_struct.get_residues()])
present = sorted(present)
if sys.argv[5] != "all":
	num_modes = int(sys.argv[5])
else:
	num_modes = 3*len(list(unbound_struct.get_residues())) - 6

boch = bound_struct[0]
ubch = unbound_struct[0]
for rid in present:
	boch_rid = boch[rid[0]][rid[1]]
	ubch_rid = ubch[rid[0]][rid[1]]
	if boch_rid.get_resname() != ubch_rid.get_resname():
		print("[warning] the residues corresponding to %s are not the same in both files" % (str(rid)))
		print("[warning] resolving resn conflict by deleting nonoverlapping atoms")

	dcntbo = 0
	dcntub = 0
	common_anames = set(a.get_name() for a in boch_rid.get_atoms()) & set(a.get_name() for a in ubch_rid.get_atoms())
	for a in list(boch_rid.get_atoms()):
		if a.get_name() not in common_anames:
			boch_rid.detach_child(a.get_name())
			dcntbo += 1
	for a in list(ubch_rid.get_atoms()):
		if a.get_name() not in common_anames:
			ubch_rid.detach_child(a.get_name())
			dcntub += 1
	if dcntbo > 0 or dcntub > 0:
		print("[warning] for residue %s, deleted %d atoms from boch and %d atoms from ubch" % (str(rid), dcntbo, dcntub))

ubdict = {get_res_crid(x):x for x in unbound_struct.get_residues()}
bdict = {get_res_crid(x):x for x in bound_struct.get_residues()}

for fres in unbound_struct.get_residues():
	if not (get_res_crid(fres) in bdict):
		print("not found:")
		print(fres)
	elif bdict[get_res_crid(fres)].get_resname() != ubdict[get_res_crid(fres)].get_resname():
		print("resname mismatch:")
		print(ubdict[get_res_crid(fres)])
		print(bdict[get_res_crid(fres)])


#present = sorted(list(present), key=lambda x: x[1])
# only want residues present in both
# only want CA atoms
io = PDBIO()
pdbio = PDBIO()

def saveStructureToFile(struct, outfile, select=None):
	pdbio.set_structure(struct)
	if select is None:
		pdbio.save(outfile)
	else:
		pdbio.save(outfile, select)

def convertStructToPrody(struct, select=None):
	t = tempfile.mkdtemp()
	f = t+"/struct.pdb"
	
	saveStructureToFile(struct, f, select=select)

	pstruct = prody.parsePDB(f)

	shutil.rmtree(t)

	return pstruct

def calcRTBModesFromStruct(struct, blocks, cutoff=15.0, gamma=1.0, n_modes=None, select=None):
	tstart = time.time()
	if isinstance(struct, Bio.PDB.Entity.Entity):
		pstruct = convertStructToPrody(struct, select)
	else:
		pstruct = struct
	print("    prody-ified struct in %f seconds" % (time.time() - tstart))
	
	rtb = RTB("RTB analysis")
	tstart = time.time()
	rtb.buildHessian(pstruct, blocks, cutoff=cutoff)
	print("    ran rtb.buildHessian in %f seconds" % (time.time() - tstart))
	tstart = time.time()
	rtb.calcModes(n_modes)
	print("    ran rtb.calcModes in %f seconds" % (time.time() - tstart))

	return rtb

class CAAndPresSelect(Select):
	def accept_atom(self, atom):
		if atom.element != "H" and (get_res_crid(atom.get_parent()) in present):
			return 1
		else:
			return 0

class CASelect(Select):
	def accept_atom(self, atom):
		if atom.element == "H":
			return 0
		else:
			return 1

tmp = tempfile.mkdtemp()

saveStructureToFile(bound_struct, tmp+"/bound.pdb", CAAndPresSelect())
saveStructureToFile(unbound_struct, tmp+"/unbound.pdb", CAAndPresSelect())
saveStructureToFile(unbound_struct, tmp+"/unbound-whole.pdb", CASelect())

unf_bound_file = tmp+"/bound.pdb"
unf_unbound_file = tmp+"/unbound.pdb"
unf_unbound_whole_file = tmp+"/unbound-whole.pdb"

# is atom present in frag
atomMask = [(get_res_crid(r.get_parent()) in present) for r in unbound_struct.get_atoms()]


present_sele_string = "resnum " + (" ".join([str(rid[1]) for rid in present]))
print(present_sele_string)

atomVecMask = []
for a in atomMask:
	atomVecMask.append(a)
	atomVecMask.append(a)
	atomVecMask.append(a)


str_bound =         convertStructToPrody(bound_struct,select=CAAndPresSelect())
str_unbound =       convertStructToPrody(unbound_struct,select=CAAndPresSelect())
str_unbound_whole = convertStructToPrody(unbound_struct,select=CASelect())
str_receptor = prody.parsePDB(receptor_file)


ch_bound = str_bound
ch_unbound = str_unbound
print("ATOMCOUNTSQQ:")
print(ch_unbound.numAtoms())
print(ch_bound.numAtoms())
ch_unbound, t = prody.superpose(ch_unbound, ch_bound)

ch_unbound_whole = str_unbound_whole
t.apply(ch_unbound_whole)

t_ch_unbound_whole = t

prody.writePDB(outpref+".unbound.pdb", ch_unbound)
prody.writePDB(outpref+".unbound-whole.pdb", ch_unbound_whole)
prody.writePDB(outpref+".bound.pdb", ch_bound)
print("ATOMCOUNTS:")
print(ch_unbound.numAtoms())
print(ch_unbound_whole.numAtoms())
print(ch_bound.numAtoms())

prody_present = set(ch_unbound_whole.getResnums()) & set(ch_bound.getResnums())
prody_present_sele_str = "protein resnum " + (" ".join([str(x) for x in prody_present]))
print("SELE STR:")
print(prody_present_sele_str)
prody_atom_mask = [(x in prody_present) for x in ch_unbound_whole.getResnums()]
prody_atom_mask_ca = [((x.getResnum() in prody_present) and (x.getName() == "CA")) for x in ch_unbound_whole]
prody_coord_mask = []
prody_coord_mask_ca = []
for x in prody_atom_mask:
	prody_coord_mask.append(x)
	prody_coord_mask.append(x)
	prody_coord_mask.append(x)
for x in prody_atom_mask_ca:
	prody_coord_mask_ca.append(x)
	prody_coord_mask_ca.append(x)
	prody_coord_mask_ca.append(x)
print("len(prody_atom_mask) = %d" % len(prody_atom_mask))
print("len(prody_coord_mask) = %d" % len(prody_coord_mask))

if False:
	rrr = prody.calcRMSD(ch_bound, ch_unbound)
	print("PRODY RMSD = %f" % rrr)
	defvec = prody.calcDeformVector(ch_unbound, ch_bound)
	print("defvec shape:")
	print(defvec.getArray().shape)
	print(defvec.getArrayNx3().shape)
	print("MYNUMM: %d" % num_modes)
	anm_unbound_whole, sel = prody.calcANM(ch_unbound_whole, n_modes=None, selstr="protein")
	
	print("anm_unbound_whole.numModes() == %d" % anm_unbound_whole.numModes())

rmsds = []

def applyMask(mask, a):
	aarr = a.getArray()
	tmp = list(a)
	a = tmp
	ret = []

	for idx, mmm in zip(range(len(a)), a):
		mode = aarr[:,idx]
		rmode = []
		for v, m in zip(list(mode), mask):
			if m:
				rmode.append(v)
		rmode = prody.dynamics.mode.Vector(rmode)
		ret.append(rmode)
	return ret

def getChainTorsionAngles(struct):
	coords = struct.select("backbone").getCoords()
	print(coords)
	nrows, ncols = coords.shape
	#for row in coords:
	#	print row
	angles = []
	print(coords.shape)
	for i in range(nrows-3):
		quad = coords[i:i+4,:]
		args = [ vec for vec in quad ]
		#print args
		#print args[0].shape
		quad = tuple([Vector(v[0], v[1], v[2]) for v in quad])
		#print quad
		print(calc_dihedral(*quad)*180/3.1415926)
	print("")
	print("")
	print(coords)

def forceBackboneDistances(reference, struct):
	# reference and struct should have the same number of atoms in the same order

	resnums = reference.getResnums()
	resnames = reference.getNames()
	
	ref = reference
	ref_backbone = ref.select("protein and backbone")
	struct_backbone = struct.select("protein and backbone")

	# everything "before" a shifted atom stays in place
	# everything "after" a shifted atom inherits the same shift
	# To be shifted:
	#	N-CA
	#	CA-R
	#	CA-C
	#	C-O
	#	C-(next)

	


c = ch_unbound_whole.copy()
factor = 1.0



from itertools import groupby
from operator import itemgetter
def consecutive_groups(iterable, ordering=lambda x: x):
	for k, g in groupby(enumerate(iterable), key=lambda x: x[0] - ordering(x[1])):
		yield map(itemgetter(1), g)

def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[i:i + n]

def fitRigid(unbound, target):
	#prody.confProDy(verbosity="warning")
	prody_present = set(unbound.getResnums()) & set(target.getResnums())
	prody_present_sele_str = "protein resnum " + (" ".join([str(x) for x in prody_present]))

	target = target.select(prody_present_sele_str)
	unbound_part = unbound.select(prody_present_sele_str)

	unbound_part, t = prody.superpose(unbound_part, target) # align unbound part to target
	t.apply(unbound) # apply the same transformation to the whole unbound structure
	return unbound

def fitWithRTB(unbound, target, nmodes=None, factor=1.0):

	tstart = time.time()
	#prody.confProDy(verbosity="warning")
	print("    did prody config in %f seconds" % (time.time() - tstart))
	
	tstart = time.time()
	prody_present = set(unbound.getResnums()) & set(target.getResnums())
	#prody_present_sele_str = "protein resnum " + (" ".join([str(x) for x in prody_present]))
	sss = " ".join([str(x) for x in prody_present])
	sarr = sorted(map(int, sss.split()))
	seles = " ".join(["to".join((str(min(x)), str(max(x)))) for x in consecutive_groups(sarr)])
	prody_present_sele_str = "protein resnum " + seles

	target = target.select(prody_present_sele_str)
	unbound_part = unbound.select(prody_present_sele_str)
	print("    did prody sele in %f seconds" % (time.time() - tstart))

	tstart = time.time()
	unbound_part, t = prody.superpose(unbound_part, target) # align unbound part to target
	t.apply(unbound) # apply the same transformation to the whole unbound structure
	print("    did superpose in %f seconds" % (time.time() - tstart))
	
	tstart = time.time()
	atom_mask = [(x in prody_present) for x in unbound.getResnums()]
	atom_mask_ca = [((x.getResnum() in prody_present) and (x.getName() == "CA")) for x in ch_unbound_whole]
	coord_mask = []
	coord_mask_ca = []
	for x in prody_atom_mask:
		coord_mask.append(x)
		coord_mask.append(x)
		coord_mask.append(x)
	for x in prody_atom_mask_ca:
		coord_mask_ca.append(x)
		coord_mask_ca.append(x)
		coord_mask_ca.append(x)

	print("    set up masks in %f seconds" % (time.time() - tstart))

	blocks = list(ch_unbound_whole.getResnums())
	btmp = []
	bcnt = 0
	bcur = blocks[0]
	for b in blocks:
		if b != bcur:
			bcnt += 1
			bcur = b
		btmp.append(bcnt)
	
	#print blocks
	#print btmp
	blocks = btmp
	
	justblocks = list(sorted(set(blocks)))
	newblocks = list(chunks(justblocks, 3))
	
	bmap = {}
	for chunk, i in zip(newblocks, range(len(newblocks))):
		for bid in chunk:
			bmap[bid] = i

	#blocks = [bmap[b] for b in blocks]
	
	numBlocks = len(set(blocks))

	class MyGamma(Gamma):
		def __init__(self, atoms):
			self.atoms = atoms
			self.resnums = atoms.getResnums()
			self.names = atoms.getNames()
		def gamma(self, dist2, i, j):
			resnums = self.resnums
			names = self.names

			if resnums[i] == resnums[j]:
				return 1.0
			elif resnums[i] == (resnums[j]+1) and names[i] == "N" and names[j] == "CA":
				return 1.0
			elif resnums[j] == (resnums[i]+1) and names[j] == "N" and names[i] == "CA":
				return 1.0
			else:
				pass
				

			return 1.0

	tstart = time.time()
	rtb = calcRTBModesFromStruct(unbound, blocks, n_modes=nmodes, gamma=MyGamma(unbound))
	print("    ran calcRTBModesFromStruct in %f seconds" % (time.time() - tstart))
	
	tstart = time.time()
	npa = numpy.array(applyMask(coord_mask, rtb[:rtb.numModes()]))
	tend = time.time()
	tdiff = tend - tstart
	print("    applied mask in %f seconds" % tdiff)

	defvec = prody.calcDeformVector(unbound_part, target)
	weights = npa * defvec

	weights *= factor # scale the weights (e.g. factor<1.0 can be used to mitigate overfitting)

	#make npa be the full modes again
	npa = numpy.array(list(rtb[:rtb.numModes()]))
	approx = sum(weights * npa)
	unbound.setCoords(unbound.getCoords() + approx.getArrayNx3())

	return unbound

def getCACADistances(struct):
	struct = struct.select("ca")
	arr = struct.getCoords() # Nx3
	(m, n) = arr.shape
	diff = arr[0:m-2,:] - arr[1:m-1,:]
	dists = [numpy.linalg.norm(vec) for vec in diff]
	return dists

c = ch_unbound_whole.copy()

orig_dists = getCACADistances(c)
factor = 0.05
use_n_modes = 20
use_n_iter = 500
mytar = tarfile.open("this.%s.rtb.tar" % os.environ["thisname"], "w")

prody.confProDy(verbosity="warning")
for i in range(use_n_iter):
	sys.stdout.flush()
	print("RTB FIT ITER %d" % i)

	orig_c = c.copy()
	
	myrec = None
		
	if i <= 100:
	#if i <= 20:
		myrec = None
	elif i % 10 == 0:
		myrec = str_receptor.copy()

	notches = [factor*x for x in [1.0]]
	if i < 5: 
		notches = [factor*x for x in [1.0]]

	if i == 5:
		orig_dists = getCACADistances(c) 
	

	for fff in notches:
		print("  invoking fitWithRTB with factor = %f..." % fff)
		tstart = time.time()
		c = fitWithRTB(c, ch_bound, nmodes=use_n_modes, factor=fff)
		tend = time.time()
		tdiff = tend - tstart
		print("    finished in %f seconds" % tdiff)
	
		print("  invoking doPhenixGeometryMinimization...")
		tstart = time.time()
		
		if True:
			mini_c = doPhenixGeometryMinimization(c, receptor=myrec)
		else:
			mini_c = c.copy()
		tend = time.time()
		tdiff = tend - tstart
		print("    finished in %f seconds" % tdiff)
		
		c = mini_c
		
		print("  realigning...")
		c = fitRigid(c, ch_bound)

		#forceBackboneDistances(orig_c, c)
		dists = getCACADistances(c)
		if any([(cur > (orig+0.1)) for (orig, cur) in zip(orig_dists, dists)]):
			#c = orig_c
			nviol = sum([1 if (cur > (orig+0.1)) else 0 for (orig, cur) in zip(orig_dists, dists)])
			print("  fitting with factor = %f violated CA-CA distance constraint (%d violations)" % (fff, nviol))
			#print ("ORIG", orig_dists)
			#print ("CURR", dists)
			sys.stdout.flush()
			continue
		else:
			break

	if use_n_modes is None:
		curfname = "%s.rtb.xfact.%.2f.with.all.modes.%03d.pdb" % (outpref, factor, i)
	else:
		curfname = "%s.rtb.xfact.%.2f.with.%d.modes.%03d.pdb" % (outpref, factor, use_n_modes, i)
	prody.writePDB(curfname, c)
	
	# always set lastiter so we dont have to regen later
	shutil.copy(curfname, "tempiter.pdb")
	os.rename("tempiter.pdb", "lastiter.pdb")

	mytar.add(curfname)
	os.remove(curfname)
	


# clean up temp dir
shutil.rmtree(tmp)

exit()

