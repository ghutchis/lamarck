#!/usr/bin/python

import sys
import random
import pybel
import openbabel

datafile = open(sys.argv[1])

def hCount(atom):
	return atom.ImplicitHydrogenCount() + atom.ExplicitHydrogenCount()

def globalopt(mol, debug=False, fast=False):
	pybel._builder.Build(mol.OBMol)
	mol.addh()

	# if fast, skip the forcefield cleanup
	if not fast:
		ff = pybel._forcefields["mmff94"]
		success = ff.Setup(mol.OBMol)
		if not success:
			ff = pybel._forcefields["uff"]
			success = ff.Setup(mol.OBMol)
		if not success:
			sys.exit("Cannot set up forcefield")

		ff.SteepestDescent(500, 1.0e-4)	
		ff.WeightedRotorSearch(200, 25)
		ff.ConjugateGradients(500, 1.0e-6)
		ff.GetCoordinates(mol.OBMol)

i = 0
for line in datafile:
		# get the smiles string
		monomer = line.rstrip()
		print "monomer: ", monomer
		
		# get the actual molecule
		try:
			molecule = pybel.readstring('smi', monomer)
		except IOError:
			continue
		globalopt(molecule) # now have 3D coords here
		


		mutationType = random.choice([1,2])
		
		if mutationType == 1: # Functional group "growing"
			numAtoms = molecule.OBMol.NumAtoms()
			# pick a random atom, not first or last
			for i in range (0,4): # 5 times:
				randomAtomIdx = random.randint(1, numAtoms -1)
				randomAtom = molecule.atoms[randomAtomIdx]
				# Make sure that the chosen atom has an open valence
				if (hCount(randomAtom.OBAtom) >= 1):
					break
			# after the for loop
			if (hCount(randomAtom.OBAtom) < 1):
				continue # move to the next smiles
							
			functionalGroups = ['(O)', '(C)', '(C#N)', '(F)', '(Cl)', '(N)', '(P)', '([Si])'] 
			functionalGroups += ['(NC)', '(SC)', '(CC)', '(S)', '(=O)', '(=S)', '(N(=O)=O)']
			functionalGroups += ['(C(F)(F)F)', '(OC)', '(C=O)', '(C=C)', '(C#C)']
			functionalGroups += ['(c(s1)ccc1)', '(c(n1)ccc1)', '(c(o1)ccc1)', '(C(=O)N)']
			functionalGroups += ['(c1ccccc1)', '(C(=O)O)', '(C1C=CC=C1)']
			# get a molecule corresponding to the random functional group
			group = pybel.readstring('smi', random.choice(functionalGroups))
			globalopt(group) # now have 3D coords here
			# we'll want to attach to atom #1 in the group, and randomAtom in the molecule
			molecule.OBMol += group.OBMol
			
			# Open Babel indexes atoms from 1 (sigh!)
			pybel._builder.Connect(molecule.OBMol, randomAtomIdx + 1, numAtoms + 1, 1)
			# we probably have 2 extra hydrogens at the attachment point
			molecule.removeh()
			molecule.addh()
			print  "SMILES: ", molecule.write()
			
		elif mutationType == 2: # change an atom
			numAtoms = molecule.OBMol.NumAtoms()
			# pick a random atom, not first or last
			randomAtomIdx = random.randint(1, numAtoms -1)
			randomAtom = molecule.atoms[randomAtomIdx]
			
			totalBonds = randomAtom.OBAtom.BOSum()
			
			if (totalBonds == 4):
				fourValentElements = [6, 14, 32]
				randomElement = random.choice(fourValentElements)
				# Make sure that the random element you selected is actually 
				# a different element than is what currently in molecule
				while (randomAtom == randomElement):
					 randomElement = random.choice(fourValentElements)
				randomAtom.OBAtom.SetAtomicNum(randomElement)
				print molecule
			elif (totalBonds == 3):
				threeValentElements = [6, 7, 14, 15, 32] # skip 33 b/c of toxicity
				randomElement = random.choice(threeValentElements)
				while (randomAtom == randomElement):
					 randomElement = random.choice(threeValentElements)
				randomAtom.OBAtom.SetAtomicNum(randomElement)			
				print molecule
			elif (totalBonds == 2):
				twoValentElements = [6, 7, 8, 14, 15, 16, 32, 34]
				randomElement = random.choice(twoValentElements)
				while (randomAtom == randomElement):
					 randomElement = random.choice(twoValentElements)
				randomAtom.OBAtom.SetAtomicNum(randomElement) 
				print molecule
			elif (totalBonds == 1):
				oneValentElements = [6, 7, 8, 9, 14, 15, 16, 17, 32, 34, 35]
				randomElement = random.choice(oneValentElements)
				while (randomAtom == randomElement):
					 randomElement = random.choice(oneValentElements)
				randomAtom.OBAtom.SetAtomicNum(randomElement) 
				print molecule

# ##Broken: need to figure out how to get it to pick a ring atom randomly.
# 		elif mutationType == 3: # Delete an atom from a ring
# 			numAtoms = molecule.OBMol.NumAtoms()
# 			# pick a random atom, not first or last
# 			randomAtomIdx = random.randint(1, numAtoms -1)
# 			randomAtom = molecule.atoms[randomAtomIdx]
# 			
# 			# check to make sure that the molecule is large enough to delete a ring atom
# 			ringSize = randomAtom.OBAtom.CountRingBonds()
# 			while ringSize <= 5:
# 				mutationType = random.choice ([1,2])
# 			else:
# 				continue
# 			# pick a random atom, not first or last
# 			randomAtomIdx = random.randint(1, numAtoms -1)
# 			randomAtom = molecule.atoms[randomAtomIdx]
# 			# Find an atom with open valence
# 			while (hCount(randomAtom.OBAtom) < 1):
# 				randomAtom = random.choice(molecule.atoms)
# 			monomer = monomer[:randomAtom] + monomer[(randomAtom):]
# 			
# 
# 		elif mutationType == 4: # Add an atom in a ring
# 			numAtoms = molecule.OBMol.NumAtoms()
# 			# pick a random atom, not first or last
# 			randomAtomIdx = random.randint(1, numAtoms -1)
# 			randomAtom = molecule.atoms[randomAtomIdx]
# 			ringElements = [6, 7, 14, 15, 32] 
# 			randomElement = random.choice(ringElements)
# 			while (randomAtom == randomElement):
# 				randomElement = random.choice(ringElements)
# 			# Does not work...need to find how to add a atom into the molecule
# 			randomAtom.OBMol.InsertAtom(randomElement)
# 			def globalopt(mol, debug=False, fast=False):
# 					pybel._builder.Build(mol.OBMol)
# 					mol.addh()
# 			print  molecule


		try:
			molecule = pybel.readstring('smi', monomer)
		except IOError:
			print "except"
			continue

