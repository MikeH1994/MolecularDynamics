import pymold
import math
import random
import time
import sys
import getopt
from constants import *
from utilities import *
from pymold import *

def run():
	###############################################################
	PrimBox = 5.7064
	timestep = 0.01e-12/tunit
	initial_temperature = 84
	repeat = 3
	XYZfile = ""
	nsteps = 1000
	###############################################################
	try:
		opts, args = getopt.getopt(sys.argv[1:],"ht:n:d:r:x:",["help","temperature=","nsteps=","timestep=","box-repeat=","xyz-file"])
	except getopt.GetoptError:
		print 'pymold.py -h -t <temperature> -n <nsteps> -d <timestep> -r <box repeat> -x <XYZ file>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'pymold.py -h -t <temperature> -n <nsteps> -d <timestep> -r <box repeat> -x <XYZ file>'
			sys.exit()
		elif opt in ("-t", "--temperature"):
			initial_temperature = float(arg)
		elif opt in ("-n", "--nsteps"):
			nsteps = int(arg)
		elif opt in ("-d", "--timestep"):
			timestep = float(arg) * 1.0e-12 / tunit
			print "Timestep set to:",timestep, "(",tunit,")"
		elif opt in ("-r", "--box-repeat"):
			repeat = int(arg)
		elif opt in ("-x", "--xyz-file"):
			XYZfile = arg	
	###############################################################
	particleOptionsA = ParticleOptions("Ar")
	particleOptionsB = ParticleOptions("Kr")
	###############################################################
	ratioA = 1
	ratioB = 0
	runOptions = RunOptions(PrimBox = PrimBox, timestep = timestep, 
		repeat = repeat, XYZfile = XYZfile, nsteps = nsteps,ratioA = ratioA,
		ratioB = ratioB, partAName = particleOptionsA._name, partBName = particleOptionsB._name)
	###############################################################
	coords,particles=	GenerateSupercell(runOptions,particleOptionsA,
																	particleOptionsB = particleOptionsB,ratioA = ratioA, ratioB = ratioB)
	kBT = kB/e*initial_temperature
	Pot_Energy, forces = LJ_Energy(coords, particles, runOptions._Lbox)
	random.seed(2.0)
	velocity = Maxwell_Boltzmann(coords, particles, kBT)
	n=len(coords)
	print "Lennard-Jones energy = ", Pot_Energy, "eV = ",Pot_Energy*kjmol, "kJ/mol"

	MD_Loop(coords,velocity,particles,runOptions,kBT)

if __name__ == "__main__":
	run()
