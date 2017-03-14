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
	nbins = 50
	gr_histogram = [ 0 for i in range(nbins+1)]
	###############################################################
	PrimBox = 5.7064
	timestep = 0.01e-12/tunit
	initial_temperature = 84
	repeat = 3
	XYZfile = "../ComparisonNew.txt"
	nsteps = 20000
	runOptions = RunOptions(PrimBox = PrimBox, timestep = timestep, 
		repeat = repeat, XYZfile = XYZfile, nsteps = nsteps)
	ratioA = 1
	ratioB = 0
	###############################################################
	particleOptionsA = ParticleOptions("Ar")
	particleOptionsB = ParticleOptions("Kr")
	###############################################################
	coords_t,particles = GenerateSupercell(runOptions, particleOptionsA)
	Pot_Energy, forces = LJ_Energy(coords_t, particles, runOptions._Lbox)
	random.seed(2.0)
	Velocity = Maxwell_Boltzmann(coords_t,particles,kB/e*initial_temperature)
	coords_t_minus = Verlet_Init_Backstep(timestep,coords_t,Velocity, forces, particles)
	gr_histogram = [ 0 for i in range(nbins+1) ]
	print "Lennard-Jones energy = ", Pot_Energy, "eV = ",Pot_Energy*kjmol, "kJ/mol"
	MD_Loop(coords_t, coords_t_minus, particles, gr_histogram, runOptions)
"""
coords=GenerateSupercell(PrimBox, repeat, primfrac)
n=len(coords)
forces=[]
Pot_Energy, forces = LJ_Energy(coords, Lbox, epsilon, sigma)
Velocity = Maxwell_Boltzmann(coords, mass, kB/e*initial_temperature)
coords_t_minus = Verlet_Init_Backstep(timestep, mass, coords, Velocity, forces)
print "Epsilon = ",epsilon," eV;   sigma = ", sigma, "A"
print "Lennard-Jones energy = ", Pot_Energy, "eV = ",Pot_Energy*kjmol, "kJ/mol"
gr_histogram = [ 0 for i in range(nbins+1)]
MD_Loop(mass, timestep, nsteps, coords, coords_t_minus, gr_histogram)
Calculate_RDF(len(coords), Lbox, gr_histogram)
"""
"""
def GenerateSupercell(runOptions,particleOptionsA,mixRatio = 1,
																	particleOptionsB = None):
	return particles,coords
def Maxwell_Boltzmann(coords, particles, kBT) :
	return Velocity
def LJ_Energy(coords, particles, Lbox) :
	return energy, forces
def Verlet_Step(timestep, coords_t, coords_t_minus, forces,particles) :
	return coords_t_plus
def Verlet_Init_Backstep(timestep, coords_t, velocity, forces,particles) :
	return coords_t_minus
def Verlet_Velocity(timestep, coords_t_plus, coords_t_minus) :
	return Velocity
def Kinetic_Energy(velocity,particles) :
	return Energy
def Accumulate_RDF(coords, Lbox, histogram) :
def Calculate_RDF(n, Lbox, histogram) :
def WriteXYZ(particles, coords, Energy, funit) :
def MD_Loop(coords_t, coords_t_minus, particles,gr_histogram,options) :
"""	

def commandLineArgs():
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


if __name__ == "__main__":
	commandLineArgs()
	run()
