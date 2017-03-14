import math
import random
import time
import sys
import getopt
from constants import *
from utilities import *

def GenerateSupercell(runOptions,particleOptionsA,mixRatio = 1,
																	particleOptionsB = None):
	#PrimBox, repeat, primfrac
	coords = []
	if particleOptionsB == None:
		particleOptionsB = particleOptionsA
	fractionA = mixRatio
	fractionB = 1-mixRatio
	PrimBox = runOptions._PrimBox
	repeat = runOptions._repeat
	particles=[]
	n = 0
	ncell = 0
	xorigin = 0
	mass,sigma,epsilon,name = 0,0,0,""
	typeA = True
	totCells = repeat**3 * len(primfrac)
	totA = mixRatio*totCells
	for i in range(repeat) :
		yorigin = 0.0
		for j in range(repeat) :
			zorigin = 0.0
			for k in range(repeat) :
				ncell = ncell + 1
				if n<(totA):
					name = particleOptionsA._name
					mass = particleOptionsA._mass
					sigma = particleOptionsA._sigma
					epsilon = particleOptionsA._epsilon
				else:
					name = particleOptionsB._name
					mass = particleOptionsB._mass
					sigma = particleOptionsB._sigma
					epsilon = particleOptionB._epsilon
				for w in range(len(primfrac)):
					xcoord = PrimBox*primfrac[w][0] + xorigin
					ycoord = PrimBox*primfrac[w][1] + yorigin
					zcoord = PrimBox*primfrac[w][2] + zorigin
					particles.append(Particle(mass,sigma,epsilon,name = name))
					coords.append([xcoord,ycoord,zcoord])
					n = n + 1
				zorigin += PrimBox
			yorigin += PrimBox
		xorigin += PrimBox
					
	print "SuperCell Generated containing ",ncell," primitive cells and ",n," atoms in box side "
	return coords,particles
	
"""def MD_Loop(particles, coords_t, coords_t_minus, gr_histogram, options) :
	print '# {0:4s}	{1:12s}	{2:12s}	{3:12s}	{4:6s}'.format("Step","Pot.Energy","Kin.Energy","Tot.Energy","T(Inst)")
	n = len(particles)
	XYZfile = options._XYZfile
	nsteps = options._nsteps
	timestep = options._timestep
	Lbox = options._Lbox	
"""

def Maxwell_Boltzmann(coords, particles, kBT) :
		n = len(coords)
		Velocity = [[0,0,0] for i in range(n) ]
		for i in range(n) :
				Root_kTm = math.sqrt(kBT/particles[i]._mass)
				for k in range(3) :
						Velocity[i][k] =	Root_kTm * random.gauss(0.0, 1.0)

		return Velocity

#
# Compute energy and forces of system of Lennard-Jones atoms in cubic box
#
def LJ_Energy(coords, particles, Lbox) :
		n = len(coords)
		energy	= 0.0
		rij=[0,0,0]
		forces = [[0,0,0] for i in range(n) ]
		for i in range(n) :
				for j in range(n):
					if not i==j:
						dist_sq = 0.0
						for k in range(3) :
								rcomp = coords[j][k] - coords[i][k]
								rcomp -= Lbox*((rcomp/Lbox+0.5)//1)
								rij[k] = rcomp
								dist_sq += rcomp**2
						if( dist_sq <= Lbox**2/4.0 ) :
								rdist = (particles[i]._sigma**2/dist_sq)**3
								energy += 4.0*particles[i]._epsilon*(rdist**2 - rdist)/2.0
								force = 	24.0*particles[i]._epsilon * (2.0*rdist**2 -rdist)/dist_sq/2.0
								for k in range(3) :
										forces[i][k] -= force*rij[k]
										forces[j][k] += force*rij[k]
								#print "Force=",force
		return energy, forces

#
# Step the atomic co-ordinates using the Verlet algorithm
#
def Verlet_Step(timestep, coords_t, coords_t_minus, forces,particles) :

		n = len(coords_t)
		coords_t_plus = [ [0,0,0] for i in range(0, n)]
		for i in range(n) :
				for k in range(3) :
						coords_t_plus[i][k] = 2.0*coords_t[i][k] -	coords_t_minus[i][k] + timestep**2*forces[i][k]/particles[i]._mass
		return coords_t_plus
#
# Calculate an initial set of atomic coordinates at (t-delta t) to start Verlet integration
#
def Verlet_Init_Backstep(timestep, coords_t, velocity, forces,particles):
		n = len(coords_t)
		coords_t_minus = [ [0,0,0] for i in range(0, n)]
		for i in range(n) :
				for k in range(3) :
						coords_t_minus[i][k] = coords_t[i][k] -	timestep*velocity[i][k] + 0.5*timestep**2*forces[i][k]/particles[i]._mass
		return coords_t_minus
#
# Compute the velocity in the Verlet algorithm
#
def Verlet_Velocity(timestep, coords_t_plus, coords_t_minus):
		n = len(coords_t_plus)
		Velocity = [[0,0,0] for i in range(n) ]
		for i in range(n) :
				for k in range(3) :
						Velocity[i][k] =	(coords_t_plus[i][k] - coords_t_minus[i][k]) / (2.0 * timestep)
		return Velocity
#
# Compute the instantaneous kinetic energy from the velocities
#
def Kinetic_Energy(velocity,particles) :
	n = len(velocity)
	Energy = 0.0
	for i in range(n) :
		Energy += 0.5*particles[i]._mass*(velocity[i][0]**2 + velocity[i][1]**2 + velocity[i][2]**2)
	return Energy
		
def Accumulate_RDF(coords, particles, Lbox, histogram,typeA = "",typeB = ""):
		n = len(coords)
		nbins = len(histogram)
		for i in range(n) :
			if (particles[i]._name == typeA or (typeA == "" and typeB == "") ):
				for j in range(i+1,n):
					if (particles[i]._name == typeB or (typeA == "" and typeB == "") )
						dist_sq = 0.0
						for k in range(3) :
								rcomp = coords[j][k] - coords[i][k]
								rcomp = rcomp - Lbox*((rcomp/Lbox+0.5)//1)
								dist_sq += rcomp**2

						if( dist_sq < Lbox**2/4.0 ) :
								bin = int(math.sqrt(4.0*dist_sq/Lbox**2) * nbins)
								try:
										histogram[bin] += 1
								except:
										print "RDF Error - r=", math.sqrt(dist_sq),",	bin=", bin

def Calculate_RDF(n, Lbox, histogram) :
		
		nbins = len(histogram)-1
		bin_width = Lbox/(2*nbins)
		rho = n/Lbox**3

		sum = 0
		for i in range(nbins) :
				sum += histogram[i]

		fileh = open('gr.dat','w')

		fileh.write( '# {0:6s}	{1:9s}\n'.format("r","g(r)"))
		for i in range(nbins) :
			r = (i+0.5)*Lbox/(2.0*nbins)
			rdf = histogram[i] /(4.0*math.pi*r**2 * bin_width * rho) * (n-1) / sum
			fileh.write( '{0:6f}	{1:9.6f}\n'.format(r,rdf))

#
# Run "nsteps" MD steps 
#
def MD_Loop(coords_t, coords_t_minus, particles,gr_histogram,options) :
	 	XYZfile = options._XYZfile
		nsteps = options._nsteps
		timestep = options._timestep
		Lbox = options._Lbox	
		n = len(coords_t)
		XYZfh = None
		print '# {0:4s}  {1:12s}  {2:12s}  {3:12s}  {4:6s}'.format("Step","Pot.Energy","Kin.Energy","Tot.Energy","T(Inst)")
		if XYZfile != ""	:
				XYZfh = open(XYZfile,'w')

		for istep in range(nsteps) :
				Pot_Energy, forces = LJ_Energy(coords_t, particles, Lbox)
				Accumulate_RDF(coords_t, Lbox, gr_histogram)
				coords_t_plus = Verlet_Step(timestep, coords_t, coords_t_minus, forces, particles)
				Velocity = Verlet_Velocity(timestep, coords_t_plus, coords_t_minus)
				Kin_Energy = Kinetic_Energy(Velocity,particles)
				T_Inst = evk* 2.0*Kin_Energy/(3*n - 3)
				coords_t_minus = coords_t
				coords_t = coords_t_plus
				if XYZfile != "":
					XYZfh.write( '{0:4d}	{1:12.7f}	{2:12.7f}	{3:12.7f}	{4:6.2f}'.format(istep,
															Pot_Energy,Kin_Energy,Pot_Energy+Kin_Energy,T_Inst)
				print '{0:4d}	{1:12.7f}	{2:12.7f}	{3:12.7f}	{4:6.2f}'.format(istep,Pot_Energy,Kin_Energy,Pot_Energy+Kin_Energy,T_Inst)
		if XYZfile != "":
				XYZfh.close()
				
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
