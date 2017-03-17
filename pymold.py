import math
import random
import time
import sys
import getopt
from constants import *
from utilities import *

def GenerateSupercell(runOptions,particleOptionsA,
																	particleOptionsB = None,ratioA = 1, ratioB = 0):
	coords = []
	if particleOptionsB == None:
		particleOptionsB = particleOptionsA
	PrimBox = runOptions._PrimBox
	repeat = runOptions._repeat
	particles=[]
	n = 0
	ncell = 0
	xorigin = 0
	mass,sigma,epsilon,name = 0,0,0,""
	totParticles = repeat**3 * len(primfrac)
	swapIndex = 0
	swapEvery = ratioA
	typeA = True
	for i in range(repeat) :
		yorigin = 0.0
		for j in range(repeat) :
			zorigin = 0.0
			for k in range(repeat) :
				ncell = ncell + 1
				for w in range(len(primfrac)):
					if swapIndex>=swapEvery and not ratioB==0:
						typeA = not typeA
						swapIndex = 0
						if typeA or ratioB:
							swapEvery = ratioA
						else:
							swapEvery = ratioB
					if typeA:
						name = particleOptionsA._name
						mass = particleOptionsA._mass
						sigma = particleOptionsA._sigma
						epsilon = particleOptionsA._epsilon
					else:
						name = particleOptionsB._name
						mass = particleOptionsB._mass
						sigma = particleOptionsB._sigma
						epsilon = particleOptionsB._epsilon
					xcoord = PrimBox*primfrac[w][0] + xorigin
					ycoord = PrimBox*primfrac[w][1] + yorigin
					zcoord = PrimBox*primfrac[w][2] + zorigin
					particles.append(Particle(mass,sigma,epsilon,name = name))
					coords.append([xcoord,ycoord,zcoord])
					
					n = n + 1
					swapIndex+=1
				zorigin += PrimBox
			yorigin += PrimBox
		xorigin += PrimBox
	nA = 0
	nB = 0		
	for i in range(n):
		#print(i,particles[i]._name)
		if particles[i]._name == particleOptionsA._name:
			nA+=1
		else:
			nB+=1
	print('Particle {}: {}/{} Particle {}: {}/{}'.format(particleOptionsA._name,nA,totParticles, particleOptionsB._name,nB,totParticles))
	print "SuperCell Generated containing ",ncell," primitive cells and ",n," atoms in box side "
	return coords,particles
	

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
	#				rcomp -= Lbox*math.floor(rcomp/Lbox+0.5)
					rcomp -= Lbox*((rcomp/Lbox+0.5)//1)
					rij[k] = rcomp
					dist_sq += rcomp**2

				if( dist_sq <= Lbox**2/4.0 ) :
					rdist = (particles[i]._sigma**2/dist_sq)**3
					energy += 4.0*particles[i]._epsilon*(rdist**2 - rdist)/2.0
					force = 24.0*particles[i]._epsilon * (2.0*rdist**2 -rdist)/dist_sq/2.0
					for k in range(3) :
						forces[i][k] -= force*rij[k]
						forces[j][k] += force*rij[k]
					#print "Force=",force
	return energy, forces

#
# Step the atomic co-ordinates and Nose-Hoover variables using the Velocity Verlet algorithm
#
def VVerlet_NH_Step_1(timestep, particles, NHQ, kBT, KE, coords, velocities, forces, zeta, lns) :

	n = len(coords)
	coords_t_plus = [ [0,0,0] for i in range(0, n)]
	velocities_t_half = [ [0,0,0] for i in range(0, n)]
	for i in range(n) :
		for k in range(3) :
			coords_t_plus[i][k] = coords[i][k] + timestep*velocities[i][k] + timestep**2*(forces[i][k] - zeta*velocities[i][k])/(2.0*particles[i]._mass)
			velocities_t_half[i][k] = velocities[i][k] + timestep*(forces[i][k] - zeta*velocities[i][k])/(2.0*particles[i]._mass)

	zeta_t_half = zeta + timestep/(2.0*NHQ)*(KE - (3*n)*kBT/2.0)
	lns_t_plus = lns + timestep*zeta_t_half + timestep**2/NHQ*(KE - (3*n)*kBT/2.0)
		
	return lns_t_plus, zeta_t_half, coords_t_plus, velocities_t_half
#
# Step the atomic co-ordinates and Nose-Hoover variables using the Velocity Verlet algorithm
# An approximate scheme is used for the second velocity half-step.
#
def VVerlet_NH_Step_2(timestep, particles, NHQ, kBT, KE, velocities, forces_t_plus, zeta) :
	
	n = len(velocities)
	velocities_t_plus = [ [0,0,0] for i in range(0, n)]
	for i in range(n) :
		for k in range(3) :
			velocities_t_plus[i][k] = (velocities[i][k] + timestep*forces_t_plus[i][k]/(2.0*particles[i]._mass)) \
										/ (1 + (timestep*zeta/2.0))


	zeta_t_plus = zeta + timestep/(2.0*NHQ)*(KE - (3*n)*kBT/2.0)

	return zeta_t_plus, velocities_t_plus
#
# Compute the instantaneous kinetic energy from the velocities
#
def Kinetic_Energy(particles, velocity) :

	n = len(velocity)

	Energy = 0.0
	for i in range(n) :
		Energy += 0.5*particles[i]._mass*(velocity[i][0]**2 + velocity[i][1]**2 + velocity[i][2]**2)	

	return Energy
		
def Accumulate_RDF(coords, Lbox, histogram):
	n = len(coords)
	nbins = len(histogram)
	for i in range(n) :
		for j in range(i+1,n):
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

def Partial_RDF(coords, particles,Lbox, histogram,fromAtom, toAtom):
	n = len(coords)
	nbins = len(histogram)
	increment = 1
	if fromAtom==toAtom:
		increment = 0.5
	for i in range(n) :
		if particles[i]._name == fromAtom:
			for j in range(n):
				if not i==j:
					if particles[j]._name == toAtom:
						dist_sq = 0.0
						for k in range(3) :
							rcomp = coords[j][k] - coords[i][k]
							rcomp = rcomp - Lbox*((rcomp/Lbox+0.5)//1)
							dist_sq += rcomp**2
							if( dist_sq < Lbox**2/4.0 ) :
								bin = int(math.sqrt(4.0*dist_sq/Lbox**2) * nbins)
								try:
									histogram[bin] += increment
								except:
									print "RDF Error - r=", math.sqrt(dist_sq),",	bin=", bin


def Calculate_RDF(n, Lbox, histogram,fileh) :
	fileh = open(fileh,'w')
	nbins = len(histogram)-1
	bin_width = Lbox/(2*nbins)
	rho = n/Lbox**3
	sum = 0
	for i in range(nbins) :
		sum += histogram[i]
	fileh.write( '# {0:6s}	{1:9s}\n'.format("r","g(r)"))
	for i in range(nbins) :
		r = (i+0.5)*Lbox/(2.0*nbins)
		rdf = histogram[i] /(4.0*math.pi*r**2 * bin_width * rho) * (n-1) / sum
		fileh.write( '{0:6f}	{1:9.6f}\n'.format(r,rdf))
	fileh.close()
		
def Calculate_Partial_RDF(n, Lbox, nA,nB,t,histogram,fileh) :
	fileh = open(fileh,'w')
	nbins = len(histogram)-1
	bin_width = Lbox/(2*nbins)
	totParticles = runOptions._repeat**3 * len(primfrac)
	v = Lbox**3
	sum = 0
	for i in range(nbins) :
		sum += histogram[i]
	fileh.write( '# {0:6s}	{1:9s}\n'.format("r","g(r)"))
	for i in range(nbins) :
		r = (i+0.5)*Lbox/(2.0*nbins)
		rdf = 3*v*histogram[i] /(4.0*math.pi*nA*nB*t*((r+bin_width)**3-r**3))
		fileh.write( '{0:6f}	{1:9.6f}\n'.format(r,rdf))
	fileh.close()
		
#
# Run "nsteps" MD steps 
#
#def MD_Loop(mass, timestep, nsteps, coords, velocity, gr_histogram, NHQ, kBT) :
def MD_Loop(coords,velocity,particles,options,kBT) :
	    print '# {0:4s}  {1:12s}  {2:12s}  {3:9s}  {4:9s}  {5:12s}  {6:6s}'.format("Step",
           "Pot.Energy","Kin.Energy","NHzeta","s","Tot.Energy","T(Inst)")
	n = len(coords)
	nbins = 50
	gr_histogram = [ 0 for i in range(nbins+1)]	
	gr_aa_histogram = [ 0 for i in range(nbins+1)]	
	gr_ab_histogram = [ 0 for i in range(nbins+1)]
	gr_bb_histogram = [ 0 for i in range(nbins+1)]		
	###############################################
	#get run options
	XYZfile = options._XYZfile
	nsteps = options._nsteps
	timestep = options._timestep
	Lbox = options._Lbox
	nameA = options._particleAName
	nameB = options._particleBName
	NHQ = options._NHQ
	###############################################
	XYZfh = None
	if XYZfile != "":
		XYZfh = open(XYZfile,'w')
		XYZfh.write('{0:4s}	{1:12s}	{2:12s}	{3:12s}	{4:6s}'.format("Step","Pot.Energy","Kin.Energy","Tot.Energy","T(Inst)"))
	###############################################	
	zeta = 0
	lns = 0
	#
	# Initialize output file for xyz if requested
	#
	if XYZfile != ""	:
			XYZfh = open(XYZfile,'w')

	Pot_Energy, forces = LJ_Energy(coords, particles, Lbox)
	Kin_Energy = Kinetic_Energy(particles, velocity)

	for istep in range(nsteps) :
		lns, zeta, coords_t_plus, velocity_t_plus =	VVerlet_NH_Step_1(timestep, particles, NHQ, kBT, Kin_Energy,coords, velocity, forces, zeta, lns)
		Accumulate_RDF(coords, Lbox, gr_histogram)				
		Pot_Energy, forces_t_plus = LJ_Energy(coords_t_plus,particles, Lbox)
		Kin_Energy = Kinetic_Energy(particles, velocity_t_plus)
		zeta, velocity_t_plus = VVerlet_NH_Step_2(timestep, particles, NHQ, kBT, Kin_Energy,velocity_t_plus, forces_t_plus, zeta )
		T_Inst = evk* 2.0*Kin_Energy/(3*n - 3)
		forces = forces_t_plus
		velocity = velocity_t_plus
		coords = coords_t_plus
		Tot_Energy = Pot_Energy + Kin_Energy + 0.5*NHQ*zeta**2 + (3*n)*kBT*math.exp(lns)

		if istep%5==0:
			print '{0:4d}	{1:12.7f}	{2:12.7f}	{3:9.6f}	{4:9.6f}	{5:12.7f}	{6:6.2f}'.format(istep,
			Pot_Energy,Kin_Energy, zeta, math.exp(lns), Tot_Energy,T_Inst)
			if XYZfile != "":
				XYZfh.write( '	{0:4d}	{1:12.7f}	{2:12.7f}	{3:12.7f}	{4:6.2f}'.format(istep,
																Pot_Energy,Kin_Energy,Pot_Energy+Kin_Energy,T_Inst))
	if XYZfile != "":
		XYZfh.close()
	mixRatio = options._ratioA/(options._ratioA + options._ratioB)
	totParticles = options._repeat**3 * len(primfrac)
	nA = mixRatio*totParticles
	nB = (1-mixRatio)*totParticles
	#Calculate_Partial_RDF(n, Lbox, nA,nB,nsteps*timestep,histogram,fileh)
	
"""
def GenerateSupercell(runOptions,particleOptionsA,
																	particleOptionsB = None,ratioA = 1, ratioB = 0):
	return coords,particles
def Maxwell_Boltzmann(coords, particles, kBT) :
	return Velocity
def LJ_Energy(coords, particles, Lbox) :
	return energy, forces
def VVerlet_NH_Step_1(timestep, particles, NHQ, kBT, KE, coords, velocities, forces, zeta, lns) :
	return lns_t_plus, zeta_t_half, coords_t_plus, velocities_t_half
def VVerlet_NH_Step_2(timestep, particles, NHQ, kBT, KE, velocities, forces_t_plus, zeta) :
	return zeta_t_plus, velocities_t_plus
def Kinetic_Energy(particles, velocity) :
	return Energy		
def Accumulate_RDF(coords, Lbox, histogram):
def Partial_RDF(coords, particles,Lbox, histogram,fromAtom, toAtom):
def Calculate_RDF(n, Lbox, histogram,fileh) :
def Calculate_Partial_RDF(n, Lbox, nA,nB,t,histogram,fileh) :
def MD_Loop(coords,velocity,particles,options,kBT) :
"""
