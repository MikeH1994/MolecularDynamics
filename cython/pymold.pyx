import math
import random
import time
import sys
import getopt
import copy
import os
from cpython cimport array

kB    = 1.38064852e-23        # Boltzmann's constant
NA    = 6.022140857e23        # Avogadro number
e     = 1.6021766208e-19      # Fundamental Charge
evk   = e/kB                  # Conversion factor from eV to K
kjmol = e*NA/1e3              # Conversion factor from eV to kJ/mol
tunit = 1.01805058e-14        # Unit of time in a.m.u, eV, A units system
PrimBox = 6.101
#5.7064
#
arMass = 39.948			   # Mass of Ar atom in a.m.u
arEpsilon_k = 119.8
arSigma = 3.405
arFcc = 5.7064

krMass = 83.798
krEpsilon_k = 167.0
krSigma = 3.633

neMass = 20.1797
neEpsilon_k = 36.2
neSigma = 2.800

primfrac = [ 
			[0.0,  0.0,  0.0],
			[0.0,  0.5,  0.5],
			[0.5,  0.0,  0.5],
			[0.5,  0.5,  0.0]]

particleAIndex = []
particleBIndex = []

cdef class ParticleOptions:
	cdef public float _epsilon_k
	cdef public float _epsilon
	cdef public float _sigma
	cdef public float _mass
	cdef public str _name
	def __init__(self,name = "Ar"):
		self._name = name
		if self._name == "Ar":
			self._sigma = arSigma
			self._epsilon_k = arEpsilon_k
			self._mass = arMass
		elif self._name == "Kr":
			self._sigma = krSigma
			self._epsilon_k = krEpsilon_k
			self._mass = krMass
		else:
			self._name = "Ne"
			self._sigma = neSigma
			self._epsilon_k = neEpsilon_k
			self._mass = neMass
		
		self._epsilon= self._epsilon_k*kB/e		 # LJ epsilon for Argon in units of eV
		print('Particle {} Options: sigma = {}, epsilon_k = {}, mass = {}'.format(
	 		self._name,self._sigma,self._epsilon_k,self._mass))
	 																								

cdef class RunOptions:
	cdef public float _PrimBox
	cdef public float _timestep
	cdef public int _repeat
	cdef public str _XYZfile
	cdef public int _nsteps
	cdef public float _Lbox
	cdef public float _ratioA
	cdef public float _ratioB
	cdef public str _particleAName
	cdef public str _particleBName
	cdef public float _NHQ
	def __init__(self,PrimBox = 5.7064, timestep = 0.01e-12/tunit, 
			repeat = 3, XYZfile = "", nsteps = 20, ratioA = 1,ratioB = 0,
			partAName = "", partBName = "", NHQ = 16*arMass):
		self._PrimBox = PrimBox
		self._timestep = timestep
		self._repeat = repeat
		self._XYZfile = XYZfile
		self._nsteps = nsteps
		self._Lbox = PrimBox*repeat
		self._ratioA = ratioA
		self._ratioB = ratioB
		self._particleAName = partAName
		self._particleBName = partBName
		self._NHQ    = NHQ
		print(('RunOptions: PrimBox = {}, TimeStep = {}, Repeat = {}, '+
							'XYZfile = {}, nSteps = {}, Lbox = {}').format(self._PrimBox,self._timestep,
							self._repeat,self._XYZfile,self._nsteps,self._Lbox))
		
cdef class Particle:
	cdef public float _mass
	cdef public float _sigma
	cdef public float _epsilon
	cdef public str _name
	def __init__(self,mass,sigma,epsilon,name = ""):
		self._mass = mass
		self._sigma = sigma
		self._epsilon = epsilon
		self._name = name

def GenerateSupercell(runOptions,particleOptionsA,
																	particleOptionsB = None,ratioA = 1, ratioB = 0, uniformDistribution = False):
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
	totA = float (ratioA)/(ratioA+ratioB)*totParticles
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
					if uniformDistribution:
						#if uniformDistribution, switch between A and B every swapEvery times
						if swapIndex>=swapEvery and not ratioB==0:
							typeA = not typeA
							swapIndex = 0
							if typeA or ratioB==0:
								swapEvery = ratioA
							else:
								swapEvery = ratioB
					else:
						#if not, 
						if n>=totA:
							typeA = False
					if typeA:
						name = particleOptionsA._name
						mass = particleOptionsA._mass
						sigma = particleOptionsA._sigma
						epsilon = particleOptionsA._epsilon
						particleAIndex.append(n)
					else:
						name = particleOptionsB._name
						mass = particleOptionsB._mass
						sigma = particleOptionsB._sigma
						epsilon = particleOptionsB._epsilon
						particleBIndex.append(n)
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
cpdef LJ_Energy(list coords,list particles, Lbox) :
	n = len(coords)
	cdef float energy	= 0.0
	rij=[0,0,0]
	forces = [[0,0,0] for i in range(n) ]
	cdef float sigma = 0.0
	cdef float dist_sq = 0.0
	cdef float rcomp = 0.0
	cdef float force = 0.0
	cdef Particle particlei
	cdef Particle particlej
	for i in range(n) :
		particlei = particles[i]
		for j in range(i+1,n):
			particlej = particles[j]
			sigma = (particlei._sigma + particlej._sigma)/2.0
			epsilon = math.sqrt(particlei._epsilon*particlej._epsilon)
			dist_sq = 0.0
			for k in range(3) :
				rcomp = coords[j][k] - coords[i][k]
				rcomp -= Lbox*((rcomp/Lbox+0.5)//1)
				rij[k] = rcomp
				dist_sq += rcomp**2
			if( dist_sq <= Lbox**2/4.0 ) :
				rdist = (sigma**2/dist_sq)**3
				energy += 4.0*epsilon*(rdist**2 - rdist)
				force = 24.0*epsilon * (2.0*rdist**2 -rdist)/dist_sq
				for k in range(3) :
					forces[i][k] -= force*rij[k]
					forces[j][k] += force*rij[k]
				#print "Force=",force
	return energy, forces

#
# Step the atomic co-ordinates and Nose-Hoover variables using the Velocity Verlet algorithm
#
cpdef VVerlet_NH_Step_1(float timestep,list particles,float NHQ,float kBT,float KE,list coords,list velocities,list forces,float zeta,float lns) :
	n = len(coords)
	cdef list coords_t_plus = [ [0,0,0] for i in range(0, n)]
	cdef list velocities_t_half = [ [0,0,0] for i in range(0, n)]
	cdef Particle particlei
	cdef float coord
	cdef float veloc
	cdef float force
	for i in range(n) :
		particlei = particles[i]
		for k in range(3) :
			coord = coords[i][k]
			veloc = velocities[i][k]
			force = forces[i][k]
			coords_t_plus[i][k] = coord + timestep*veloc + timestep**2*(force - zeta*veloc)/(2.0*particlei._mass)
			velocities_t_half[i][k] = veloc + timestep*(force - zeta*veloc)/(2.0*particlei._mass)

	zeta_t_half = zeta + timestep/(2.0*NHQ)*(KE - (3*n)*kBT/2.0)
	lns_t_plus = lns + timestep*zeta_t_half + timestep**2/NHQ*(KE - (3*n)*kBT/2.0)
		
	return lns_t_plus, zeta_t_half, coords_t_plus, velocities_t_half
#
# Step the atomic co-ordinates and Nose-Hoover variables using the Velocity Verlet algorithm
# An approximate scheme is used for the second velocity half-step.
#
cdef VVerlet_NH_Step_2(float timestep, list particles,float NHQ,float kBT,float KE,list velocities,list forces_t_plus,float zeta) :
	n = len(velocities)
	cdef list velocities_t_plus = [ [0,0,0] for i in range(0, n)]
	cdef Particle particlei
	cdef float veloc
	cdef float force
	for i in range(n) :
		particlei = particles[i]
		for k in range(3) :
			veloc = velocities[i][k]
			force = forces_t_plus[i][k]
			velocities_t_plus[i][k] = (veloc + timestep*force/(2.0*particlei._mass)) \
										/ (1 + (timestep*zeta/2.0))


	zeta_t_plus = zeta + timestep/(2.0*NHQ)*(KE - (3*n)*kBT/2.0)

	return zeta_t_plus, velocities_t_plus
#
# Compute the instantaneous kinetic energy from the velocities
#
cdef float Kinetic_Energy(list particles,list velocity) :
	n = len(velocity)
	cdef float Energy = 0.0
	cdef Particle particlei
	for i in range(n) :
		particlei = particles[i]
		Energy += 0.5*particlei._mass*(velocity[i][0]**2 + velocity[i][1]**2 + velocity[i][2]**2)	
	return Energy
		
cpdef getVACF(list velocity_t,list velocity_t0):
	n = len(velocity_t0)
	cdef float magT = 0.0
	cdef float magT0 = 0.0
	cdef float vacf = 0.0
	cdef float vacfI = 0.0
	cdef float veloct = 0.0
	cdef float veloc0 = 0.0
	for i in range(n):
		magT = 0.0
		magT0 = 0.0
		vacfI = 0.0
		for k in range(3):
			veloct = velocity_t[i][k]
			veloct0 = velocity_t0[i][k]
			magT+=veloct*veloct
			magT0+=veloct0*veloct0
			vacfI+=veloct0*veloct
		magT = math.sqrt(magT)
		magT0 = math.sqrt(magT0)
		vacf += vacfI/(magT*magT0)
	vacf/=n
	return vacf
	
def getMSD(list coords_t,list coords_t0):
	cdef float msdAB = 0.0
	cdef float msdA = 0.0
	cdef float msdB = 0.0
	cdef float coordt = 0
	cdef float coordt0 = 0
	n = len(coords_t0)
	for i in range(n):
		for k in range(3):
			coordt = coords_t[i][k]
			coordt0 = coords_t0[i][k]
			msdAB+=(coordt0-coordt)**2
	msdAB/=n
	
	n = len(particleAIndex)
	for i in particleAIndex:
		for k in range(3):
			coordt = coords_t[i][k]
			coordt0 = coords_t0[i][k]
			msdA+=(coordt0-coordt)**2
	msdA/=n
	
	n = len(particleBIndex)
	for i in particleBIndex:
		for k in range(3):
			coordt = coords_t[i][k]
			coordt0 = coords_t0[i][k]
			msdB+=(coordt0-coordt)**2
	msdB/=n
	
	return msdAB,msdA,msdB
	
def Accumulate_RDF(list coords, Lbox, histogram):
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

cdef Partial_RDF(list coords,list particles,float Lbox,list histogram,fromAtom,toAtom):
	
	n = len(coords)
	nbins = len(histogram)
	cdef float rcomp = 0
	cdef float dist_sq
	particleFromList,particleToList = [],[]
	if fromAtom=="A":
		particleFromList = particleAIndex
	else:
		particleFromList = particleBIndex
		
	if toAtom=="A":
		particleToList = particleAIndex
	else:
		particleToList = particleBIndex

	startingFrom = 0
	boolean = (fromAtom==toAtom)
	cdef int fromIndex,toIndex	
	for i in range(len(particleFromList)) :
		fromIndex = particleFromList[i]
		startingFrom = 0
		if boolean:
			startingFrom = i+1
		for j in range(startingFrom,len(particleToList)):
			toIndex = particleToList[j]
			dist_sq = 0.0
			for k in range(3) :
				rcomp = coords[toIndex][k] - coords[fromIndex][k]
				rcomp = rcomp - Lbox*((rcomp/Lbox+0.5)//1)
				dist_sq += rcomp**2
			if( dist_sq < Lbox**2/4.0 ) :
				bin = int(math.sqrt(4.0*dist_sq/Lbox**2) * nbins)
				try:
					histogram[bin] += 1
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
		
def Calculate_Partial_RDF(n, options, nA,nB,t,nSteps,histogram,fileh) :
	if nB==0:
		nB=nA
	fileh = open(fileh,'w')
	nbins = len(histogram)-1
	Lbox = options._Lbox
	bin_width = Lbox/(2*nbins)
	totParticles = options._repeat**3 * len(primfrac)
	v = Lbox**3
	sum = 0
	for i in range(nbins) :
		sum += histogram[i]
	fileh.write( '# {0:6s}\t{1:9s}\t{2:9s}\n'.format("r","g(r)","delta"))
	for i in range(nbins) :
		r = (i+0.5)*Lbox/(2.0*nbins)
		rdf = 3*v*histogram[i] /(4.0*math.pi*nA*nB*t*((r+bin_width)**3-r**3))
		rdfError = 3*v*math.sqrt(histogram[i]) /(4.0*math.pi*nA*nB*t*((r+bin_width)**3-r**3))
		fileh.write( '{0:6f}\t{1:9.6f}\t{2:9.6f}\n'.format(r,rdf,rdfError))
	fileh.close()
		
def MD_Loop(coords,velocity,particles,options,kBT) :
	print '# {0:4s}\t{1:12s}\t{2:12s}\t{3:9s}\t{4:9s}\t{5:12s}\t{6:6s}\t{7:4s}\t{8:3s}\t{9:3s}'.format("Step",
           "Pot.Energy","Kin.Energy","NHzeta","s","Tot.Energy","T(Inst)","MSD_AB","MSD_A","MSD_B")
	n = len(coords)
	nbins = 150
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

	mixRatio = float (options._ratioA)/(options._ratioA + options._ratioB)
	totParticles = options._repeat**3 * len(primfrac)
	nA = mixRatio*totParticles
	nB = (1-mixRatio)*totParticles
	percentageA = '{0:3.0f}'.format(nA/(nA+nB)*100)
	percentageB = '{0:3.0f}'.format(nB/(nA+nB)*100)
	folderName = '{}{}_{}_{}_{}_{}'.format(nameA,nameB,percentageA,percentageB,totParticles,nsteps)
	if not os.path.exists("../" + folderName):
		os.makedirs("../" + folderName)

	###############################################
	XYZfh = None
	if XYZfile != "":
		XYZfile = "../" + folderName + "/" + XYZfile
		XYZfh = open(XYZfile,'w')
		XYZfh.write('# {0:4s}\t{1:12s}\t{2:12s}\t{3:9s}\t{4:9s}\t{5:12s}\t{6:6s}\t{7:4s}\t{8:3s}\t{9:3s}\n'.format("Step",
           "Pot.Energy","Kin.Energy","NHzeta","s","Tot.Energy","T(Inst)","MSD_AB","MSD_A","MSD_B"))
	###############################################	
	zeta = 0
	lns = 0
	
	coords_t0 = copy.deepcopy(coords)
	velocity_t0 = copy.deepcopy(velocity)
	Pot_Energy, forces = LJ_Energy(coords, particles, Lbox)
	Kin_Energy = Kinetic_Energy(particles, velocity)
	D = 0
	msd = 0
	for istep in range(nsteps) :
		lns, zeta, coords_t_plus, velocity_t_plus =	VVerlet_NH_Step_1(timestep, particles, NHQ, kBT, Kin_Energy,coords, velocity, forces, zeta, lns)
		Pot_Energy, forces_t_plus = LJ_Energy(coords_t_plus,particles, Lbox)
		Kin_Energy = Kinetic_Energy(particles, velocity_t_plus)
		zeta, velocity_t_plus = VVerlet_NH_Step_2(timestep, particles, NHQ, kBT, Kin_Energy,velocity_t_plus, forces_t_plus, zeta )
		T_Inst = evk* 2.0*Kin_Energy/(3*n - 3)
		forces = forces_t_plus
		velocity = velocity_t_plus
		coords = coords_t_plus
		Tot_Energy = Pot_Energy + Kin_Energy + 0.5*NHQ*zeta**2 + (3*n)*kBT*math.exp(lns)
		
		if istep%5==0:
			Partial_RDF(coords, particles,Lbox, gr_aa_histogram,"A", "A")
			Partial_RDF(coords, particles,Lbox, gr_ab_histogram,"A", "B")
			Partial_RDF(coords, particles,Lbox, gr_bb_histogram,"B", "B")
			msdAB,msdA,msdB = getMSD(coords,coords_t0)
			print '{0:4d}\t{1:12.7f}\t{2:12.7f}\t{3:9.6f}\t{4:9.6f}\t{5:12.7f}\t{6:6.2f}\t{7:1.6f}\t{8:4.6f}\t{9:4.6f}'.format(istep,
			Pot_Energy,Kin_Energy, zeta, math.exp(lns), Tot_Energy,T_Inst,msdAB,msdA,msdB)
			if XYZfile != "":
				XYZfh.write( '{0:4d}\t{1:12.7f}\t{2:12.7f}\t{3:9.6f}\t{4:9.6f}\t{5:12.7f}\t{6:6.2f}\t{7:4.6f}\t{8:4.6f}\t{9:4.6f}\n'.format(istep,
			Pot_Energy,Kin_Energy, zeta, math.exp(lns), Tot_Energy,T_Inst,msdAB,msdA,msdB))
	
	if XYZfile != "":
		XYZfh.close()
	Calculate_Partial_RDF(n, options, nA,nA,nsteps*timestep,nsteps,gr_aa_histogram,'../{}/grAA.dat'.format(folderName))
	Calculate_Partial_RDF(n, options, nA,nB,nsteps*timestep,nsteps,gr_ab_histogram,'../{}/grAB.dat'.format(folderName))
	Calculate_Partial_RDF(n, options, nB,nB,nsteps*timestep,nsteps,gr_bb_histogram,'../{}/grBB.dat'.format(folderName))
	

def run():
	###############################################################
	timestep = 0.01e-12/tunit
	initial_temperature = 115.8
	repeat = 4
	XYZfile = "data.dat"
	nsteps = 50000
	###############################################################
	ratioA = 1
	ratioB = 1
	try:
		opts, args = getopt.getopt(sys.argv[1:],"ht:n:d:r:x:a:b:b",["help","temperature=","nsteps=","timestep=","box-repeat=","xyz-file","ratioA=","ratioB="])
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
		elif opt in ("-a", "-ra"):
			ratioA = int (arg)
		elif opt in ("-b", "-rb"):
			ratioB = int (arg)
	###############################################################
	particleOptionsA = ParticleOptions("Ar")
	particleOptionsB = ParticleOptions("Kr")
	###############################################################
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

run()
