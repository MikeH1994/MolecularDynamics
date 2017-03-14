from constants import *

class ParticleOptions():
	def __init__(self,name = "Ar"):
		self._epsilon_k = 0
		self._epsilon = 0
		self._sigma = 0
		self._mass = 0
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
	 																								

class RunOptions():
	def __init__(self,PrimBox = 5.7064, timestep = 0.01e-12/tunit, 
			repeat = 3, XYZfile = "", nsteps = 20):
		self._PrimBox = PrimBox
		self._timestep = timestep
		self._repeat = repeat
		self._XYZfile = XYZfile
		self._nsteps = nsteps
		self._Lbox = PrimBox*repeat
		print(('RunOptions: PrimBox = {}, TimeStep = {}, Repeat = {}, '+
							'XYZfile = {}, nSteps = {}, Lbox = {}').format(PrimBox,timestep,
							repeat,XYZfile,nsteps,self._Lbox))
		
class Particle():
	def __init__(self,mass,sigma,epsilon,name = ""):
		self._mass = mass
		self._sigma = sigma
		self._epsilon = epsilon
		self._name = name
