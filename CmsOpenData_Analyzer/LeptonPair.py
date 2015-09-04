# Name: LeptonPair.py
#
# CMS Open Data
#
# Description: sums pairs of leptons and gets their mass, energy, momentum and transverse momentum
#
# Returns: 


__author__ = "Palmerina Gonzalez Izquierdo"
__copyright__ = "Copyright (C) 2015 Palmerina G. I."
__license__ = "Public Domain"
__version__ = "1.0"
__maintainer__ = "Palmerina Gonzalez"
__email__ = "pgi25@alumnos.unican.es"

import ROOT
import math
import numpy as np

class LeptonPair(object):

	def __init__(self, l1, l2):
	 	"""
	 	l1, l2: leptons (getMuons() from GetData.py) or l1 = Z boson and l2 = None
	 	It sums the four-momentums of l1 and l2 and gets their mass, their energy and their transverse momentum
	 	"""

 		self.px = l1.px()+l2.px() #GeV/c
 		self.py = l1.py()+l2.py()
 		self.pz = l1.pz()+l2.pz()
		self.energy = l1.energy() + l2.energy() #GeV

		self.squareEnergy = np.power(self.energy,2)			
	 	self.squareP = np.power(self.px,2) + np.power(self.py,2) + np.power(self.pz,2)

	 	#l.px(),l.py(),l.pz() and l.energy() are ROOT methods


	def mass(self):
	 	#invariant mass
		squareM = self.squareEnergy-self.squareP
		if squareM < 0:
			return np.sqrt(-squareM) # con un - en el codigo de TLorentzVector
		else:
			return np.sqrt(squareM) # in natural units (GeV/c**2)
		#return self.p4.M()
		

	def pt(self):
	 	#Transverse momentum
	 	return np.sqrt(np.power(self.px,2) + np.power(self.py,2))
