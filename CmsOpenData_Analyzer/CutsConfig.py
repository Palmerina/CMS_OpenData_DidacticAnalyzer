# Name: LCutsConfig.py
#
# CMS Open Data
#
# Description: Cuts configuration to select the good muons.
#
# Returns: 


__author__ = "Palmerina Gonzalez Izquierdo"
__copyright__ = "Copyright (C) 2015 Palmerina G. I."
__license__ = "Public Domain"
__version__ = "1.0"
__maintainer__ = "Palmerina Gonzalez"
__email__ = "pgi25@alumnos.unican.es"



class CutsConfig(object):
	"""

	"""

	def __init__(self, pt_min = 5, eta_max = 2.4, distance = 0.2, dB_min = 0.02, isolation = 0.15):
		
		self.pt_min = pt_min 
		self.eta_max = eta_max
		self.distance = distance
		self.dB_min = dB_min # cm. dB=impact parameter
		#normChi2_max = 10
		self.isolation = isolation #dimensionless. (sumPt+emEnergy+hadEnergy)/muon.pt = maxima energia antes de considerarlo como un jet de particulas.
<<<<<<< HEAD
		#SIP variable
=======
		#SIP variable?
>>>>>>> 99aa6b109d672ca3c4df52eb5a73dd5f4596431c
