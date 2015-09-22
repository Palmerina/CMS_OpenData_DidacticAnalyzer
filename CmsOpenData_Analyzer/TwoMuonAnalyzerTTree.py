# Name: TwoMuonAnalyzer.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 


__author__ = "Palmerina Gonzalez Izquierdo"
__copyright__ = "Copyright (C) 2015 Palmerina G. I."
__license__ = "Public Domain"
__version__ = "1.0"
__maintainer__ = "Palmerina Gonzalez"
__email__ = "pgi25@alumnos.unican.es"

import ROOT
from LeptonPair import LeptonPair #class LeptonPair inside LeptonPair.py
from CutsConfig import CutsConfig
import numpy as n
from scipy.stats import norm
#from scipy import optimize
#import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('QT4agg')
import matplotlib.pylab as P
import array

class TwoMuonAnalyzer(object):
	"""
	Analyzes the properties of the muons in every event
	and selects those coming from the Z boson decay
	"""
	
	def __init__(self, cutsConfig):
		"""
		TwoMuonAnalyzer initializer
		"""
		self.f = ROOT.TFile("mytree.root", "read")
		self.tree = self.f.Get("test")

		self.zMass = []
		self.badZMass = []
		self.Pt = array.array("d", [])
		self.badZPt = []
		self.zPt1 = []
		self.badZPt1 = []
		self.zPt2 = []
		self.badZPt2 = []
		self.eta = []
		self.badEta = []
		self.chi2 = []
		self.badChi2 = []
		self.numValidHits = []
		self.badNumValidHits = []
		self.dB = []
		self.distance = []


		

	def getMuons(self):
		"""
		event: one element of self.events
		
		returns:
		"""


		Muon_pt = array.array("d", [0.])
		Muon_eta = array.array("d", [0.])
		Muon_px = array.array("d", [0.])
		Muon_py = array.array("d", [0.])
		Muon_energy = array.array("d", [0.])
		Muon_vertex_z = array.array("d", [0.])
		Muon_isGlobalMuon = array.array("i", [0])
		Muon_isTrackerMuon = array.array("i", [0])
		Muon_dB = array.array("d", [0.])
		Muon_edB = array.array("d", [0.])
		Muon_isolation_sumPt = array.array("d", [0.])
		Muon_isolation_emEt = array.array("d", [0.])
		Muon_isolation_hadEt = array.array("d", [0.])
		Muon_numberOfValidHits = array.array("i", [0])
		Muon_normChi2 = array.array("d", [0.])
		Muon_charge = array.array("d", [0.])

		self.tree.SetBranchAddress("Muon_pt", Muon_pt)
		self.tree.SetBranchAddress("Muon_eta", Muon_eta)
		self.tree.SetBranchAddress("Muon_px", Muon_px)
		self.tree.SetBranchAddress("Muon_py", Muon_py)
		self.tree.SetBranchAddress("Muon_energy", Muon_energy)
		self.tree.SetBranchAddress("Muon_vertex_z", Muon_vertex_z)
		self.tree.SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon)
		self.tree.SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon)
		self.tree.SetBranchAddress("Muon_dB", Muon_dB)
		self.tree.SetBranchAddress("Muon_isolation_sumPt", Muon_isolation_sumPt)
		self.tree.SetBranchAddress("Muon_isolation_emEt", Muon_isolation_emEt)
		self.tree.SetBranchAddress("Muon_isolation_hadEt", Muon_isolation_hadEt)
		self.tree.SetBranchAddress("Muon_numberOfValidHits", Muon_numberOfValidHits)
		self.tree.SetBranchAddress("Muon_normChi2", Muon_normChi2)
		self.tree.SetBranchAddress("Muon_charge", Muon_charge)



	def getVertex(self):
		"""
		event: one element of self.events
		
		returns:
		"""
		
		Vertex_z = array.array("d", [0.])

		self.tree.SetBranchAddress("Vertex_z", Vertex_z)


	def selectMuons(self):
		"""
		muon:
		vertex:
	
		returns: boolean
		"""

	        #muon=getMuons(), vertex=getVertex()
		#The muon must be detected by both the tracker and the muon chambers
		if not (Muon_isGlobalMuon and Muon_isTrackerMuon):

			return False
	
		# Minimum transverse momentum (pt) and maximum eta angle
		if Muon_pt < self.cutsConfig.pt_min or abs(Muon_eta) > self.cutsConfig.eta_max:
			return False

		# Maximum distance of the muon respect to the vertex
		if abs(Muon_vertex_z - Vertex_z) > self.cutsConfig.distance:
	       		return False

	       	# Maximum impact parameter
	       	if Muon_dB > self.cutsConfig.dB_max:
	       		return False


	       	# I_trk + I_ECAL + I_HCAL
	       	# sumPt = suma de los momentos transversos
	       	# emEt = electromagnetic energy
	       	# hadEt = hadronic energy
		# Maximum energy content in that region before consider the "muon" as a jet of particle
		if (Muon_isolation_sumPt + Muon_isolation_emEt + Muon_isolation_hadEt) / Muon_pt > self.cutsConfig.isolation:
			return False

	       	# muon SIP variable # Symmetrized Impact Parameter in 2010?
	       	if (Muon_dB / Muon_edB) > 4:
	       		return False
	    
		if not Muon_normChi2 == 0.0:
	       		# Maximum chi2
	       		if Muon_normChi2 > self.cutsConfig.chi2:
	       			return False
	
		if not Muon_numberOfValidHits == 0:

	       		# Minimum number of valid hits on the global track. 
	       		if Muon_numberOfValidHits < self.cutsConfig.numValidHits:
	       			return False

	       	return True


	def plotter1(self):
		"""
		Plots the transverse momentum 

		"""

	#	P.figure()
	#	P.hist(self.badZPt, bins = 100, normed=1, alpha=0.5)
	#	P.xlim(0, 500)
	#	P.xlabel("Total pt (GeV/c)")
	#	P.ylabel("frequency")

		fig1 = P.figure()
		ax_1=fig1.add_subplot(211)
		ax_1.hist(self.badZPt1, bins = 100, alpha=0.5)
		ax_1.set_xlim(0, 200)
		ax_1.set_xlabel("pt_1 (GeV/c)")
		ax_1.set_ylabel("frequency")
		ax_1.set_title("Transverse momentum of each component of the muon pair")

		ax_2=fig1.add_subplot(212)
		ax_2.hist(self.badZPt2, bins = 300, alpha=0.5)
		ax_2.set_xlim(0, 1000)
		ax_2.set_xlabel("pt_2 (GeV/c)")
		ax_2.set_ylabel("frequency")

		fig2 = P.figure()
		ax_3 = fig2.add_subplot(211)
		ax_3.hist(self.badChi2, bins = 100, alpha=0.5)
		ax_3.set_xlabel("Chi**2")
		ax_3.set_ylabel("frequency")


		ax_4 = fig2.add_subplot(212)
		ax_4.hist(self.badNumValidHits, bins = 20, alpha=0.5)
		ax_4.set_xlabel("Number of valid hits")
		ax_4.set_ylabel("frequency")

		fig3 = P.figure()
		ax_5 = fig3.add_subplot(211)
		ax_5.hist(self.dB, bins = 100, alpha=0.5, log=True)
		ax_5.set_xlabel("Impact parameter")
		ax_5.set_ylabel("frequency")
		ax_5.set_title("Distance to the primary vertex")


		ax_6 = fig3.add_subplot(212)
		ax_6.hist(self.distance, bins = 100, alpha=0.5, log=True)
		ax_6.set_xlabel("Dz to PV")
		ax_6.set_ylabel("frequency")

		P.show()


	def plotter(self):

		"""
		Plots the histograms
		"""
		
		P.figure()
		P.hist(self.Pt, bins = 50, log = True)
		P.xlabel("pt")
		P.ylabel("frequency")


		P.show()


		


	def process(self, maxEv = -1):
		"""
		maxEv: maximum number of processed events
		       maxEv=-1 runs over all the events

		It selects the good muons applying the cut configuration
		and paires up them creating objects of the class LeptonPair.
		It gets the mass of every pair and adds the one which approaches 
		the most to the Z boson's mass to the list self.zMass.  

		"""
		Muon_pt = array.array("d", [0.])
		Muon_eta = array.array("d", [0.])
		Muon_px = array.array("d", [0.])
		Muon_py = array.array("d", [0.])
		Muon_energy = array.array("d", [0.])
		Muon_vertex_z = array.array("d", [0.])
		Muon_isGlobalMuon = array.array("i", [0])
		Muon_isTrackerMuon = array.array("i", [0])
		Muon_dB = array.array("d", [0.])
		Muon_edB = array.array("d", [0.])
		Muon_isolation_sumPt = array.array("d", [0.])
		Muon_isolation_emEt = array.array("d", [0.])
		Muon_isolation_hadEt = array.array("d", [0.])
		Muon_numberOfValidHits = array.array("i", [0])
		Muon_normChi2 = array.array("d", [0.])
		Muon_charge = array.array("d", [0.])

		self.tree.SetBranchAddress("Muon_pt", Muon_pt)
		self.tree.SetBranchAddress("Muon_eta", Muon_eta)
		self.tree.SetBranchAddress("Muon_px", Muon_px)
		self.tree.SetBranchAddress("Muon_py", Muon_py)
		self.tree.SetBranchAddress("Muon_energy", Muon_energy)
		self.tree.SetBranchAddress("Muon_vertex_z", Muon_vertex_z)
		self.tree.SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon)
		self.tree.SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon)
		self.tree.SetBranchAddress("Muon_dB", Muon_dB)
		self.tree.SetBranchAddress("Muon_isolation_sumPt", Muon_isolation_sumPt)
		self.tree.SetBranchAddress("Muon_isolation_emEt", Muon_isolation_emEt)
		self.tree.SetBranchAddress("Muon_isolation_hadEt", Muon_isolation_hadEt)
		self.tree.SetBranchAddress("Muon_numberOfValidHits", Muon_numberOfValidHits)
		self.tree.SetBranchAddress("Muon_normChi2", Muon_normChi2)
		self.tree.SetBranchAddress("Muon_charge", Muon_charge)

		
		Vertex_z = array.array("d", [0.])
		self.tree.SetBranchAddress("Vertex_z", Vertex_z)

		h_Muon_pt = ROOT.TH1D("Muon_pt", "Transverse momentum", 100, 0, 1000)
	

		for i in range(0, self.tree.GetEntries()):


			selectedMuons = []


			
			self.tree.GetEntry(i)
		
			print Muon_pt
	
	#		h_Muon_pt.Fill(Muon_pt[0])
			
			self.Pt.append(Muon_pt[0])		

	#	h_Muon_pt.Draw()
