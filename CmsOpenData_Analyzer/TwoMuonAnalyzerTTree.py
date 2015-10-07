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
		self.tree = self.f.Get("muons")

		self.cutsConfig = cutsConfig
		

		# Initialization of variables		
		
		self.Muon_pt = array.array("f")
		self.Muon_eta = array.array("f", [0.]*50)
		self.Muon_px = array.array("f", [0.]*50)
		self.Muon_py = array.array("f", [0.]*50)
		self.Muon_pz = array.array("f", [0.]*50)
		self.Muon_energy = array.array("f", [0.]*50)
		self.Muon_vertex_z = array.array("f", [0.]*50)
		self.Muon_isGlobalMuon = array.array("i", [0]*50)
		self.Muon_isTrackerMuon = array.array("i", [0]*50)
		self.Muon_dB = array.array("f", [0.]*50)
		self.Muon_edB = array.array("f", [0.]*50)
		self.Muon_isolation_sumPt = array.array("f", [0.]*50)
		self.Muon_isolation_emEt = array.array("f", [0.]*50)
		self.Muon_isolation_hadEt = array.array("f", [0.]*50)
		self.Muon_numberOfValidHits = array.array("i", [0]*50)
		self.Muon_normChi2 = array.array("f", [0.]*50)
		self.Muon_charge = array.array("i", [0]*50)

		self.Vertex_z = array.array("f", [0.])
		self.npart = array.array("i", [0])

		


		# Arrays where the variables are going to be stored
	
		self.allMuons_pt = array.array("f")
		self.allMuons_eta = array.array("f", [0.])	
		self.allMuons_px = array.array("f", [0.])
		self.allMuons_py = array.array("f", [0.])
		self.allMuons_pz = array.array("f", [0.])
		self.allMuons_energy = array.array("f", [0.])
		self.allMuons_isGlobalMuon = array.array("i", [0])
		self.allMuons_isTrackerMuon = array.array("i", [0])
		self.allMuons_vertex_z = array.array("f", [0.])
		self.allMuons_dB = array.array("f", [0.])
		self.allMuons_edB = array.array("f", [0.])
		self.allMuons_isolation_sumPt = array.array("f", [0.])
		self.allMuons_isolation_emEt = array.array("f", [0.])
		self.allMuons_isolation_hadEt = array.array("f", [0.])
		self.allMuons_numberOfValidHits = array.array("i", [0])
		self.allMuons_normChi2 = array.array("f", [0.])
		self.allMuons_charge = array.array("i", [0])
		self.mass = array.array("f", [0.])

		self.event_vertex_z = array.array("f", [0.])

		self.goodMuons_pt = array.array("f")
		self.goodMuons_eta = array.array("f", [0.])	
		self.goodMuons_px = array.array("f", [0.])
		self.goodMuons_py = array.array("f", [0.])
		self.goodMuons_pz = array.array("f", [0.])
		self.goodMuons_energy = array.array("f", [0.])
		self.goodMuons_isGlobalMuon = array.array("i", [0])
		self.goodMuons_isTrackerMuon = array.array("i", [0])
		self.goodMuons_vertex_z = array.array("f", [0.])
		self.goodMuons_dB = array.array("f", [0.])
		self.goodMuons_edB = array.array("f", [0.])
		self.goodMuons_isolation_sumPt = array.array("f", [0.])
		self.goodMuons_isolation_emEt = array.array("f", [0.])
		self.goodMuons_isolation_hadEt = array.array("f", [0.])
		self.goodMuons_numberOfValidHits = array.array("i", [0])
		self.goodMuons_normChi2 = array.array("f", [0.])
		self.goodMuons_charge = array.array("i", [0])
		self.good_mass = array.array("f", [0.])

	def selectMuons(self, iMuon):
		"""
		muon:
		vertex:
	
		returns: boolean
		"""

		

	        #muon=getMuons(), vertex=getVertex()
		#The muon must be detected by both the tracker and the muon chambers
		if not (self.Muon_isGlobalMuon[iMuon] and self.Muon_isTrackerMuon[iMuon]):

			return False
	
		# Minimum transverse momentum (pt) and maximum eta angle
		if self.Muon_pt[iMuon] < self.cutsConfig.pt_min or abs(self.Muon_eta[iMuon]) > self.cutsConfig.eta_max:
			return False

		# Maximum distance of the muon respect to the vertex
		if abs(self.Muon_vertex_z[iMuon] - self.Vertex_z[0]) > self.cutsConfig.distance:
	       		return False

	       	# Maximum impact parameter
	       	if self.Muon_dB[iMuon] > self.cutsConfig.dB_max:
	       		return False


	       	# I_trk + I_ECAL + I_HCAL
	       	# sumPt = suma de los momentos transversos
	       	# emEt = electromagnetic energy
	       	# hadEt = hadronic energy
		# Maximum energy content in that region before consider the "muon" as a jet of particle
		if (self.Muon_isolation_sumPt[iMuon] + self.Muon_isolation_emEt[iMuon] + self.Muon_isolation_hadEt[iMuon]) / self.Muon_pt[iMuon] > self.cutsConfig.isolation:
			return False

	       	# muon SIP variable # Symmetrized Impact Parameter in 2010?
	       	#if (self.Muon_dB[0] / self.Muon_edB[0]) > 4:
	       	#	return False
	    
		if not self.Muon_normChi2[iMuon] == 0.0:
	       		# Maximum chi2
	       		if self.Muon_normChi2[iMuon] > self.cutsConfig.chi2:
	       			return False
	
		if not self.Muon_numberOfValidHits[iMuon] == 0:

	       		# Minimum number of valid hits on the global track. 
	       		if self.Muon_numberOfValidHits[iMuon] < self.cutsConfig.numValidHits:
	       			return False

	       	return True


	


		


	def process(self):
		"""
		maxEv: maximum number of processed events
		       maxEv=-1 runs over good the events
		It selects the good muons applying the cut configuration
		and paires up them creating objects of the class LeptonPair.
		It gets the mass of every pair and adds the one which approaches 
		the most to the Z boson's mass to the list self.zMass.  
		"""
	
		print "arranca process"	
		# Address arrays to the TTree's branches

		self.tree.SetBranchAddress("Muon_pt", self.Muon_pt)
		self.tree.SetBranchAddress("Muon_eta", self.Muon_eta)
		self.tree.SetBranchAddress("Muon_px", self.Muon_px)
		self.tree.SetBranchAddress("Muon_py", self.Muon_py)
		self.tree.SetBranchAddress("Muon_pz", self.Muon_pz)
		self.tree.SetBranchAddress("Muon_energy", self.Muon_energy)
		self.tree.SetBranchAddress("Muon_vertex_z", self.Muon_vertex_z)
		self.tree.SetBranchAddress("Muon_isGlobalMuon", self.Muon_isGlobalMuon)
		self.tree.SetBranchAddress("Muon_isTrackerMuon", self.Muon_isTrackerMuon)
		self.tree.SetBranchAddress("Muon_dB", self.Muon_dB)
		self.tree.SetBranchAddress("Muon_isolation_sumPt", self.Muon_isolation_sumPt)
		self.tree.SetBranchAddress("Muon_isolation_emEt", self.Muon_isolation_emEt)
		self.tree.SetBranchAddress("Muon_isolation_hadEt", self.Muon_isolation_hadEt)
		self.tree.SetBranchAddress("Muon_numberOfValidHits", self.Muon_numberOfValidHits)
		self.tree.SetBranchAddress("Muon_normChi2", self.Muon_normChi2)
		self.tree.SetBranchAddress("Muon_charge", self.Muon_charge)

		self.tree.SetBranchAddress("Vertex_z", self.Vertex_z)
		self.tree.SetBranchAddress("npart", self.npart)

		numEntries = self.tree.GetEntries()


		# Loop over the events

		for i in range(0, numEntries): 
			self.tree.GetEntry(i)  # Muon_* arrays are filled for each event
			print self.Muon_pt
			# Select events with at least two muons
			if self.npart[0]<2:
				continue
			
			print self.npart[0]

			# ALL MUONS
			# Loop over the muons
			for iMuon in range(0, self.npart[0] - 1):  
              
				# Histograms to plot             
				self.allMuons_pt.append(self.Muon_pt[iMuon])
				self.allMuons_eta.append(self.Muon_eta[iMuon])	
				self.allMuons_energy.append(self.Muon_energy[iMuon])
				self.allMuons_vertex_z.append(self.Muon_vertex_z[iMuon])
				self.allMuons_dB.append(self.Muon_dB[iMuon])
				self.allMuons_edB.append(self.Muon_edB[iMuon])

				if not self.Muon_numberOfValidHits[iMuon] == 0:
					self.allMuons_numberOfValidHits.append(self.Muon_numberOfValidHits[iMuon])
				if not self.Muon_normChi2[iMuon] == 0.0:
					self.allMuons_normChi2.append(self.Muon_normChi2[iMuon])

				# Muon's four-momentum 
				outerMuon = ROOT.TLorentzVector(self.Muon_px[iMuon], self.Muon_py[iMuon], self.Muon_pz[iMuon], self.Muon_energy[iMuon])
				outerMuon_charge = self.Muon_charge[iMuon]

				# Selec the good muons
				if self.selectMuons(iMuon):

					# Histograms to plot
					self.goodMuons_pt.append(self.Muon_pt[iMuon])
					self.goodMuons_eta.append(self.Muon_eta[iMuon])	
					self.goodMuons_energy.append(self.Muon_energy[iMuon])
					self.goodMuons_vertex_z.append(self.Muon_vertex_z[iMuon])
					self.goodMuons_dB.append(self.Muon_dB[iMuon])
					self.goodMuons_edB.append(self.Muon_edB[iMuon])

					# GOOD MUONS
					# Loop over the muons
					for jMuon in range(0, self.npart[0] - 1):

						innerGoodMuon_charge = self.Muon_charge[jMuon]

						if innerGoodMuon_charge * outerMuon_charge >=0:
							continue

						# Selec the good muons
						if self.selectMuons(jMuon):

							innerGoodMuon = ROOT.TLorentzVector(self.Muon_px[jMuon], self.Muon_py[jMuon], self.Muon_pz[jMuon], self.Muon_energy[jMuon])

							goodMass = (outerMuon + innerGoodMuon).M()

							if not (goodMass > self.cutsConfig.mass_min and goodMass < 120):
								continue

							self.z_mass.append(goodMass)


				#ALL MUONS
				# Loop over all muons				
				for kMuon in range(0, self.npart[0] - 1):

					innerMuon_charge = self.Muon_charge[kMuon]

					if innerMuon_charge * outerMuon_charge >=0:
						continue

					innerMuon = ROOT.TLorentzVector(self.Muon_px[kMuon], self.Muon_py[kMuon], self.Muon_pz[kMuon], self.Muon_energy[kMuon])

					mass = (outerMuon + innerMuon).M()

					if not (mass > self.cutsConfig.mass_min and mass < 120):
						continue

					self.mass.append(mass)




	
		
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
		P.hist(self.allMuons_pt, bins = 50, log = True)
		P.xlabel("pt")
		P.ylabel("frequency")


		P.show()


		


	def plotHistos(self):


	#	c1 = ROOT.TCanvas("pt", "All muons' transverse momentum")
	#	c2 = ROOT.TCanvas("good pt", "Good muons' transverse momentum")

	#	self.h_Muon_pt.Draw()	
	#	self.h_goodMuon_pt.Draw()


		"""
		Plots the histograms
		"""
		
		P.figure()
		P.hist(self.h_Muon_pt, bins = 50, log = True)
		P.xlabel("pt")
		P.ylabel("frequency")

		P.figure()
		P.hist(self.h_goodMuon_pt, bins = 50, log = True)
		P.xlabel("pt")
		P.ylabel("frequency")

		P.show()


