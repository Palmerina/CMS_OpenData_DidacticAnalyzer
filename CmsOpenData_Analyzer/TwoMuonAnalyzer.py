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
from DataFormats.FWLite import Events, Handle
from CutsConfig import CutsConfig
import numpy as n
from scipy.stats import norm
import matplotlib.mlab as mlab
import matplotlib
matplotlib.use('QT4agg')
import matplotlib.pylab as P


class TwoMuonAnalyzer(object):
	"""
	Analyzes the properties of the muons in every event
	and selects those coming from the Z boson decay
	"""
	
	def __init__(self, cutsConfig, data_files):
		"""
		TwoMuonAnalyzer initializer
		"""

		self.muonHandle = Handle('std::vector<pat::Muon>')
		self.vertexHandle = Handle('std::vector<reco::Vertex>')	
		self.cutsConfig = cutsConfig
		self.events = Events(data_files)
		self.zMass = []
		self.badZMass = []
		self.zPt = []
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


	def getMuons(self, event):
		"""
		event: one element of self.events
		
		returns:
		"""

		event.getByLabel('patMuons', self.muonHandle)
		muons = self.muonHandle.product()
		return muons

	def getVertex(self, event):
		"""
		event: one element of self.events
		
		returns:
		"""

		event.getByLabel('offlinePrimaryVertices', self.vertexHandle)
		vertex = self.vertexHandle.product()[0] #it only takes the first element which corresponds to the primary vertex
		return vertex


	def selectMuons(self, muon, vertex):
		"""
		muon:
		vertex:
		
		returns: boolean
		"""

	        #muon=getMuons(), vertex=getVertex()
		#The muon must be detected by both the tracker and the muon chambers
		if not (muon.isGlobalMuon() and muon.isTrackerMuon()):
			return False
	
		# Minimum transverse momentum (pt) and maximum eta angle
		if muon.pt() < self.cutsConfig.pt_min or abs(muon.eta()) > self.cutsConfig.eta_max:
			return False

		# Maximum distance of the muon respect to the vertex
		if abs(muon.vertex().z() - vertex.z()) > self.cutsConfig.distance:
	       		return False

	       	# Maximum impact parameter
	       	if muon.dB(muon.PV3D) > self.cutsConfig.dB_min:
	       		return False


	       	# I_trk + I_ECAL + I_HCAL
	       	# sumPt = suma de los momentos transversos
	       	# emEt = electromagnetic energy
	       	# hadEt = hadronic energy
		# Maximum energy content in that region before consider the "muon" as a jet of particle
		if (muon.isolationR03().sumPt + muon.isolationR03().emEt + muon.isolationR03().hadEt) / muon.pt() > self.cutsConfig.isolation:
			return False

	       	# muon SIP variable # Symmetrized Impact Parameter in 2010?
	       	if (muon.dB(muon.PV3D) / muon.edB(muon.PV3D)) > 4:
	       		return False
	    
	       	# Maximum chi2
	       	if muon.normChi2() > 10:
	       		return False

	       	# Minimum number of valid hits on the global track. 
	       	if muon.numberOfValidHits() < 10:
	       		return False

	       	return True

	def plotter(self):
		"""
		Plots the histograms
		"""

		P.figure()
		P.hist(self.zMass, bins = 100, normed=1, alpha=0.5, label="Good Muons")
		P.hist(self.badZMass, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
		P.xlabel("Invariant mass (GeV/c2)")
		P.ylabel("frequency")
		P.legend(loc='upper right')

		P.figure()
		P.hist(self.zPt, bins = 100, normed=1, alpha=0.5, label="Good Muons")
		P.hist(self.badZPt, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
		P.xlabel("Total pt (GeV/c)")
		P.ylabel("frequency")
		P.legend(loc='upper right')

		P.figure()
		P.hist(self.zPt1, bins = 100, normed=1, alpha=0.5, label="Good Muons")
		P.hist(self.badZPt1, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
		P.xlabel("pt_1 (GeV/c)")
		P.ylabel("frequency")
		P.legend(loc='upper right')

		P.figure()
		P.hist(self.zPt2, bins = 100, normed=1, alpha=0.5, label="Good Muons")
		P.hist(self.badZPt2, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
		P.xlabel("pt_2 (GeV/c)")
		P.ylabel("frequency")
		P.legend(loc='upper right')


		P.figure()
		P.hist(self.eta, bins = 100, normed=1, alpha=0.5, label="Good Muons")
		P.hist(self.badEta, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
		P.xlabel("Eta")
		P.ylabel("frequency")
		P.legend(loc='upper right')

	#	P.figure()
	#	P.hist(self.chi2, bins = 100, normed=1, alpha=0.5, label="Good Muons")
	#	P.hist(self.badChi2, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
	#	P.xlabel("Chi**2")
	#	P.ylabel("frequency")
	#	P.legend(loc='upper right')


		P.figure()
		P.hist(self.numValidHits, bins = 100, normed=1, alpha=0.5, label="Good Muons")
		P.hist(self.badNumValidHits, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
		P.xlabel("Number of valid hits")
		P.ylabel("frequency")
		P.legend(loc='upper right')
		P.show()


	def gaussianFit(self):
		
		#bins = 100
		(mu, sigma) = norm.fit(self.zMass)

		n,bins,patches=P.hist(self.zMass, 100, normed=1)

		y = mlab.normpdf( bins, mu, sigma)
		l = P.plot(bins, y, 'r--', linewidth=2)

		P.xlabel("Invariant mass (GeV/c2)")
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


		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			selectedMuons = []


			muons = self.getMuons(event)
			vertex = self.getVertex(event)
			

			
			
			for muon in muons: # it applies the cutsConfig

			#	self.chi2.append(muon.normChi2())

				if self.selectMuons(muon, vertex) == True:
					selectedMuons.append(muon)
			#		self.badChi2.append(muon.normChi2())
				else:
					continue 

			numMuons = len(selectedMuons)
			if selectedMuons < 2: continue

			for outer in xrange(numMuons-1): #outer loop
				outerMuon=selectedMuons[outer]

				for inner in xrange(outer+1, numMuons): #inner loop
					innerMuon=selectedMuons[inner]

					if outerMuon.charge() * innerMuon.charge() >= 0:
						continue

					muPair = LeptonPair(innerMuon, outerMuon) #sum of the four-momentums of both muons
					

					if not ((muPair.mass() > self.cutsConfig.mass_min) and (muPair.mass() < 120)):
						continue

					self.zMass.append(muPair.mass())
					self.zPt.append(muPair.pt())
					self.zPt2.append(muPair.pt2())
					self.zPt1.append(muPair.pt1())
					self.eta.append(muPair.eta1())
					self.eta.append(muPair.eta2())
				#	self.chi2.append(muPair.chi2())
				#	self.chi2.append(muPair.chi1())
					self.numValidHits.append(muPair.numValidHits1())
					self.numValidHits.append(muPair.numValidHits2())

					print muPair.mass()


			# Without selecting the good muons:
			numBadMuons=len(muons)
			for outer in xrange(numBadMuons-1): #outer loop
				outerMuon=muons[outer]

				for inner in xrange(outer+1, numBadMuons): #inner loop
					innerMuon=muons[inner]

					if outerMuon.charge() * innerMuon.charge() >= 0:
						continue

					badMuPair = LeptonPair(innerMuon, outerMuon) #sum of the four-momentums of both muons
					

					if not ((badMuPair.mass() > self.cutsConfig.mass_min) and (badMuPair.mass() < 120)):
						continue

					self.badZMass.append(badMuPair.mass())
					self.badZPt.append(badMuPair.pt())
					self.badZPt2.append(badMuPair.pt2())
					self.badZPt1.append(badMuPair.pt1())
					self.badEta.append(badMuPair.eta1())
					self.badEta.append(badMuPair.eta2())
				#	self.badChi2.append(badMuPair.chi2())
				#	self.badChi2.append(badMuPair.chi1())
			#		self.badNumValidHits.append(badMuPair.numValidHits1())
			#		self.badNumValidHits.append(badMuPair.numValidHits2())


					print badMuPair.mass()
					print ""




