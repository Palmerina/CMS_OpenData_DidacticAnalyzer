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
#from scipy import optimize
#import matplotlib.mlab as mlab
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
		self.dB = []
		self.distance = []
		self.charge=[]
		self.edB = []


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
	       	if muon.dB(muon.PV3D) > self.cutsConfig.dB_max:
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
	       	if muon.normChi2() > self.cutsConfig.chi2:
	       		return False

	       	# Minimum number of valid hits on the global track. 
	       	if muon.numberOfValidHits() < self.cutsConfig.numValidHits:
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

		fig1 = P.figure()
	
		ax_1 = fig1.add_subplot(211)
		ax_1.hist(self.zMass, bins = 80, alpha=0.5, label="Good Muons")
		ax_1.hist(self.badZMass, bins = 80, alpha=0.5, label="All Muons")
		ax_1.set_xlim(40,120)
		ax_1.set_xlabel("Invariant mass (GeV/c2)")
		ax_1.set_ylabel("frequency")
		ax_1.legend(loc='upper right')

		ax_2 = fig1.add_subplot(212)
		ax_2.hist(self.zPt, bins = 20, alpha=0.5, label="Good Muons", log=True)
		ax_2.hist(self.badZPt, bins = 500, alpha=0.5, label="All Muons", log = True)
		ax_2.set_xlim(0, 500)
		ax_2.set_xlabel("Total pt (GeV/c)")
		ax_2.set_ylabel("frequency")
		ax_2.legend(loc='upper right')

		fig2 = P.figure()

		ax_3 = fig2.add_subplot(211)
		ax_3.hist(self.zPt1, bins = 50, alpha=0.5, label="Good Muons")
		#ax_3.hist(self.badZPt1, bins = 100, normed=1, alpha=0.5, label="All Muons")
		ax_3.set_xlim(0, 50)
		ax_3.set_xlabel("pt_1 (GeV/c)")
		ax_3.set_ylabel("frequency")
		ax_3.legend(loc='upper right')

		ax_4 = fig2.add_subplot(212)
		ax_4.hist(self.zPt2, bins = 50, alpha=0.5, label="Good Muons")
		#ax_4.hist(self.badZPt2, bins = 300, normed=1, alpha=0.5, label="All Muons")
		ax_4.set_xlim(0, 200)
		ax_4.set_xlabel("pt_2 (GeV/c)")
		ax_4.set_ylabel("frequency")
		ax_4.legend(loc='upper right')


		P.figure()
		P.hist(self.zMass, bins = 80, alpha=0.5, label="Good Muons")
		P.xlabel("Invariant mass (GeV/c2)")
		P.ylabel("frequency")
		P.xlim(60, 120)


		P.figure()
		P.hist(self.eta, bins = 50, alpha=0.5, label="Good Muons", log = True)
		P.hist(self.badEta, bins = 50, alpha=0.5, label="All Muons", log = True)
		P.xlabel("Eta")
		P.ylabel("frequency")
		P.legend(loc='upper right')


		P.figure()
		P.hist(self.edB, bins = 50, log = True)
		P.xlabel("Charge")
		P.ylabel("frequency")


		P.show()


	#	P.figure()
	#	P.hist(self.chi2, bins = 100, normed=1, alpha=0.5, label="Good Muons")
	#	P.hist(self.badChi2, bins = 100, normed=1, alpha=0.5, label="Bad Muons")
	#	P.xlabel("Chi**2")
	#	P.ylabel("frequency")
	#	P.legend(loc='upper right')


	#	P.figure()
	#	P.hist(self.numValidHits, bins = 50, normed=1, alpha=0.5, label="Good Muons")
	#	P.hist(self.badNumValidHits, bins = 50, normed=1, alpha=0.5, label="Bad Muons")
	#	P.xlabel("Number of valid hits")
	#	P.ylabel("frequency")
	#	P.legend(loc='upper right')




	def gaussianFit(self):
		
		bins = 50
		mu, sigma = norm.fit(self.zMass)

		P.hist(self.zMass, bins, normed=1, alpha = 0.6)

		xmin, xmax = P.xlim()
		x = n.linspace(xmin, xmax, 100)

		y = norm.pdf(x, mu, sigma)
		P.plot(x, y, 'r--', linewidth=2)

		P.xlabel("Invariant mass (GeV/c2)")
		P.ylabel("frequency")
				
		P.show()


	def gaussianFit2(self):

		data = P.hist(self.zMass, bins = 100)

		# Equation for Gaussian
		def f(x, a, b, c):
			return a * P.exp(-(x - b)**2.0 / (2 * c**2))

		# Generate data from bins as a set of points 
		x = [0.5 * (data[1][i] + data[1][i+1]) for i in xrange(len(data[1])-1)]
		y = data[0]

		popt, pcov = optimize.curve_fit(f, x, y)

		x_fit = P.linspace(x[0], x[-1], 100)
		y_fit = f(x_fit, *popt)

		P.plot(x_fit, y_fit, lw=4, color="r")

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

		print self.events
		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			selectedMuons = []


			muons = self.getMuons(event)
			vertex = self.getVertex(event)
			

			
			
			for muon in muons: # it applies the cutsConfig
				
				self.charge.append(muon.charge())
				self.edB.append(muon.edB(muon.PV3D))
				if not muon.globalTrack().isNull():

					self.chi2.append(muon.normChi2())
					self.numValidHits.append(muon.numberOfValidHits())

				if self.selectMuons(muon, vertex) == True:
					selectedMuons.append(muon)

					if not muon.globalTrack().isNull():

						self.badChi2.append(muon.normChi2())
						self.badNumValidHits.append(muon.numberOfValidHits())
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

					muPair = LeptonPair(innerMuon, outerMuon, vertex) #sum of the four-momentums of both muons
					

					if not ((muPair.mass() > self.cutsConfig.mass_min) and (muPair.mass() < 120)):
						continue

					self.zMass.append(muPair.mass())
					self.zPt.append(muPair.pt())
					self.zPt2.append(muPair.pt2())
					self.zPt1.append(muPair.pt1())
					self.eta.append(muPair.eta1())
					self.eta.append(muPair.eta2())

			#		print muPair.mass()


			# Without selecting the good muons:
			numBadMuons=len(muons)
			for outer in xrange(numBadMuons-1): #outer loop

				if muon.pt() > 10 and abs(muon.eta()) < self.cutsConfig.eta_max:
				
					outerMuon=muons[outer]

				else:
					continue


				for inner in xrange(outer+1, numBadMuons): #inner loop

					if muon.pt() > self.cutsConfig.pt_min and abs(muon.eta()) < self.cutsConfig.eta_max:
				
						innerMuon=muons[inner]

					else:
						continue



					if outerMuon.charge() * innerMuon.charge() >= 0:
						continue

					badMuPair = LeptonPair(innerMuon, outerMuon, vertex) #sum of the four-momentums of both muons
					

					if not ((badMuPair.mass() > self.cutsConfig.mass_min) and (badMuPair.mass() < 120)):
						continue

					self.badZMass.append(badMuPair.mass())
					self.badZPt.append(badMuPair.pt())
					self.badZPt2.append(badMuPair.pt2())
					self.badZPt1.append(badMuPair.pt1())
					self.badEta.append(badMuPair.eta1())
					self.badEta.append(badMuPair.eta2())

					if abs(badMuPair.dB1())<10:
						self.dB.append(badMuPair.dB1())

					if abs(badMuPair.dB2())<10:
						self.dB.append(badMuPair.dB2())
					self.distance.append(badMuPair.distance1())
					self.distance.append(badMuPair.distance2())

			#		print badMuPair.mass()
		#			print ""


