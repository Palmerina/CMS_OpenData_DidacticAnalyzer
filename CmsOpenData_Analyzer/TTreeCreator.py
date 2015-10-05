# Name: TTreeCreator.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 


__author__ = "Palmerina Gonzalez Izquierdo"
__copyright__ = "Copyright (C) 015 Palmerina G. I."
__license__ = "Public Domain"
__version__ = "1."
__maintainer__ = "Palmerina Gonzalez"
__email__ = "pgi2@alumnos.unican.es"

import ROOT
from DataFormats.FWLite import Events, Handle
import array




class TTreeCreator(object):

	def __init__(self, data_files):

		self.muonHandle = Handle('std::vector<pat::Muon>')
		self.vertexHandle = Handle('std::vector<reco::Vertex>')	
		self.electronHandle = Handle('std::vector<pat::Electron>')

		self.events = Events(data_files)

		self.f = ROOT.TFile("mytree.root","RECREATE")
		self.tree=ROOT.TTree("muons","muons tree")
		
		self.npart = array.array("I")
		
		self.Muon_pt = array.array("f")
		self.Muon_eta = array.array("f")
		self.Muon_px = array.array("f")
		self.Muon_py = array.array("f")
		self.Muon_pz = array.array("f")
		self.Muon_energy = array.array("f")
		self.Muon_vertex_z = array.array("f")
		self.Muon_isGlobalMuon = array.array("i")
		self.Muon_isTrackerMuon = array.array("i")
		self.Muon_dB = array.array("f")
		self.Muon_edB = array.array("f")
		self.Muon_isolation_sumPt = array.array("f")
		self.Muon_isolation_emEt = array.array("f")
		self.Muon_isolation_hadEt = array.array("f")
		self.Muon_numberOfValidHits = array.array("i")
		self.Muon_normChi2 = array.array("f")
		self.Muon_charge = array.array("f") 

		self.Vertex_z = array.array("f")


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


	def getElectrons(self, event):
		"""
		event: one element of self.events
		
		returns:
		"""

		event.getByLabel('patElectrons', self.electronHandle)
		#electrons = self.electronHandle.product()
		return electrons



	def process(self, maxEv = -1):
		"""
		maxEv: maximum number of processed events
		       maxEv=-1 runs over all the events

		It selects the good muons applying the cut configuration
		and paires up them creating objects of the class LeptonPair.
		It gets the mass of every pair and adds the one which approaches 
		the most to the Z boson's mass to the list self.zMass.  

		"""

		self.tree.Branch("Vertex_z", self.Vertex_z, "Vertex_z/F")

		self.tree.Branch("npart", self.npart, "npart/I")

		self.tree.Branch("Muon_isGlobalMuon", self.Muon_isGlobalMuon, "Muon_isGlobalisMuon/B")


		self.tree.Branch("Muon_pt", self.Muon_pt, "Muon_pt/F")

		self.tree.Branch("Muon_eta", self.Muon_eta, "Muon_eta/F")

		self.tree.Branch("Muon_px", self.Muon_px, "Muon_px/F")

		self.tree.Branch("Muon_py", self.Muon_py, "Muon_py/F")
		self.tree.Branch("Muon_pz", self.Muon_pz, "Muon_pz/F")
		self.tree.Branch("Muon_vertex_z", self.Muon_vertex_z, "Muon_vertex_z/F")

		self.tree.Branch("Muon_energy", self.Muon_energy, "Muon_energy/F")

		self.tree.Branch("Muon_isTrackerMuon", self.Muon_isTrackerMuon, "Muon_isTrackerMuon/F")

		self.tree.Branch("Muon_dB", self.Muon_dB, "Muon_dB/F")

		self.tree.Branch("Muon_edB", self.Muon_edB, "Muon_edB/F")

		self.tree.Branch("Muon_isolation_sumPt", self.Muon_isolation_sumPt, "Muon_isolation_sumPt/F")

		self.tree.Branch("Muon_isolation_emEt", self.Muon_isolation_emEt, "Muon_isolation_emEt/F")

		self.tree.Branch("Muon_isolation_hadEt", self.Muon_isolation_hadEt, "Muon_isolation_hadEt/F")

		self.tree.Branch("Muon_numberOfValidHits", self.Muon_numberOfValidHits, "Muon_numberOfValidHits/I")

		self.tree.Branch("Muon_normChi2", self.Muon_normChi2, "Muon_normChi2/F")

		self.tree.Branch("Muon_charge", self.Muon_charge, "Muon_charge/F")


		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			

			muons = self.getMuons(event)
			vertex = self.getVertex(event)
			

			self.Vertex_z.append(vertex.z())
			self.npart.append(len(muons))		

	
			for i, muon in enumerate(muons): 
				
				print i 
	

			
				self.Muon_pt.append(muon.pt())
				self.Muon_eta.append(muon.eta())
				self.Muon_px.append(muon.px())
				self.Muon_py.append(muon.py())
				self.Muon_pz.append(muon.pz())
				self.Muon_energy.append(muon.energy())
				self.Muon_vertex_z.append(muon.vertex().z())
				self.Muon_isGlobalMuon.append(muon.isGlobalMuon())
				self.Muon_isTrackerMuon.append(muon.isTrackerMuon())
				self.Muon_dB.append(muon.dB(muon.PV3D))
				self.Muon_edB.append(muon.edB(muon.PV3D))
				self.Muon_isolation_sumPt.append(muon.isolationR03().sumPt)
				self.Muon_isolation_emEt.append(muon.isolationR03().emEt)
				self.Muon_isolation_hadEt.append(muon.isolationR03().hadEt)
				self.Muon_charge.append(muon.charge())

		
		
				if not muon.globalTrack().isNull():

					self.Muon_numberOfValidHits.append(muon.numberOfValidHits())
					self.Muon_normChi2.append(muon.normChi2())

				else:
					self.Muon_numberOfValidHits.append(0)

					self.Muon_normChi2.append(0.0)

			print self.Muon_pt


		#	for j in range(len(muons)-1, 49):

			#	print j
			#	self.Muon_px.pop()
			#	self.Muon_py.pop()
			#	self.Muon_pz.pop()
			#	self.Muon_eta.pop()
			#	self.Muon_energy.pop()
			#	self.Muon_vertex_z.pop()
			#	self.Muon_isGlobalMuon.pop()
			#	self.Muon_isTrackerMuon.pop()
			#	self.Muon_dB.pop()
			#	self.Muon_edB.pop()
			#	self.Muon_isolation_sumPt.pop()
			#	self.Muon_isolation_emEt.pop()
			#	self.Muon_isolation_hadEt.pop()
			#	self.Muon_charge.pop()
			#	self.Muon_numberOfValidHits.pop()
			#	self.Muon_normChi2.pop()
			
			self.tree.Fill()

			for k in range(len(self.Muon_pt)):

				self.Muon_pt.pop()


		print "Write"

		self.f.Write()

