# Name: TTreeCreatorMyStruct.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 


__author__ = "Palmerina Gonzalez Izquierdo"
__copyright__ = "Copyright (C) 5015 Palmerina G. I."
__license__ = "Public Domain"
__version__ = "1.0"
__maintainer__ = "Palmerina Gonzalez"
__email__ = "pgi25@alumnos.unican.es"

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
		
		self.Muon_pt = array.array("d", [-999.0]*50)
		self.Muon_eta = array.array("d", [-999.0]*50)
		self.Muon_px = array.array("d", [-999.0]*50)
		self.Muon_py = array.array("d", [-999.0]*50)
		self.Muon_pz = array.array("d", [-999.0]*50)
		self.Muon_energy = array.array("d", [-999.0]*50)
		self.Muon_vertex_z = array.array("d", [-999.0]*50)
		self.Muon_isGlobalMuon = array.array("i", [-999]*50)
		self.Muon_isTrackerMuon = array.array("i", [-999]*50)
		self.Muon_dB = array.array("d", [-999.0]*50)
		self.Muon_edB = array.array("d", [-999.0]*50)
		self.Muon_isolation_sumPt = array.array("d", [-999.0]*50)
		self.Muon_isolation_emEt = array.array("d", [-999.0]*50)
		self.Muon_isolation_hadEt = array.array("d", [-999.0]*50)
		self.Muon_numberOfValidHits = array.array("i", [0]*50)
		self.Muon_normChi2 = array.array("d", [-999.0]*50)
		self.Muon_charge = array.array("i", [-999]*50)

		self.Vertex_z = array.array("d", [-999.0])
		self.npart = array.array("i", [-999])


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

		self.tree.Branch("Vertex_z", self.Vertex_z, "Vertex_z[1]/D")

		self.tree.Branch("npart", self.npart, "npart[1]/I")

		self.tree.Branch("Muon_isGlobalMuon", self.Muon_isGlobalMuon, "Muon_isGlobalisMuon[50]/B")


		self.tree.Branch("Muon_pt", self.Muon_pt, "Muon_pt[50]/D")

		self.tree.Branch("Muon_eta", self.Muon_eta, "Muon_eta[50]/D")

		self.tree.Branch("Muon_px", self.Muon_px, "Muon_px[50]/D")

		self.tree.Branch("Muon_py", self.Muon_py, "Muon_py[50]/D")
		self.tree.Branch("Muon_pz", self.Muon_pz, "Muon_pz[50]/D")
		self.tree.Branch("Muon_vertex_z", self.Muon_vertex_z, "Muon_vertex_z[50]/D")

		self.tree.Branch("Muon_energy", self.Muon_energy, "Muon_energy[50]/D")

		self.tree.Branch("Muon_isTrackerMuon", self.Muon_isTrackerMuon, "Muon_isTrackerMuon[50]/B")

		self.tree.Branch("Muon_dB", self.Muon_dB, "Muon_dB[50]/D")

		self.tree.Branch("Muon_edB", self.Muon_edB, "Muon_edB[50]/D")

		self.tree.Branch("Muon_isolation_sumPt", self.Muon_isolation_sumPt, "Muon_isolation_sumPt[50]/D")

		self.tree.Branch("Muon_isolation_emEt", self.Muon_isolation_emEt, "Muon_isolation_emEt[50]/D")

		self.tree.Branch("Muon_isolation_hadEt", self.Muon_isolation_hadEt, "Muon_isolation_hadEt[50]/D")

		self.tree.Branch("Muon_numberOfValidHits", self.Muon_numberOfValidHits, "Muon_numberOfValidHits[50]/I")

		self.tree.Branch("Muon_normChi2", self.Muon_normChi2, "Muon_normChi2[50]/D")

		self.tree.Branch("Muon_charge", self.Muon_charge, "Muon_charge[50]/I")


		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			

			muons = self.getMuons(event)
			vertex = self.getVertex(event)
			

			self.Vertex_z[0] = vertex.z()
			self.npart[0] = len(muons)		

	
			for i, muon in enumerate(muons): 
				
	
				self.Muon_pt[i]=muon.pt()
				self.Muon_eta[i]=muon.eta()
				self.Muon_px[i]=muon.px()
				self.Muon_py[i]=muon.py()
				self.Muon_pz[i]=muon.pz()
				self.Muon_energy[i]=muon.energy()
				self.Muon_vertex_z[i]=muon.vertex().z()
				self.Muon_isGlobalMuon[i]=muon.isGlobalMuon()
				self.Muon_isTrackerMuon[i]=muon.isTrackerMuon()
				self.Muon_dB[i]=muon.dB(muon.PV3D)
				self.Muon_edB[i]=muon.edB(muon.PV3D)
				self.Muon_isolation_sumPt[i]=muon.isolationR03().sumPt
				self.Muon_isolation_emEt[i]=muon.isolationR03().emEt
				self.Muon_isolation_hadEt[i]=muon.isolationR03().hadEt
				self.Muon_charge[i]=muon.charge()
				
				if not muon.globalTrack().isNull():

					self.Muon_numberOfValidHits[i]=muon.numberOfValidHits()
					self.Muon_normChi2[i]=muon.normChi2()

				else:
					self.Muon_numberOfValidHits[i]= -999
					self.Muon_normChi2[i]= -999.0

							

			self.tree.Fill()


			# Empty the arrays for the next event
			for i, muon in enumerate(muons): 
				
	
				self.Muon_pt[i]=-999.0
				self.Muon_eta[i]=-999.0
				self.Muon_px[i]=-999.0
				self.Muon_py[i]=-999.0
				self.Muon_pz[i]=-999.0
				self.Muon_energy[i]=-999.0
				self.Muon_vertex_z[i]=-999.0
				self.Muon_isGlobalMuon[i]=-999
				self.Muon_isTrackerMuon[i]=-999
				self.Muon_dB[i]=-999.0
				self.Muon_edB[i]=-999.0
				self.Muon_isolation_sumPt[i]=-999.0
				self.Muon_isolation_emEt[i]=-999.0
				self.Muon_isolation_hadEt[i]=-999.0
				self.Muon_charge[i]=-999
				
				if not muon.globalTrack().isNull():

					self.Muon_numberOfValidHits[i]=-999
					self.Muon_normChi2[i]=-999.0

		print "Write"
		self.f.Write()
		self.f.Close()
