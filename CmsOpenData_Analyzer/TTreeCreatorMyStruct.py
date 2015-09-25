# Name: TTreeCreatorMyStruct.py
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
from DataFormats.FWLite import Events, Handle
import array


npart_max = 1000

ROOT.gROOT.ProcessLine(
"struct MyStruct {\
   Double_t     pt;\
   Double_t     eta;\
   Double_t     px;\
   Double_t     py;\
   Double_t     pz;\
   Double_t     energy;\
   Double_t     vertex_z;\
   UChar_t       isGlobal;\
   UChar_t       isTracker;\
   Double_t     dB;\
   Double_t     edB;\
   Double_t     isolation_sumPt;\
   Double_t     isolation_emEt;\
   Double_t     isolation_hadEt;\
   UChar_t        numberOfValidHits;\
   Double_t     normChi2;\
   Float_t     charge;\
};" );


ROOT.gROOT.ProcessLine(
"struct MyStructEvent {\
   Double_t     z;\
   Int_t     npart;\
};" );



class TTreeCreator(object):

	def __init__(self, data_files):

		self.muonHandle = Handle('std::vector<pat::Muon>')
		self.vertexHandle = Handle('std::vector<reco::Vertex>')	
		self.electronHandle = Handle('std::vector<pat::Electron>')

		self.events = Events(data_files)

		self.mystruct_muons = ROOT.MyStruct()
#		self.mystruct_electrons = ROOT.MyStruct()
		self.mystruct_event = ROOT.MyStructEvent()

		self.f = ROOT.TFile("mytree.root","RECREATE")
		self.tree=ROOT.TTree("test","test tree")
		
		self.Muon_pt = array.array("d", [1000])
		self.Muon_eta = array.array("d", [1000])
		self.Muon_px = array.array("d", [1000])
		self.Muon_py = array.array("d", [1000])
		self.Muon_pz = array.array("d", [1000])
		self.Muon_energy = array.array("d", [1000])
		self.Muon_vertex_z = array.array("d", [1000])
		self.Muon_isGlobalMuon = array.array("B", [1000])
		self.Muon_isTrackerMuon = array.array("B", [1000])
		self.Muon_dB = array.array("d", [1000])
		self.Muon_edB = array.array("d", [1000])
		self.Muon_isolation_sumPt = array.array("d", [1000])
		self.Muon_isolation_emEt = array.array("d", [1000])
		self.Muon_isolation_hadEt = array.array("d", [1000])
		self.Muon_numberOfValidHits = array.array("b", [1000])
		self.Muon_normChi2 = array.array("d", [1000])
		self.Muon_charge = array.array("d", [1000])

		self.Vertex_z = array.array("d", [1000])
		self.npart = array.array("I", [1000])


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

	#	self.tree.Branch("Muon", self.mystruct_muons, "pt/D:eta/D:px/D:py/D:energy/D:vertex_z/D:isGlobal/B:isTracker/B:edB/D:isolation_sumPt/D:isolation_emEt/D:isolation_hadEt/D:numberOfValidHits/I:normChi2/D:charge/F")


		self.tree.Branch("Muon_pt", ROOT.AddressOf(self.mystruct_muons, "pt"), "pt[1000]/D")
		self.tree.Branch("Muon_eta", ROOT.AddressOf(self.mystruct_muons, "eta"), "eta[1000]/D")
		self.tree.Branch("Muon_px", ROOT.AddressOf(self.mystruct_muons, "px"), "px[1000]/D")
		self.tree.Branch("Muon_py", ROOT.AddressOf(self.mystruct_muons, "py"), "py[1000]/D")
		self.tree.Branch("Muon_pz", ROOT.AddressOf(self.mystruct_muons, "pz"), "pz[1000]/D")
		self.tree.Branch("Muon_energy", ROOT.AddressOf(self.mystruct_muons, "energy"), "energy[1000]/D")
		self.tree.Branch("Muon_vertex_z", ROOT.AddressOf(self.mystruct_muons, "vertex_z"), "vertex_z[1000]/D")
		self.tree.Branch("Muon_isGlobalMuon", ROOT.AddressOf(self.mystruct_muons, "isGlobal"), "isGlobal[1000]/b")
		self.tree.Branch("Muon_isTrackerMuon", ROOT.AddressOf(self.mystruct_muons, "isTracker"), "isTracker[1000]/b")
		self.tree.Branch("Muon_dB", ROOT.AddressOf(self.mystruct_muons, "dB"), "dB[1000]/D")
		self.tree.Branch("Muon_edB", ROOT.AddressOf(self.mystruct_muons, "edB"), "edB[1000]/D")
		self.tree.Branch("Muon_isolation_sumPt", ROOT.AddressOf(self.mystruct_muons, "isolation_sumPt"), "isolation_sumPt[1000]/D")
		self.tree.Branch("Muon_isolation_emEt", ROOT.AddressOf(self.mystruct_muons, "isolation_emEt"), "isolation_emEt[1000]/D")
		self.tree.Branch("Muon_isolation_hadEt", ROOT.AddressOf(self.mystruct_muons, "isolation_hadEt"), "isolation_hadEt[1000]/D")
		self.tree.Branch("Muon_numberOfValidHits", ROOT.AddressOf(self.mystruct_muons, "numberOfValidHits"), "numberOfValidHits[1000]/b")
		self.tree.Branch("Muon_normChi2", ROOT.AddressOf(self.mystruct_muons, "normChi2"), "normChi2[1000]/D")
		self.tree.Branch("Muon_charge", ROOT.AddressOf(self.mystruct_muons, "charge"), "charge[1000]/F")

	#	self.tree.Branch("Electron", self.mystruct_electrons, "pt/D:eta/D:px/D:py/D:energy/D:vertex_z/D:isGlobal/B:isTracker/B:dB/D:edB/D:isolation_sumPt/D:isolation_emEt/D:isolation_hadEt/D:numberOfValidHits/I:normChi2/D:charge/F")
		
		self.tree.Branch("Vertex_z", ROOT.AddressOf(self.mystruct_event,"z"), "z/D")
		self.tree.Branch("event_npart", ROOT.AddressOf(self.mystruct_event, "npart"), "npart/I")



		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			

			muons = self.getMuons(event)
			#electrons = self.getElectrons(event)
			vertex = self.getVertex(event)
			
			self.mystruct_event.z = vertex.z()

		#	self.tree.Fill()
			
			for i, muon in enumerate(muons): 
				
					
				self.mystruct_muons.pt[i]=muon.pt()
				self.mystruct_muons.eta[i]=muon.eta()
				self.mystruct_muons.px[i]=muon.px()
				self.mystruct_muons.py[i]=muon.py()
				self.mystruct_muons.pz[i]=muon.pz()
				self.mystruct_muons.energy[i]=muon.energy()
				self.mystruct_muons.vertex_z[i]=muon.vertex().z()
				self.mystruct_muons.isGlobal[i]=muon.isGlobalMuon()
				self.mystruct_muons.isTracker[i]=muon.isTrackerMuon()
				self.mystruct_muons.dB[i]=muon.dB(muon.PV3D)
				self.mystruct_muons.edB[i]=muon.edB(muon.PV3D)
				self.mystruct_muons.isolation_sumPt[i]=muon.isolationR03().sumPt
				self.mystruct_muons.isolation_emEt[i]=muon.isolationR03().emEt
				self.mystruct_muons.isolation_hadEt[i]=muon.isolationR03().hadEt
				self.mystruct_muons.charge[i]=muon.charge()
				
				if not muon.globalTrack().isNull():

					self.mystruct_muons.numberOfValidHits[i]=muon.numberOfValidHits()
					self.mystruct_muons.normChi2[i]=muon.normChi2()

				else:
					self.mystruct_muons.numberOfValidHits[i]= 0
					self.mystruct_muons.normChi2[i]= 0.0

							
			self.mystruct_event.npart = self.mystruct_muons.pt.size()

			self.tree.Fill()



		self.f.Write()
		self.f.Close()

		







	def process2(self, maxEv = -1):
		"""
		maxEv: maximum number of processed events
		       maxEv=-1 runs over all the events

		It selects the good muons applying the cut configuration
		and paires up them creating objects of the class LeptonPair.
		It gets the mass of every pair and adds the one which approaches 
		the most to the Z boson's mass to the list self.zMass.  

		"""
		self.tree.Branch("Vertex_z", vertex_z, "Vertex_z/D")

		self.tree.Branch("npart", npart, "npart/I")

		self.tree.Branch("Muon_isGlobalMuon", Muon_isGlobalMuon, "Muon_isGlobalisMuon[1000]/B")


		self.tree.Branch("Muon_pt", Muon_pt, "Muon_pt[1000]/D")

		self.tree.Branch("Muon_eta", Muon_eta, "Muon_eta[1000]/D")

		self.tree.Branch("Muon_px", Muon_px, "Muon_px[1000]/D")

		self.tree.Branch("Muon_py", Muon_py, "Muon_py[1000]/D")

		self.tree.Branch("Muon_energy", Muon_energy, "Muon_energy[1000]/D")

		self.tree.Branch("Muon_isTrackerMuon", Muon_isTrackerMuon, "Muon_isTrackerMuon[1000]/B")

		self.tree.Branch("Muon_dB", Muon_dB, "Muon_dB[1000]/D")

		self.tree.Branch("Muon_edB", Muon_edB, "Muon_edB[1000]/D")

		self.tree.Branch("Muon_isolation_sumPt", Muon_isolation_sumPt, "Muon_isolation_sumPt[1000]/D")

		self.tree.Branch("Muon_isolation_emEt", Muon_isolation_emEt, "Muon_isolation_emEt[1000]/D")

		self.tree.Branch("Muon_isolation_hadEt", Muon_isolation_hadEt, "Muon_isolation_hadEt[1000]/D")

		self.tree.Branch("Muon_numberOfValidHits", Muon_numberOfValidHits, "Muon_numberOfValidHits[1000]/I")

		self.tree.Branch("Muon_normChi2", Muon_normChi2, "Muon_normChi2[1000]/D")

		self.tree.Branch("Muon_charge", Muon_charge, "Muon_charge[1000]/F")


		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			

			muons = self.getMuons(event)
			#electrons = self.getElectrons(event)
			vertex = self.getVertex(event)
			
			self.mystruct_event.z = vertex.z()

		#	self.tree.Fill()
			
			for i, muon in enumerate(muons): 
				
					
				self.mystruct_muons.pt[i]=muon.pt()
				self.mystruct_muons.eta[i]=muon.eta()
				self.mystruct_muons.px[i]=muon.px()
				self.mystruct_muons.py[i]=muon.py()
				self.mystruct_muons.pz[i]=muon.pz()
				self.mystruct_muons.energy[i]=muon.energy()
				self.mystruct_muons.vertex_z[i]=muon.vertex().z()
				self.mystruct_muons.isGlobal[i]=muon.isGlobalMuon()
				self.mystruct_muons.isTracker[i]=muon.isTrackerMuon()
				self.mystruct_muons.dB[i]=muon.dB(muon.PV3D)
				self.mystruct_muons.edB[i]=muon.edB(muon.PV3D)
				self.mystruct_muons.isolation_sumPt[i]=muon.isolationR03().sumPt
				self.mystruct_muons.isolation_emEt[i]=muon.isolationR03().emEt
				self.mystruct_muons.isolation_hadEt[i]=muon.isolationR03().hadEt
				self.mystruct_muons.charge[i]=muon.charge()
				
				if not muon.globalTrack().isNull():

					self.mystruct_muons.numberOfValidHits[i]=muon.numberOfValidHits()
					self.mystruct_muons.normChi2[i]=muon.normChi2()

				else:
					self.mystruct_muons.numberOfValidHits[i]= 0
					self.mystruct_muons.normChi2[i]= 0.0

							
			self.mystruct_event.npart = self.mystruct_muons.pt.size()

			self.tree.Fill()



		self.f.Write()
		self.f.Close()
