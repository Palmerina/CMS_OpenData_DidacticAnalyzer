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
from DataFormats.FWLite import Events, Handle
import array


ROOT.gROOT.ProcessLine(
"struct MyStruct {\
   Double_t     pt;\
   Double_t     eta;\
   Double_t     px;\
   Double_t     py;\
   Double_t     energy;\
   Double_t     vertex_z;\
   Bool_t       isGlobal;\
   Bool_t       isTracker;\
   Double_t     dB;\
   Double_t     edB;\
   Double_t     isolation_sumPt;\
   Double_t     isolation_emEt;\
   Double_t     isolation_hadEt;\
   Int_t        numberOfValidHits;\
   Double_t     normChi2;\
   Float_t     charge;\
};" );


ROOT.gROOT.ProcessLine(
"struct MyStructVertex {\
   Double_t     z;\
};" );



class TTreeCreator(object):

	def __init__(self, data_files):

		self.muonHandle = Handle('std::vector<pat::Muon>')
		self.vertexHandle = Handle('std::vector<reco::Vertex>')	
		self.electronHandle = Handle('std::vector<pat::Electron>')

		self.events = Events(data_files)

		self.mystruct_muons = ROOT.MyStruct()
		self.mystruct_electrons = ROOT.MyStruct()
		self.mystruct_vertex = ROOT.MyStructVertex()

		self.f = ROOT.TFile("mytree.root","RECREATE")
		self.tree=ROOT.TTree("test","test tree")


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
		electrons = self.electronHandle.product()
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


		self.tree.Branch("Muon_pt", ROOT.AddressOf(self.mystruct_muons, "pt"), "pt/D")
		self.tree.Branch("Muon_eta", ROOT.AddressOf(self.mystruct_muons, "eta"), "eta/D")
		self.tree.Branch("Muon_px", ROOT.AddressOf(self.mystruct_muons, "px"), "px/D")
		self.tree.Branch("Muon_py", ROOT.AddressOf(self.mystruct_muons, "py"), "py/D")
		self.tree.Branch("Muon_energy", ROOT.AddressOf(self.mystruct_muons, "energy"), "energy/D")
		self.tree.Branch("Muon_vertex_z", ROOT.AddressOf(self.mystruct_muons, "vertex_z"), "vertex_z/D")
		self.tree.Branch("Muon_isGlobalMuon", ROOT.AddressOf(self.mystruct_muons, "isGlobal"), "isGlobal/B")
		self.tree.Branch("Muon_isTrackerMuon", ROOT.AddressOf(self.mystruct_muons, "isTracker"), "isTracker/B")
		self.tree.Branch("Muon_dB", ROOT.AddressOf(self.mystruct_muons, "dB"), "dB/D")
		self.tree.Branch("Muon_edB", ROOT.AddressOf(self.mystruct_muons, "edB"), "edB/D")
		self.tree.Branch("Muon_isolation_sumPt", ROOT.AddressOf(self.mystruct_muons, "isolation_sumPt"), "isolation_sumPt/D")
		self.tree.Branch("Muon_isolation_emEt", ROOT.AddressOf(self.mystruct_muons, "isolation_emEt"), "isolation_emEt/D")
		self.tree.Branch("Muon_isolation_hadEt", ROOT.AddressOf(self.mystruct_muons, "isolation_hadEt"), "isolation_hadEt/D")
		self.tree.Branch("Muon_numberOfValidHits", ROOT.AddressOf(self.mystruct_muons, "numberOfValidHits"), "numberOfValidHits/I")
		self.tree.Branch("Muon_normChi2", ROOT.AddressOf(self.mystruct_muons, "normChi2"), "normChi2/D")
		self.tree.Branch("Muon_charge", ROOT.AddressOf(self.mystruct_muons, "charge"), "charge/F")

		self.tree.Branch("Electron", self.mystruct_electrons, "pt/D:eta/D:px/D:py/D:energy/D:vertex_z/D:isGlobal/B:isTracker/B:dB/D:edB/D:isolation_sumPt/D:isolation_emEt/D:isolation_hadEt/D:numberOfValidHits/I:normChi2/D:charge/F")
		
		self.tree.Branch("Vertex_z", self.mystruct_vertex, "z/D")


		for N, event in enumerate(self.events):

			if maxEv >= 0 and (N + 1) >= maxEv:
				break

			

			muons = self.getMuons(event)
			electrons = self.getElectrons(event)
			vertex = self.getVertex(event)
			
			self.mystruct_vertex.z = vertex.z()

			self.tree.Fill()
			
			for muon in muons: 
					
				self.mystruct_muons.pt=muon.pt()
				self.mystruct_muons.eta=muon.eta()
				self.mystruct_muons.px=muon.px()
				self.mystruct_muons.py=muon.py()
				self.mystruct_muons.energy=muon.energy()
				self.mystruct_muons.vertex_z=muon.vertex().z()
				self.mystruct_muons.isGlobal=muon.isGlobalMuon()
				self.mystruct_muons.isTracker=muon.isTrackerMuon()
				self.mystruct_muons.dB=muon.dB(muon.PV3D)
				self.mystruct_muons.edB=muon.edB(muon.PV3D)
				self.mystruct_muons.isolation_sumPt=muon.isolationR03().sumPt
				self.mystruct_muons.isolation_emEt=muon.isolationR03().emEt
				self.mystruct_muons.isolation_hadEt=muon.isolationR03().hadEt
				self.mystruct_muons.charge=muon.charge()
				
				if not muon.globalTrack().isNull():

					self.mystruct_muons.numberOfValidHits=muon.numberOfValidHits()
					self.mystruct_muons.normChi2=muon.normChi2()

					
				self.tree.Fill()


			for electron in electrons: 
					
				self.mystruct_electrons.pt=electron.pt()
				self.mystruct_electrons.eta=electron.eta()
				self.mystruct_electrons.px=electron.px()
				self.mystruct_electrons.py=electron.py()
				self.mystruct_electrons.energy=electron.energy()
				self.mystruct_electrons.vertex_z=electron.vertex().z()

				self.tree.Fill()



				#if not muon.globalTrack().isNull():
					#self.chi2.append(muon.normChi2())
					#self.numValidHits.append(muon.numberOfValidHits())
				#else:
					#continue 


				#	self.zMass.append(muPair.mass())
				#	self.zPt.append(muPair.pt())
				#	self.zPt2.append(muPair.pt2())
				#	self.zPt1.append(muPair.pt1())
				#	self.eta.append(muPair.eta1())
				#	self.eta.append(muPair.eta2())


		self.f.Write()
		self.f.Close()

		

		return self.tree







