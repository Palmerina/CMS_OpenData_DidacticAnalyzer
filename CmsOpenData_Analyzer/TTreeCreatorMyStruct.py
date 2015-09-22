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
   Double_t     z;\
   Double_t     vertex_z;\
};" );


class TTreeCreator(object):

	def __init__(self, data_files):

		self.muonHandle = Handle('std::vector<pat::Muon>')
		self.vertexHandle = Handle('std::vector<reco::Vertex>')	
		self.electronHandle = Handle('std::vector<pat::Electron>')
		self.events = Events(data_files)
		self.zMass = []
		self.badZMass = []
		self.Pt = array.array("d",[0.])
		self.badZPt = []
		self.zPt1 = []
		self.badZPt1 = []
		self.zPt2 = []
		self.badZPt2 = []
		self.eta = array.array("d",[0.])
		self.badEta = []
		self.chi2 = []
		self.badChi2 = []
		self.numValidHits = []
		self.badNumValidHits = []
		self.dB = []
		self.distance = []


		self.mystruct_muons = ROOT.MyStruct()
		self.mystruct_electrons = ROOT.MyStruct()
		self.mystruct_vertex = ROOT.MyStruct()
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

		self.tree.Branch("Muon", self.mystruct_muons, "pt/D:eta/D:px/D:py/D:energy/D:vertex_z/D")
		self.tree.Branch("Electron", self.mystruct_electrons, "pt/D:eta/D:px/D:py/D:energy/D:vertex_z/D")
		
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







