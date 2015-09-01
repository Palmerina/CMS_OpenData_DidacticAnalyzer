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


class TwoMuonAnalyzer(object):
	"""
	Ni zorra
	"""

	def __init__(self, cutsConfig, data_files):

        self.muonHandle = Handle('std::vector<pat::Muon>')
        self.vertexHandle = Handle('std::vector<reco::Vertex>')
        self.cutsConfig = cutsConfig
        self.events = Events(data_files)

    def getMuons(self, event):

        event.getByLabel('patMuons', self.muonHandle)
        muons = muonHandle.product()
        return muons

	def getVertex(self, event):

        event.getByLabel('offlinePrimaryVertices', self.vertexHandle)
        vertex = vertexHandle.product()[0] #it only takes the first element which corresponds to the primary vertex
        return vertex


	def selectMuons(self, muon, vertex): #muon=getMuons(), vertex=getVertex()

		#The muon must be detected by both the tracker and the muon chambers
		if not (muon.isGlobalMuon() and muon.isTrackerMuon()):
            return False

        # Minimum transverse momentum (pt) and maximum eta angle
        if muon.pt() < cutsConfig.pt_min or abs(muon.eta()) > cutsConfig.eta_max:
            return False

        # Maximum distance of the muon respect to the vertex
        if abs(muon.vertex().z() - vertex.z()) > cutsConfig.distance:
            return False

        # Maximum impact parameter
        if muon.dB(muon.PV3D) > cutsConfig.dB_min:
            return False

        
        # I_trk + I_ECAL + I_HCAL
        # sumPt = suma de los momentos transversos
        # emEt = electromagnetic energy
        # hadEt = hadronic energy

        # Maximum energy content in that region before consider the "muon" as a jet of particles
        if (muon.isolationR03().sumPt +
                muon.isolationR03().emEt +
                muon.isolationR03().hadEt) / muon.pt() > cutsConfig.isolation:
            return False

        # muon SIP variable # Symmetrized Impact Parameter in 2010?
        if (muon.dB(muon.PV3D) / muon.edB(muon.PV3D)) > 4:
            return False

        # Maximum chi2
        if muon.normChi2() > 10:
            return False

        # Minimum number of hits
        if muon.numberOfValidHits() < 10:
            return False

        return True


	def plotter(self):


	def eventLoop(self, maxEv = 100):


		for event in self.events:

			muons = self.getMuons(event)
			vertex = self.getVertex(event)

			selectedMuons = self.selectMuons(muon, vertex) # it applies the cutsConfig

			numMuons = len(selectedMuons)
			if selectedMuons < 2: continue

			for outer in xrange(numMuons-1): #outer loop
				outerMuon=muons[outer]

				for inner in xrange(outer+1, numMuons): #inner loop
					innerMuon=muons[inner]

					if outerMuon.charge() * innerMuon.charge() >= 0:
						continue

					p4 = LeptonPair(innerMuon, outerMuon)

					#zmassHist.Fill(p4.M())