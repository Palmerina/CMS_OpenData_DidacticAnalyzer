# Name: execute.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 

import ROOT
from TwoMuonAnalyzer import TwoMuonAnalyzer
from CutsConfig import CutsConfig

# CMS data:

data_files = [
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_1.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_3.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_4.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_5.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_6.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Electron/PATtuples/Electron_PAT_data_500files_1.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Electron/PATtuples/Electron_PAT_data_500files_2.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Electron/PATtuples/Electron_PAT_data_500files_3.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Electron/PATtuples/Electron_PAT_data_500files_4.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Electron/PATtuples/Electron_PAT_data_500files_5.root',
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Electron/PATtuples/Electron_PAT_data_500files_6.root'
]

pt_min = 5 
eta_max = 2.4
distance = 0.2
dB_max = 0.02  # cm. dB=impact parameter
isolation = 0.15 #dimensionless. (sumPt+emEnergy+hadEnergy)/muon.pt = maxima energia antes de considerarlo como un jet de particulas.
mass_min=60
chi2 = 5
numValidHits = 20


maxEv = 100000 #number of processed events. maxEvents = -1 runs over all of them


cutsConfig = CutsConfig(pt_min, eta_max, distance, dB_max, isolation, mass_min)
analyzer = TwoMuonAnalyzer(cutsConfig, data_files) # creates an object of the TwoMuonAnalyzer class


analyzer.process(maxEv) 

# Plots the transverse momentum, the eta angle, the number of valid hits and the chi**2 of all muons
# and the transverse momentum and the invariant mass of all pairs of muons

analyzer.plotter1()


