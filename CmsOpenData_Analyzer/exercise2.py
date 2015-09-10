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

# cutsConfig parameters:

# These are the cuts applied to the muons in order to select the good ones (see TwoMuonAnalyzer.py)
# Change the values to get the Z mass boson peak

pt_min = 0 
eta_max = 100
distance = 1
dB_max = 1  # cm. dB=impact parameter
isolation = 1 #dimensionless. (sumPt+emEnergy+hadEnergy)/muon.pt = maxima energia antes de considerarlo como un jet de particulas.
mass_min=60
chi2 = 50
numValidHits = 0

maxEv = 100000 #number of processed events. maxEvents = -1 runs over all of them


cutsConfig = CutsConfig(pt_min, eta_max, distance, dB_max, isolation, mass_min)
analyzer = TwoMuonAnalyzer(cutsConfig, data_files) # creates an object of the TwoMuonAnalyzer class


analyzer.process(maxEv) 
analyzer.plotter()


