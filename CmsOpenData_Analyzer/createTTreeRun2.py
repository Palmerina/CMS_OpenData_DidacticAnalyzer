# Name: execute.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 

import ROOT
#from TTreeCreator import TTreeCreator
from TTreeCreator2 import TTreeCreator

# CMS data:

data_files = [
        'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/PATtuples/Mu_PAT_data_500files_1.root',
]

maxEv = 100000 #number of processed events. maxEvents = -1 runs over all of them

treeCreator = TTreeCreator(data_files)
tree=treeCreator.process(maxEv) 




