# Name: execute.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 

import ROOT
from TwoMuonAnalyzerTTree2 import TwoMuonAnalyzer
from CutsConfig import CutsConfig




cutsConfig = CutsConfig()

analyzer = TwoMuonAnalyzer(cutsConfig) # creates an object of the TwoMuonAnalyzer class


analyzer.process()


# Exercise 2
analyzer.plotter()
analyzer.plotHistos()

# Exercise 3
#analyzer.gaussianFit()


