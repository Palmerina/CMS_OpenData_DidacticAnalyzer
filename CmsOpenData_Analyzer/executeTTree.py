# Name: execute.py
#
# CMS Open Data
#
# Description: 
#
# Returns: 

import ROOT
from TwoMuonAnalyzerTTree import TwoMuonAnalyzer
from CutsConfig import CutsConfig




cutsConfig = CutsConfig()

analyzer = TwoMuonAnalyzer(cutsConfig) # creates an object of the TwoMuonAnalyzer class


analyzer.process()

# Exercise 1
analyzer.plotter1()

# Exercise 2
analyzer.plotter2()

# Exercise 3
#analyzer.gaussianFit()


