# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 13:13:38 2026

@author: madis
"""

def checkCapstoneConstraints(LoD,tf,pdrop,X):
  """
  

  Parameters
  ----------
  LoD : Float
    length over diameter of reactor
  tf : Float
    residence time in mins
  pdrop : float
    Pressure drop in %
  X : float
    Average conversion across reactor.

  Returns
  -------
  string explaining what conditions are 

  """
  #CONSTRAINSTS FOR CAPSTONE SOLVE
  LoverDmax= 2#length/diameter
  tfmax=15#mins
  pdropmax=8/100 #8% pressure drop
  Xmin = 76/100 #minimum Conversion
  
  
  
  results ={"LoDResult": [LoverDmax>LoD,LoD],"tfResult":[tfmax>tf,tf],"pdropResult":[pdropmax>pdrop,pdrop],"XResult":[X>Xmin,X]}
  
  falseResults = 0 #to count our false results
  for key, value in results.items():
    falseResults = 0
    if value[0] == False:
      print(f"{key} failed: value is {value[1]}")
      falseResults+=1
  if falseResults ==0:
    print("All values within simulation converge")
  return
 
checkCapstoneConstraints(3,4,15,0.15,0.25)

  