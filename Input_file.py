import numpy as np
import math
"""Created a Class of parameters that govern a Monte Carlo simulation of light transport in tissue"""
class Mc_variables:
  def __init__(self):
    self.absv = 0.0 
    self.dr = 0.05
    self.dz = 0.02
    self.Rmax = 1.0    
    self.Hmax = 2.0
    self.rb = 0.0    #Backscattered reflection
    self.rt = 0.0    #transmitted reflection
    self.pi = math.pi
    self.ua = 10.0    #Absorption coefficient
    self.us = 90.0    #Scattering coefficient
    self.g = 0.0    #Anisotropy
    self.n = 1.5    #referactive index of tissue medium
    self.rs=0.04    #Specular reflection
    self.heatrb = np.zeros((20))
    self.heatrt = np.zeros((20))
    self.heata = np.zeros((100,20))
    self.ca= math.sqrt(1.0-(1.0/(self.n*self.n)))    #cos thc
    self.thc=math.acos(self.ca)    #Critical angle
    self.alb=self.us/(self.us+self.ua)    #Albedo
    self.photons=100000