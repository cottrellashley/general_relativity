import sympy as smp
import numpy as np
import itertools as it


class NotableTensors:
 
 
    def __init__(self, Coordinates):
        self.Coordinates = Coordinates
 
    
    def Reissner_Nordstrom(self,S_radius : object, Charge : object ):
        C = self.Coordinates
        P = S_radius
        Q = Charge
        
        return smp.MutableDenseNDimArray([[-(1 - (P)/(C[1]) + (Q**2)/(C[1]**2)),0,0,0],[0,1/(1 - (P)/(C[1])  + (Q**2)/(C[1]**2)),0,0],[0,0,C[1]**2,0], [0,0,0,C[1]**2*smp.sin(C[2])**2]])

    def Schwarzschild(self, S_radius : object ):
        
        C = self.Coordinates
        P = S_radius
        
        return smp.MutableDenseNDimArray([[-(1 - (P)/(C[1])),0,0,0],[0,1/(1 - (P)/(C[1])),0,0],[0,0,C[1]**2,0],[0,0,0,C[1]**2*smp.sin(C[2])**2]])
       