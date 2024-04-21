# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:49:52 2024

@author: weron
"""
import numpy as np

a = 6378137
e2 = 0.00669438002290

class Transformacje_wspolrzednych:
    def __init__(self):
        pass
    
    def xyz_to_blh(self, X,Y,Z,elipsoida):
        p = np.sqrt(X**2 + Y**2)
        lam = np.arctan2(Y,X)
        phi = np.arctan(Z/(p*(1-e2)))
        while True:
            N = a/(1 - e2 * np.sin(phi)**2)**(1/2)
            h = p/np.cos(phi) - N
            phi_poprzednie = phi
            phi = np.arctan(Z / (p * (1 - (e2 * (N / (N + h))))))
            if abs(phi - phi_poprzednie) < (0.000001/206265):
             break
        
        return phi, lam, h
    
    
    def blh_to_xyz(self, fi, lam, h, elipsoida):
        pass
    
    
    def xyz_to_neu(self, x, y, z, x_ref, y_ref, z_ref):
        pass

    def bl_to_2000(self, fi, lam, h, elipsoida):
        # Metoda przekształcająca współrzędne BL na 2000
        pass

    def bl_to_1992(self, fi, lam, h, elipsoida):
        
        pass