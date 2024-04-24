# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:49:52 2024

@author: weron
"""
import numpy as np
from math import *

class Transformacje_wspolrzednych:
    def __init__(self, model: str = "1"):
        """
        Okresla uklad w jakim beda liczone kolejne funkcje ["wgs84" / "grs80" / "Krasowski"]
        """
        if model == "1":
            self.a = 6378137.0 
            self.b = 6356752.3142 
        elif model == "2":
            self.a = 6378137.0
            self.b = 6356752.3141
        elif model == "3":
            self.a = 6378245 
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.e2 = (2 * self.flattening - self.flattening ** 2)
        

    
    def xyz2blh(self, xyz, format = '1'):
        X,Y,Z = xyz
        p = np.sqrt(X**2 + Y**2)
        lam = np.arctan2(Y,X)
        phi = np.arctan(Z/(p*(1-self.e2)))
        while True:
            N = self.a/(1 - self.e2 * np.sin(phi)**2)**(1/2)
            h = p/np.cos(phi) - N
            phi_poprzednie = phi
            phi = np.arctan(Z / (p * (1 - (self.e2 * (N / (N + h))))))
            if abs(phi - phi_poprzednie) < (0.000001/206265):
             break
         
        if format == '1':
             phi = np.degrees(phi)
             lam = np.degrees(lam)
      
        elif format == '2':
          phi_deg = np.degrees(phi)
          phi_min, phi_sec = divmod(phi_deg * 3600, 60)
          phi = (int(phi_min), int(phi_sec))
          
          lam_deg = np.degrees(lam)
          lam_min, lam_sec = divmod(lam_deg * 3600, 60)
          lam = (int(lam_min), int(lam_sec))
      
        result = [phi, lam, h]
        return result
    
    def blh2xyz(self, plh, format = '3'):
        phi, lam, h = plh
        N = self.a / np.sqrt(1 - self.e2 * np.sin(phi)**2)
        X = (N + h) * np.cos(phi) * np.cos(lam)
        Y = (N + h) * np.cos(phi) * np.sin(lam)
        Z = (N * (1-self.e2) + h) * np.sin(phi)
        result = [X, Y, Z]
        return result    
    
    def sigma(self, f):
        A0 = 1 - self.e2/4 - 3*self.e2**2/64 - 5*self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + 3*self.e2**3/4)
        A6 = (35*self.e2**3)/3072
        si = self.a * (A0 * f - A2 * np.sin(2*f) + A4 * np.sin(4*f) - A6 * np.sin(6*f))
        return(si)
    
    def Np(self,f):
        N = self.a / np.sqrt(1 - self.e2 * np.sin(f)**2)
        return(N)


    def Mp(self, f):
        M = (self.a * (1 - self.e2)) / np.sqrt((1 - self.e2 * np.sin(f)**2)**3)
        return(M)
    
    def bl2PL1992(self, plh, format = '4'):
        m1992 = 0.9993
        lam0 = np.deg2rad(19)
        phi,lam,h = plh
        b2 = self.a**2 * (1-self.e2)
        e22 = (self.a**2 - b2) / b2
        dl = lam - lam0
        t = tan(f)
        eta2 = e22 * (cos(f))**2
        N = self.Np(phi)
        si = self.sigma(f)
        xgk = si + (dl**2/2) * N * np.sin(phi) * np.cos(phi) * ((1 + (dl**2/12)*(np.cos(phi))**2 * (5 -t**2 +9*eta2 + 4*eta2**2) + (dl**4/360) * np.cos(phi)**4 * (61 - 58 * t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2)))
        ygk = dl * N * np.cos(phi) * (1 + (dl**2/6) * np.cos(phi)**2 * (1 - t**2 + eta2) + (dl**4/120) * cos(phi)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*eta2*t**2))
        x1992 = xgk * m1992 - 5300000
        y1992 = ygk * m1992 + 500000
        result = [x1992, y1992]
        return result
        
        
    def PL1992tobl(self,x92y92, format = '5'):
        x92, y92 = x92y92
        m1992 = 0.9993
        xgk = (x92 +5300000) / m1992
        ygk = (y92 -500000) / m1992
        
        A0 = 1 - self.e2/4 - 3*self.e2**2/64 - 5*self.e2**3/256
        f1 = xgk / (self.a * A0)
        while True:
            fs = f1
            s = self.sigma(f1)
            f1 = f1 + (xgk - s)/(self.a * A0)
            if abs(f1 - fs) < (0.000001/206265):
                break
        return(f1)
        
        b2 = self.a**2 * (1-self.e2)
        e22 = (self.a**2 - b2) / b2
        M1 = self.Mp(f1) 
        N1 = self.Np(f1)
        t1 = tan(f1)
        eta21 = e22 * (cos(f1))**2 
        phi = f1 - (ygk**2 * t1 / (2*M1*N1)) * (1-(ygk**2/(12*N1**2)) * (5 + 3 * t1**2 + eta21 - 9 * eta21 * t1**2 - 4 * eta21**2) + (ygk**4/(360 * N1**4)) * (61 + 90 * t1**2 + 45 * t1**4))
        lam = l0 + (ygk / (N1 * np.cos(f1))) * ((1 - (ygk**2 / (6 *N1**2)) * (1 + 2*t1**2 + eta21) + (ygk**4 / (120*N1**4)) * (5 + 28*t1**2 +24*t1**4 +6*eta21 + 8*eta21*t1**2)))
        result = [phi,lam]
        return result
    
    def bl2PL2000(self, plh, format = '6'):
        m2000=0.999923
        phi, lam, h = plh
        if lam >=np.deg2rad(13.5) and lam <= np.deg2rad(16.5):
            strefa = 5
            lam0 = np.deg2rad(15)
        elif lam >np.deg2rad(16.5) and lam <= np.deg2rad(19.5):
            strefa = 6
            lam0 = np.deg2rad(18)
        elif lam >np.deg2rad(19.5) and lam <= np.deg2rad(22.5):
            strefa =7
            lam0 = np.deg2rad(21)
        elif lam >np.deg2rad(22.5) and lam <= np.deg2rad(25.5):
            strefa = 8
            lam0 = np.deg2rad(24)
        b2 = self.a**2 * (1-self.e2)
        e22 = (self.a**2 - b2) / b2
        dl = lam - lam0
        t = tan(f)
        eta2 = e22 * (cos(f))**2
        N = self.Np(phi)
        si = self.sigma(phi)
        xgk = si + (dl**2/2) * N * np.sin(f) * np.cos(f) * ((1 + (dl**2/12)*(np.cos(f))**2 * (5 -t**2 +9*eta2 + 4*eta2**2) + (dl**4/360) * np.cos(f)**4 * (61 - 58 * t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2)))
        ygk = dl * N * np.cos(f) * (1 + (dl**2/6) * np.cos(f)**2 * (1 - t**2 + eta2) + (dl**4/120) * cos(f)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*eta2*t**2))
        x2000 = xgk * m2000
        y2000 = ygk * m2000 + strefa * 1000000 + 500000
        result = [x2000, y2000]
        return result
    
    
    def PL2000tobl(self,x2000y2000, format = '7'):
        m2000 = 0.999923
        x2000, y2000 = x2000y2000
        for i in y2000:
            if str(i).startswith('7'):
                nr = 7
            elif str(i).startswith('6'):
                nr = 6
            elif str(i).startswith('8'):
                nr = 8
            elif str(i).startswith('5'):
                nr = 5
            else:
                print("Nie można określić strefy")
            
        xgk = x2000 / m2000
        ygk = (y2000 - nr * 1000000 - 500000) / m2000
        
        A0 = 1 - self.e2/4 - 3*self.e2**2/64 - 5*self.e2**3/256
        f1 = xgk / (self.a * A0)
        while True:
            fs = f1
            s = self.sigma(f1)
            f1 = f1 + (xgk - s)/(self.a * A0)
            if abs(f1 - fs) < (0.000001/206265):
                break
        return(f1)
        
        b2 = self.a**2 * (1-self.e2)
        e22 = (self.a**2 - b2) / b2
        M1 = self.Mp(f1) 
        N1 = self.Np(f1)
        t1 = tan(f1)
        eta21 = e22 * (cos(f1))**2 
        phi = f1 - (ygk**2 * t1 / (2*M1*N1)) * (1-(ygk**2/(12*N1**2)) * (5 + 3 * t1**2 + eta21 - 9 * eta21 * t1**2 - 4 * eta21**2) + (ygk**4/(360 * N1**4)) * (61 + 90 * t1**2 + 45 * t1**4))
        lam = l0 + (ygk / (N1 * np.cos(f1))) * ((1 - (ygk**2 / (6 *N1**2)) * (1 + 2*t1**2 + eta21) + (ygk**4 / (120*N1**4)) * (5 + 28*t1**2 +24*t1**4 +6*eta21 + 8*eta21*t1**2)))
        result = [phi,lam]
        return result
    
    
    
    
    def perform_transform(self, transform_type, input_file, output_file, format='2'):
        """
        Wykonuje wybraną transformację na danych z pliku wejściowego i zapisuje wyniki do pliku wyjściowego.
        transform_type: Typ transformacji ('xyz2blh', 'blh2xyz', 'xyz2neu', 'neu2xyz')
        input_file: Nazwa pliku wejściowego
        output_file: Nazwa pliku wyjściowego
        """
        results = []
        with open(input_file, 'r') as file:
            found_start = False
            for line in file:
                if line.startswith('#'):
                    found_start = True
                    continue
                if not found_start:
                    continue
                data = list(map(float, line.strip().split(',')))
                if transform_type == '1':
                    result = self.xyz2blh(data)
                elif transform_type == '2':
                    result = self.blh2xyz(data)
                elif transform_type == '3':
                    result = self.bl2PL1992(data)
                elif transform_type == '4':
                    result = self.PL1992tobl(data)
                elif transform_type == '5':
                    result = self.bl2PL2000(data)   
                elif transform_type == '6':
                    result = self.PL2000tobl(data)
                else:
                    print("Niepoprawny typ transformacji.")
                    return
                
                results.append(result)
                

            with open(output_file, 'a') as file:
                # Nagłówek dla każdej funkcji
                file.write("\nWyniki transformacji:\n")

                if format_choice == '3':  # blh2xyz
                    # Nagłówek "b l h" w odpowiednich kolumnach
                    file.write("{:>6}{:>15}{:>17}\n".format("X [m]", "Y [m]", "Z [m]"))
                elif format_choice == '2' or format_choice == '1':  # xyz2blh lub inne transformacje
                    # Nagłówek "b l h" w odpowiednich kolumnach
                    file.write("{:>6}{:>15}{:>17}\n".format("b", "l", "h"))
                elif format_choice == '4' or format_choice == '6':
                    file.write("{:>6}{:>15}\n".format("X [m]", "Y [m]"))
                elif format_choice == '5' or format_choice == '7':
                    file.write("{:>6}{:>15}\n".format("b", "l"))
                

    # Wyniki
                for result in results:
                   if format_choice == '3':  # xyz2blh
                       X_str = f"{result[0]:.3f}"
                       Y_str = f"{result[1]:.3f}"
                       Z_str = f"{result[2]:.3f}"
                       file.write("{:<15}, {:<15}, {:<15}\n".format(X_str, Y_str, Z_str))
                   elif format_choice == '2' or format_choice == '1':  # xyz2blh lub inne transformacje
                       phi, lam, h = result
                       if format_choice == '2':  # dms
                           phi_deg, phi_rem = divmod(abs(phi), 1)
                           phi_min, phi_sec = divmod(phi_rem * 60, 1)
                           phi_sec *= 60
                           if phi < 0:
                               phi_str = f"-{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\""
                           else:
                               phi_str = f"{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\""
        
                           lam_deg, lam_rem = divmod(abs(lam), 1)
                           lam_min, lam_sec = divmod(lam_rem * 60, 1)
                           lam_sec *= 60
                           if lam < 0:
                               lam_str = f"-{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\""
                           else:
                               lam_str = f"{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\""
                           h_str = f"{h:.3f}"
                           file.write("{:<20}{:<20}{:<20}\n".format(phi_str, lam_str, h_str))
                       elif format_choice == '1':  # degrees_decimal
                           phi_str = f"{phi:.8f}"
                           lam_str = f"{lam:.8f}"
                           h_str = f"{h:.8f}"
                           file.write("{:<20}{:<20}{:<20}\n".format(phi_str, lam_str, h_str))
                   elif format_choice == '4' or format_choice == '6':
                       X_str = f"{result[0]:.3f}"
                       Y_str = f"{result[1]:.3f}"
                       file.write("{:<15}, {:<15}\n".format(X_str, Y_str))
                   elif format_choice == '5' or format_choice == '7':
                       phi, lam = result
                       if format_choice == '2':  # dms
                           phi_deg, phi_rem = divmod(abs(phi), 1)
                           phi_min, phi_sec = divmod(phi_rem * 60, 1)
                           phi_sec *= 60
                           if phi < 0:
                               phi_str = f"-{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\""
                           else:
                               phi_str = f"{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\""
        
                           lam_deg, lam_rem = divmod(abs(lam), 1)
                           lam_min, lam_sec = divmod(lam_rem * 60, 1)
                           lam_sec *= 60
                           if lam < 0:
                               lam_str = f"-{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\""
                           else:
                               lam_str = f"{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\""
                           file.write("{:<20}{:<20}{:<20}\n".format(phi_str, lam_str))
                       

    
if __name__ =="__main__":
    #Wczytanie modelu 
    geo = Transformacje_wspolrzednych()
    
    model = input('wybierz model: 1 - wgs84 / 2 - grs80 / 3 - Krasowski: ')
    transform_type = input('Wybierz transformację (1 - xyz2blh, 2 - blh2xyz, 3 - bl2PL1992, 4 - PL1992tobl, 5 - bl2PL2000, 6 - PL2000tobl): ')
    input_file = input('Podaj nazwę pliku ze współrzędnymi: ')
    output_file = input('Podaj nazwę pliku pod jakim chcesz zapisać transformowane współrzędne: ')
    format_choice = input('Wybierz format wyników (1 - degrees_decimal, 2 - dms, 3 - XYZ, ): ')
    transformer = Transformacje_wspolrzednych(model)
    transformer.perform_transform(transform_type,input_file, output_file)
        
        
        
        
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # def xyz_to_neu(self, x, y, z, x_ref, y_ref, z_ref):
        
        
    #     pass

    # def bl_to_2000(self, fi, lam, h, elipsoida):
    #     # Metoda przekształcająca współrzędne BL na 2000
    #     pass

    # def bl_to_1992(self, fi, lam, h, elipsoida):
        
    #     pass