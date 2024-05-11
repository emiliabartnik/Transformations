# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:49:52 2024

@author: weron
"""
import numpy as np


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
        

    
    def xyz2blh(self, xyz, forma = '1'):
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
      
        result = [np.rad2deg(phi), np.rad2deg(lam), h]
        return result
    
    def blh2xyz(self, plh, forma = '3'):
        phi1, lam1, h = plh
        phi = np.deg2rad(phi1)
        lam = np.deg2rad(lam1)
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
    
    def bl2PL1992(self, plh, forma = '4'):
        m1992 = 0.9993
        lam0 = np.deg2rad(19)
        phi1,lam1,h = plh
        phi = np.deg2rad(phi1)
        lam = np.deg2rad(lam1)
        b2 = self.a**2 * (1-self.e2)
        e22 = (self.a**2 - b2) / b2
        dl = lam - lam0
        t = np.tan(phi)
        eta2 = e22 * (np.cos(phi))**2
        N = self.Np(phi)
        si = self.sigma(phi)
        xgk = si + (dl**2/2) * N * np.sin(phi) * np.cos(phi) * ((1 + (dl**2/12)*(np.cos(phi))**2 * (5 -t**2 +9*eta2 + 4*eta2**2) + (dl**4/360) * np.cos(phi)**4 * (61 - 58 * t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2)))
        ygk = dl * N * np.cos(phi) * (1 + (dl**2/6) * np.cos(phi)**2 * (1 - t**2 + eta2) + (dl**4/120) * np.cos(phi)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*eta2*t**2))
        x1992 = xgk * m1992 - 5300000
        y1992 = ygk * m1992 + 500000
        result = [x1992, y1992]
        return result
        
        
    def PL1992tobl(self, x92y92, forma='5'):
        x92, y92 = x92y92
        m1992 = 0.9993
        l0 = np.deg2rad(19)
        xgk = (x92 + 5300000) / m1992
        ygk = (y92 - 500000) / m1992
    
        A0 = 1 - self.e2/4 - 3*self.e2**2/64 - 5*self.e2**3/256
        phi1 = xgk / (self.a * A0)
        while True:
            phis = phi1
            s = self.sigma(phi1)
            phi1 = phi1 + (xgk - s)/(self.a * A0)
            if abs(phi1 - phis) < (0.000001/206265):
                break
    
        b2 = self.a**2 * (1-self.e2)
        e22 = (self.a**2 - b2) / b2
        M1 = self.Mp(phi1) 
        N1 = self.Np(phi1)
        t1 = np.tan(phi1)
        eta21 = e22 * (np.cos(phi1))**2 
        phi = phi1 - (ygk**2 * t1 / (2*M1*N1)) * (1-(ygk**2/(12*N1**2)) * (5 + 3 * t1**2 + eta21 - 9 * eta21 * t1**2 - 4 * eta21**2) + (ygk**4/(360 * N1**4)) * (61 + 90 * t1**2 + 45 * t1**4))
        lam = l0 + (ygk / (N1 * np.cos(phi1))) * ((1 - (ygk**2 / (6 *N1**2)) * (1 + 2*t1**2 + eta21) + (ygk**4 / (120*N1**4)) * (5 + 28*t1**2 +24*t1**4 +6*eta21 + 8*eta21*t1**2)))
        result = [np.rad2deg(phi), np.rad2deg(lam)]
        return result
    


    def bl2PL2000(self, plh, forma ='6'):
        m2000 = 0.999923
        phi1, lam1, h = plh
        phi = np.deg2rad(phi1)
        lam = np.deg2rad(lam1)
        try:
            if lam >= np.deg2rad(13.5) and lam <= np.deg2rad(16.5):
                strefa = 5
                lam0 = np.deg2rad(15)
            elif lam > np.deg2rad(16.5) and lam <= np.deg2rad(19.5):
                strefa = 6
                lam0 = np.deg2rad(18)
            elif lam > np.deg2rad(19.5) and lam <= np.deg2rad(22.5):
                strefa = 7
                lam0 = np.deg2rad(21)
            elif lam > np.deg2rad(22.5) and lam <= np.deg2rad(25.5):
                strefa = 8
                lam0 = np.deg2rad(24)
            else:
                raise ValueError("Niepoprawna wartość lam")
        
            b2 = self.a**2 * (1 - self.e2)
            e22 = (self.a**2 - b2) / b2
            dl = lam - lam0
            t = np.tan(phi)
            eta2 = e22 * (np.cos(phi))**2
            N = self.Np(phi)
            si = self.sigma(phi)
            xgk = si + (dl**2/2) * N * np.sin(phi) * np.cos(phi) * ((1 + (dl**2/12) * (np.cos(phi))**2 * (5 - t**2 + 9*eta2 + 4*eta2**2) + (dl**4/360) * np.cos(phi)**4 * (61 - 58 * t**2 + t**4 + 270*eta2 - 330 * eta2 * t**2)))
            ygk = dl * N * np.cos(phi) * (1 + (dl**2/6) * np.cos(phi)**2 * (1 - t**2 + eta2) + (dl**4/120) * np.cos(phi)**4 * (5 - 18*t**2 + t**4 + 14*eta2 - 58*eta2*t**2))
            x2000 = xgk * m2000
            y2000 = ygk * m2000 + strefa * 1000000 + 500000
            result = [x2000, y2000]
            return result
        except ValueError as e:
            print(e)
            return None
    
    
    def PL2000tobl(self, x2000y2000, forma='5'):
        m2000 = 0.999923
        x2000, y2000 = x2000y2000
        l0 = np.deg2rad(19)

        # Znalezienie strefy na podstawie y2000
        if str(y2000).startswith('7'):
            nr = 7
        elif str(y2000).startswith('6'):
            nr = 6
        elif str(y2000).startswith('8'):
            nr = 8
        elif str(y2000).startswith('5'):
            nr = 5
        else:
            print("Nie można określić strefy")
            return None  # Jeśli nie można określić strefy, zwracamy None lub podejmujemy odpowiednie działania w zależności od potrzeb

        try:
            xgk = x2000 / m2000
            ygk = (y2000 - nr * 1000000 - 500000) / m2000
            
            A0 = 1 - self.e2/4 - 3*self.e2**2/64 - 5*self.e2**3/256
            phi1 = xgk / (self.a * A0)

            while True:
                phis = phi1
                s = self.sigma(phi1)
                phi1 = phi1 + (xgk - s)/(self.a * A0)
                if abs(phi1 - phis) < (0.000001/206265):
                    break

            b2 = self.a**2 * (1 - self.e2)
            e22 = (self.a**2 - b2) / b2
            M1 = self.Mp(phi1) 
            N1 = self.Np(phi1)
            t1 = np.tan(phi1)
            eta21 = e22 * (np.cos(phi1))**2 

            phi = phi1 - (ygk**2 * t1 / (2*M1*N1)) * (1-(ygk**2/(12*N1**2)) * (5 + 3 * t1**2 + eta21 - 9 * eta21 * t1**2 - 4 * eta21**2) + (ygk**4/(360 * N1**4)) * (61 + 90 * t1**2 + 45 * t1**4))
            lam = l0 + (ygk / (N1 * np.cos(phi1))) * ((1 - (ygk**2 / (6 *N1**2)) * (1 + 2*t1**2 + eta21) + (ygk**4 / (120*N1**4)) * (5 + 28*t1**2 +24*t1**4 +6*eta21 + 8*eta21*t1**2)))

            result = [np.rad2deg(phi), np.rad2deg(lam)]
            return result
        except Exception as e:
            print(f"Błąd w funkcji PL2000tobl: {e}")
            return None

    
    
    def xyz2neu (self, xyz, x0,y0,z0, forma = '8'):
        x, y, z = xyz
        x0y0z0 = [x0, y0, z0]
        phi, lam, h = self.xyz2blh(x0y0z0)
        R = np.array([[-np.sin(phi)*np.cos(lam), -np.sin(lam), np.cos(phi)*np.cos(lam)],
                      [-np.sin(phi)*np.sin(lam), np.cos(lam), np.cos(phi)*np.sin(lam)],
                      [np.cos(phi), 0, np.sin(phi)]])
        xyz_t = np.array([[x-x0],
                         [y-y0],
                         [z-z0]])
        [[E],[N],[U]] = R.T @ xyz_t
        result = [N, E, U]
        return result
    
    
    def neu2XYZ(self, NEU, x0, y0, z0, forma = '3'):
        n, e, u = NEU
        x0y0z0 = [x0, y0, z0]
        phi, lam, h = self.xyz2blh(x0y0z0)
        R = np.array([[-np.sin(phi)*np.cos(lam), -np.sin(lam), np.cos(phi)*np.cos(lam)],
                      [-np.sin(phi)*np.sin(lam), np.cos(lam), np.cos(phi)*np.sin(lam)],
                      [np.cos(phi), 0, np.sin(phi)]])
        dx = np.array([n, e, u])
        X,Y, Z = R @ dx
        result = [X, Y, Z]
        return result
    
    
    
    def perform_transform(self, transform_type, input_file, output_file, format_choice, x0 = 0, y0 = 0, z0 =0):
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
                elif transform_type == '7':
                    result = self.xyz2neu(data, x0, y0, z0)
                elif transform_type == '8':
                    result = self.neu2XYZ(data, x0, y0, z0)
                else:
                    print("Niepoprawny typ transformacji.")
                    return
                
                results.append(result)
                

            with open(output_file, 'a') as file:
                file.write("\nWyniki transformacji:\n")
            
                if format_choice == '3':  # blh2xyz
                    file.write("{:>9}{:>15}{:>17}\n".format("X [m]", "Y [m]".rjust(17), "Z [m]".rjust(15)))
                    header_written = False
                    for result in results:
                        if not header_written:  # Sprawdzamy, czy nagłówek już został zapisany
                            file.write("#-----------------------------\n")  # Zapisujemy nagłówek tylko raz
                            header_written = True
                        X_str = f"{result[0]:.3f}"
                        Y_str = f"{result[1]:.3f}".rjust(15)
                        Z_str = f"{result[2]:.3f}".rjust(15)
                        file.write("{:<10}, {:<15}, {:<15}\n".format(X_str, Y_str, Z_str))
                        
                elif format_choice == '2' or format_choice == '1':  # xyz2blh lub inne transformacje
                    file.write("{:>9}{:>15}{:>17}\n".format("b", "l".rjust(20), "h".rjust(15)))
                    header_written = False 
                    for result in results:
                        if not header_written:  # Sprawdzamy, czy nagłówek już został zapisany
                            file.write("#------------------------------------------------\n")  # Zapisujemy nagłówek tylko raz
                            header_written = True
                        phi, lam, h = result
                        if format_choice == '2':  # dms
                            phi_deg, phi_rem = divmod(abs(phi), 1)
                            phi_min, phi_sec = divmod(phi_rem * 60, 1)
                            phi_sec *= 60
                            if phi < 0:
                                phi_str = f"-{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\"".rjust(10)
                            else:
                                phi_str = f"{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\"".rjust(10)
            
                            lam_deg, lam_rem = divmod(abs(lam), 1)
                            lam_min, lam_sec = divmod(lam_rem * 60, 1)
                            lam_sec *= 60
                            if lam < 0:
                                lam_str = f"-{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\"".rjust(20)
                            else:
                                lam_str = f"{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\"".rjust(20)
                            h_str = f"{h:.3f}".rjust(12)
                            file.write("{:<10},{:<10},{:<10}\n".format(phi_str, lam_str, h_str))
                        elif format_choice == '1':  # degrees_decimal
                            phi_str = f"{phi:.8f}"
                            lam_str = f"{lam:.8f}".rjust(20)
                            h_str = f"{h:.8f}".rjust(12)
                            file.write("{:<10},{:<10},{:<10}\n".format(phi_str, lam_str, h_str))
                            
        
                
                elif format_choice == '4' or format_choice == '6':
                    file.write("{:>6}{:>15}\n".format("X [m]", "Y [m]".rjust(15)))
                    header_written = False
                    for result in results:
                        if not header_written:  # Sprawdzamy, czy nagłówek już został zapisany
                            file.write("#-----------------------------\n")  # Zapisujemy nagłówek tylko raz
                            header_written = True
                        if result is not None:
                            X_str = f"{result[0]:.3f}"
                            Y_str = f"{result[1]:.3f}".rjust(15)
                            file.write("{:<10}, {:<15}\n".format(X_str, Y_str))
                        else:
                            file.write("Error: Result is None\n")
                        
                
                elif format_choice == '5' or format_choice == '7':  # xyz2blh lub inne transformacje
                    file.write("{:>12}{:>15}\n".format("b", "l".rjust(20)))
                    header_written = False  # Zmienna do śledzenia, czy nagłówek został już zapisany
                    for result in results:
                        if not header_written:  # Sprawdzamy, czy nagłówek już został zapisany
                            file.write("#-----------------------------\n")  # Zapisujemy nagłówek tylko raz
                            header_written = True
                        if result is not None:
                            phi, lam = result
                            if format_choice == '5':  # dms
                                phi_deg, phi_rem = divmod(abs(phi), 1)
                                phi_min, phi_sec = divmod(phi_rem * 60, 1)
                                phi_sec *= 60
                                if phi < 0:
                                    phi_str = f"-{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\"".rjust(10)
                                else:
                                    phi_str = f"{int(phi_deg):02d}°{int(phi_min):02d}'{phi_sec:.5f}\"".rjust(10)

                                lam_deg, lam_rem = divmod(abs(lam), 1)
                                lam_min, lam_sec = divmod(lam_rem * 60, 1)
                                lam_sec *= 60
                                if lam < 0:
                                    lam_str = f"-{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\"".rjust(20)
                                else:
                                    lam_str = f"{int(lam_deg):02d}°{int(lam_min):02d}'{lam_sec:.5f}\"".rjust(20)
                                file.write("{:<15},{:<15}\n".format(phi_str, lam_str))
                            elif format_choice == '7':  # degrees_decimal
                                phi_str = f"{phi:.8f}"
                                lam_str = f"{lam:.8f}".rjust(15)
                                file.write("{:<15},{:<15}\n".format(phi_str, lam_str))
                    if not header_written:  # Sprawdzamy, czy nagłówek został już zapisany
                        file.write("#-----------------------------\n")  # Jeśli nie został zapisany, zapisujemy go teraz

                

                       

    
if __name__ =="__main__":
    #Wczytanie modelu 
    geo = Transformacje_wspolrzednych()
    
    model = input('wybierz model: 1 - wgs84 / 2 - grs80 / 3 - Krasowski: ')
    transform_type = input('Wybierz transformację (1 - xyz2blh, 2 - blh2xyz, 3 - bl2PL1992, 4 - PL1992tobl, 5 - bl2PL2000, 6 - PL2000tobl, 7 - xyz2neu, 8 - neu2XYZ): ')
    input_file = input('Podaj nazwę pliku ze współrzędnymi: ')
    output_file = input('Podaj nazwę pliku pod jakim chcesz zapisać transformowane współrzędne: ')
    #format_choice = input('Wybierz format wyników (1 - BLH/degrees_decimal, 2 - BLH/dms, 3 - XYZ, 4 - X92/Y92, 5 - BL/dms, 6 - X2000/Y2000, 7 - BL/degrees_decimal, 8 - NEU): ')
    
    transformer = Transformacje_wspolrzednych(model)
    if transform_type in ['7', '8']:
        x0 = float(input('Podaj wartość x0: '))
        y0 = float(input('Podaj wartość y0: '))
        z0 = float(input('Podaj wartość z0: '))
    # Wywołanie metody perform_transform z podanymi argumentami
        if transform_type == '7':
            forma = '8'
            transformer.perform_transform(transform_type, input_file, output_file, forma, x0, y0, z0)
        elif transform_type == '8':
            forma = '3'
            transformer.perform_transform(transform_type, input_file, output_file, forma, x0, y0, z0)
    if transform_type == '1':
        format_choice = input("Wybierz format wyników (1 - BLH/degrees_decimal, 2 - BLH/dms): ")
        transformer.perform_transform(transform_type, input_file, output_file, format_choice) 
    elif transform_type == '2':
        forma = '3'
        transformer.perform_transform(transform_type, input_file, output_file, forma) 
    elif transform_type == '3':
        forma = '4'
        transformer.perform_transform(transform_type, input_file, output_file, forma) 
    elif transform_type in ['4', '6']:
         format_choice = input("Wybierz format wyników (5 - BL/dms, 7 - BL/degrees_decimal): ")
         transformer.perform_transform(transform_type, input_file, output_file, format_choice) 
    elif transform_type == '5':
        forma = '6'
        transformer.perform_transform(transform_type, input_file, output_file, forma) 
    else:
        # Wywołanie metody perform_transform bez dodatkowych argumentów
        transformer.perform_transform(transform_type, input_file, output_file, format_choice)   
        
    
   
    
 