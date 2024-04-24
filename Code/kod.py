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
                if line.startswith('# -----------------------------------------------------'):
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
                    xyz_ref = list(map(float, input("Podaj XYZ referencyjne oddzielone przecinkiem (x,y,z): ").split(',')))
                    result = self.xyz2neu(xyz_ref, data)
                elif transform_type == '4':
                    xyz_ref = list(map(float, input("Podaj XYZ referencyjne oddzielone przecinkiem (x,y,z): ").split(',')))
                    result = self.neu2xyz(xyz_ref, data)
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
                           file.write("{:<15}{:<15}{:<15}\n".format(phi_str, lam_str, h_str))
                       elif format_choice == '1':  # degrees_decimal
                           phi_str = f"{phi:.8f}"
                           lam_str = f"{lam:.8f}"
                           h_str = f"{h:.8f}"
                           file.write("{:<15}{:<15}{:<15}\n".format(phi_str, lam_str, h_str))

                       

    
if __name__ =="__main__":
    #Wczytanie modelu 
    geo = Transformacje_wspolrzednych()
    
    model = input('wybierz model: 1 - wgs84 / 2 - grs80 / 3 - Krasowski: ')
    transform_type = input('Wybierz transformację (1 - xyz2blh, 2 - blh2xyz, 3 - xyz2neu, 4 - neu2xyz): ')
    input_file = input('Podaj nazwę pliku ze współrzędnymi: ')
    output_file = input('Podaj nazwę pliku pod jakim chcesz zapisać transformowane współrzędne: ')
    format_choice = input('Wybierz format wyników (1 - degrees_decimal, 2 - dms, 3 - XYZ): ')
    transformer = Transformacje_wspolrzednych(model)
    transformer.perform_transform(transform_type,input_file, output_file)
        
        
        
        
    
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # def xyz_to_neu(self, x, y, z, x_ref, y_ref, z_ref):
        
        
    #     pass

    # def bl_to_2000(self, fi, lam, h, elipsoida):
    #     # Metoda przekształcająca współrzędne BL na 2000
    #     pass

    # def bl_to_1992(self, fi, lam, h, elipsoida):
        
    #     pass