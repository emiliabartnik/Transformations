# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:49:52 2024

@author: weron
"""
import numpy as np

# a = 6378137
# e2 = 0.00669438002290

class Transformacje_wspolrzednych:
    def __init__(self, model: str = "wgs84"):
        """
        Okresla uklad w jakim beda liczone kolejne funkcje ["wgs84" / "grs80" / "Krasowski"]
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.3142 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.3141
        elif model == "Krasowski":
            self.a = 6378245 
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.e2 = (2 * self.flattening - self.flattening ** 2)
        

    
    def xyz_to_blh(self, xyz):
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
        
        return phi, lam, h
    
    
    # def read_coordinates_from_file(filename):
    #     coordinates = []
    #     try:
    #         with open(filename, "r") as file:
    #             for line in file:
    #                 values = line.strip().split()
    #                 if len(values) >= 3:
    #                     coordinates.append([float(val) for val in values[:3]])
    #                     print("Pomyślnie wczytano współrzędne z pliku.")
    #     except FileNotFoundError:
    #         print("Plik nie został znaleziony.")
    #         exit()
    #     return coordinates


    # def save_spherical_coordinates(self, coordinates, output_filename):
    #     with open(output_filename, "w") as file:
    #         for coord in coordinates:
    #             spherical = self.xyz_to_blh(*coord)
    #             file.write(str(spherical) + "\n")
    #         print("Współrzędne sferyczne zostały zapisane do pliku.")
    
    def perform_transform(self, transform_type, input_file, output_file):
        """
        Wykonuje wybraną transformację na danych z pliku wejściowego i zapisuje wyniki do pliku wyjściowego.
        transform_type: Typ transformacji ('xyz2blh', 'blh2xyz', 'xyz2neu', 'neu2xyz')
        input_file: Nazwa pliku wejściowego
        output_file: Nazwa pliku wyjściowego
        """
        results = []
        with open(input_file, 'r') as file:

        
            for line in file:
                data = list(map(float, line.strip().split(',')))
                if transform_type == 'xyz_to_blh':
                    result = self.xyz_to_blh(data)
                elif transform_type == 'blh2xyz':
                    result = self.blh2xyz(data)
                elif transform_type == 'xyz2neu':
                    xyz_ref = list(map(float, input("Podaj XYZ referencyjne oddzielone przecinkiem (x,y,z): ").split(',')))
                    result = self.xyz2neu(xyz_ref, data)
                elif transform_type == 'neu2xyz':
                    xyz_ref = list(map(float, input("Podaj XYZ referencyjne oddzielone przecinkiem (x,y,z): ").split(',')))
                    result = self.neu2xyz(xyz_ref, data)
                else:
                    print("Niepoprawny typ transformacji.")
                    return
                
                results.append(','.join(map(str, result)))

        with open(output_file, 'w') as file:
            for result in results:
                file.write(result + '\n')
    
if __name__ =="__main__":
    #Wczytanie modelu 
    geo = Transformacje_wspolrzednych()
    
    model = input('wybierz model')
    transform_type = input('Wybierz transformację (xyz_to_blh, blh2xyz, xyz2neu, neu2xyz)')
    input_file = input('Podaj nazwę pliku ze współrzędnymi')
    output_file = input('Podaj nazwę pliku pod jakim chcesz zapisać transformowane współrzędne')
    
    transformer = Transformacje_wspolrzednych(model)
    transformer.perform_transform(transform_type,input_file, output_file)
        
        
        
        
    
    # def blh_to_xyz(self, fi, lam, h, elipsoida):
    #     N = a / np.sqrt(1 - e2 * np.sin(fi)**2)
    #     X = (N + h) * np.cos(fi) * np.cos(lam)
    #     Y = (N + h) * np.cos(fi) * np.sin(lam)
    #     Z = (N * (1-e2) + h) * np.sin(fi)
    #     return(X, Y, Z)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # def xyz_to_neu(self, x, y, z, x_ref, y_ref, z_ref):
        
        
    #     pass

    # def bl_to_2000(self, fi, lam, h, elipsoida):
    #     # Metoda przekształcająca współrzędne BL na 2000
    #     pass

    # def bl_to_1992(self, fi, lam, h, elipsoida):
        
    #     pass