# Transformacje 
## Cel programu
Program służy do transponowania współrzędnych pomiędzy różnymi układami. Dostępne opcje to: 

	-{xyz2blh}: Przekształca współrzędne kartezjańskie (XYZ) na współrzędne geodezyjne (BLH),
 
	-{blh2xyz}: Przekształca współrzędne geodezyjne (BLH) na współrzędne kartezjańskie (XYZ),
 
	-{bl2PL1992}: Przekształca współrzędne geodezyjne (BL) na płaskie współrzędne (PL) w układzie 1992,
 
	-{PL1992tobl}:  Przekształca płaskie współrzędne (PL) w układzie 1992 na współrzędne geodezyjne (BL), wykona się dla wybranego modelu WGS84 oraz GRS80

	-{bl2PL2000}: Przekształca współrzędne geodezyjne (BL) na płaskie współrzędne (PL) w układzie 2000,
 
	-{PL2000tobl}: Przekształca płaskie współrzędne (PL) w układzie 2000 na współrzędne geodezyjne (BL), wykona się dla wybranego modelu WGS84 oraz GRS80
 
	-{xyz2neu}: Przekształca współrzędne kartezjańskie (XYZ) na lokalny układ współrzędnych (NEU),
 
	-{neu2XYZ}: Przekształca lokalny układ współrzędnych (NEU) na współrzędne kartezjańskie (XYZ)
 
## Instrukcja działania
Aby uruchomić program należy, na komputerze z zainstalowanym Pythonem w wersji 3.11, otworzyć wiersz poleceń. Nastepnie konieczne jest przjście do lokalizacji folderu, w którym znajduje się pobrany plik w rozszerzeniu .py. Po otworzeniu wiersza poleceń w wybranej przez nas ścieżce należy wywołać plik 'kod.py'. Można również bezpośrednio otworzyć plik bez konieczności wcześniejszego otwierania wiersza poleceń. Program uruchomi się i na początku zapyta użytkownika o wybór elipsoidy odniesienia. Dostępne są 3 opcje

![image](https://github.com/emiliabartnik/Transformations/assets/150865197/299f207f-311a-4247-9654-b947a4567d21) 

Po wybraniu odpowiedniego modelu, użytkownik będzie poproszony o wybranie rodzaju transformacji. 

![image](https://github.com/emiliabartnik/Transformations/assets/150865197/b4544802-346b-4d4c-acdc-e0896c7228c0)

Po wybraniu interesującej nas transformacji program poprosi nas o wybranie pliku tekstowego ze współrzędnymi. Plik musi znajdować się w tej samej lokalizacji co plik z naszym kodem. Plik musi mieć nagłówek, którego końcowa linijka zaczyna się od znaku "#". W pliku współrzędne muszą być umieszczone w kolejnych kolumnach. Kolejne współrzędne muszą być oddzielone od siebie przecinkiem. Wartości dziesiętne muszą być zapisywane z kropką. Należy pamiętać o dopasowaniu pliku do transformacji, z której chcemy skorzystać. Jeśli wybieramy transformacje, w której współrzędnymi wejściowymi są wartości kątowe należy podać je w stopniach, w formacie dziesiętnym. Przykładowy wygląd pliku wejściowego:

![image](https://github.com/emiliabartnik/Transformations/assets/150865197/f4573582-27b1-4113-9205-a2cd259bb0ac)

Kiedy wybraliśmy plik wejściowy, program poprosi nas o podanie nazwy pliku, do którego chcemy zapisać przeliczone współrzędne. 

Następnie przy wybraniu niektórych transformacji, program zapyta użytkownika o format w jakim chcemy zapisać współrzędne. W przypadku innych transformacji program dopasuje format wyjściowych wspołrzędnych automatycznie. 

Po wykonaniu wszystkich działań program się zamknie, a plik z naszymi wspołrzędnymi zostanie zapisany w tej samej lokalizacji, co poprzednie pliki. 

## Autorki

Za stworzenie projektu odpowiedzialne są [@emiliabartnik](https://github.com/emiliabartnik) oraz [@weronikaga](https://github.com/weronikaga).
