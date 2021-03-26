# NSSC2

Deadline Ex1: Tuesday, April 13th 2021, 8am (submit to nssc@iue.tuwien.ac.at)

Report: https://www.overleaf.com/project/604e4e82499c472f16008a0b


Hiti, 26.03.2021: ich hab jetzt die Funktionen fürs splitting neu geschrieben und integriert

Fürs erste gibt es jetzt mal local_grid_size() und border_types(). Für die doku schaut bitte in splitting.hpp
Wichtig: die lokale Größe eines Grids ist um 1 größer für jeden ghost layer um platz für die Daten vom Nachbarn zu haben.

Aktuell ist nur der 1D-Fall vorhanden. 2D kommt in den nächsten Tagen. ich habe dafür auch ein paar Unit Test Cases in test/test_splitting.cpp geschrieben. Die Ergebnisse vom splitting sollten also passen. Die unit tests sind auch gute Beispiele wie die Funktionen zu verwenden sind


Weiters habe ich ein paar globale variablen eingeführt, damit man sie nich ständig als parameter übergeben muss:

- int g_my_rank;
- int g_n_processes;
- int g_dim;
- int g_iterations;
- size_t g_resolution;

Wenn ihr splitting.hpp in einm separaten Programm verwenden müsst, müsst ihr erst

Diese Werte sind über die CLI-Argumente festgegeben und sollten sich während der Lauzfzeit ohnehin nicht ändern.



build commands:
- make : kompiliert alle versionen inkl tests und führt sie einmal aus. Tipp: ihr könnt 'make | grep "|error|"' verwenden um schnell zu sehen ob alle versionen das gleiche ausgeben
- make test : führt die unit tests aus
- make clean : löscht alle generierten datein und deren ornder



