This folder includes data from the soil moisture stations.

Two types of probes (sensors) at BFO continuously measure soil moisture:
 - SMT100
 - IMKO
 
Besides, we have two data acquisition systems at BFO:
 - Solidus  (able to handle two probes)
 - TTGOPICO (only can handle one probe)

Currently, at BFO, we have six stations as follows:
 - SMT100DataSolidus1: Two SMT100 probes with Solidus data acquisition system
 - SMT100DataSolidus2: Two SMT100 probes with Solidus data acquisition system
 - SMT100DataSolidus3: Two SMT100 probes with Solidus data acquisition system
 - ImkoDataSolidus5  : Two IMKO probes with Solidus data acquisition system
 - TTGOPICO3         : One IMKO probe with TTGOPICO data acquisition system
 - TTGOPICO4         : One IMKO probe with TTGOPICO data acquisition system

For each station mentioned above, a file in CSV format is provided in this folder:
--------------------------------------------------------------------------
For the stations:
 - SMT100DataSolidus1
 - SMT100DataSolidus2
 - SMT100DataSolidus3
 - ImkoDataSolidus5

Each station contains two probes, probe 0 and probe 1.
Each soil moisture probe measures:
 - soil moisture SM in % 
 - soil temperature T in °C. 

The files for these stations include five columns:
Column#1: Date in the UTC (Coordinated Universal Time) format
Column#2: Soil Moisture (SM) values in % for the first probe (probe 0)
Column#3: Soil Moisture (SM) values in % for the second probe (probe 1)
Column#4: Temperature (T) values in Celcius for the first probe (probe 0)
Column#5: Temperature (T) values in Celcius for the second probe (probe 1)
--------------------------------------------------------------------------
For the stations:
 - TTGOPICO3
 - TTGOPICO4

Each station contains only one probe, and the probe measures:
 - soil moisture SM in % 
 - soil temperature T in °C. 

The files for these stations include three columns:
Column#1: Date in the UTC (Coordinated Universal Time) format
Column#2: Soil Moisture (SM) values in % 
Column#3: Temperature (T) values in Celcius


The locations of the probes are: (see also the map at http://smbfo.gis.uni-stuttgart.de/)
Solidus 1 and 2: 48.3298 N 8.3274 W
Solidus 3 and TTGPICO 4  and TTGPICO 3: 48.3277 N; 8.3257 W
Solidus 5: 48.3294 N; 8.3298 W