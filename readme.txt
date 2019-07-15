Distribution System State Sstimation (DSSE) over 4-wire distribution models
July 3, 2019
Paulo M. De Oliveira (pm.deoliveiradejes@uniandes.edu.co)
N. Rodriguez
D. Celeita
G. Ramos
Universidad de los Andes (Colombia)
--------------------------------------------------------
PMUWLSDSSE.m run a general n-bus WLS-DSSE model over the Kersting NEV test system [1]
Excel Worksheet illustrates the method for 2-bus NEV system.
IterativePMUWLSDSSE.m performs a itereyive procedure to get key performance indexes from a sample
Both programs runs a base power flow (KerstingGeneric_powerflow.m) in order to build the measurement vector for a given noise level.
OpenDSS file: KersNeV2nNEV.dss (two bus example).

[1]Test case; Kersting NEV
Kersting, W.H. A three-phase unbalanced line model with grounded neutrals
through a resistance.  In Proceedings of the  2008 IEEE Power and Energy 
Society General Meeting-PESGM, Pittsburgh, PA, USA, 20--24 July 2008; 
pp.  12651-12652.

