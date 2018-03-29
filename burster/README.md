Written by Maciej Dobrzynski
Centre for Mathematics and Computer Science (CWI), Amsterdam, The Netherlands.

Gillespie simulation to analyze interarrival times of product in the simple burst-generating model (the interrupted Poisson process).

Type 'make' to compile, and './slide.x' to execute the program. 
Help message is displayed if no command-line parameters are given.
Type 'make clean' to remove compiled files.

  Reactions:
  Switch:     OFF <-> ON
  Initiation: ON -> ON + P
  Product degradation: P -> null  

Description
-----------
Product P is generated with rate constant k_ini if the switch is in the ON state. The switch flips between the ON and OFF state with k_ON and k_OFF rate constants. etting the k_OFF constant to zero results in the switch remaining always in the ON state. Product is degraded with rate constant k_deg.

Program output
--------------
There are four output files created:
output*.par  - contains all simulation parameters, interarrival statistics of the product, statistics of the number of initiations per ON state, and time averages of all the species,
output*.prel - intervals between consecutive product releases (1 column),
output*.bl   - numbers of initiations for every on state,
output*.spec - time average of the ON state and the product P.

