Written by Maciej Dobrzynski
Centre for Mathematics and Computer Science (CWI), Amsterdam, The Netherlands.

Gillespie simulation to analyze interarrival times of motor proteins progressing along the biopolymer chain.

Type 'make' to compile, and './slide.x' to execute the program. 
Help message is displayed if no command-line parameters are given.
Type 'make clean' to remove compiled files.

  Reactions:
  Switch:     OFF <-> ON
  Initiation: ON -> ON + a1
  Elongation: a1 -> a2 -> ... -> aN
  Product release: aN -> P
  Product degradation: P -> null  
  + pausing (at every site )
  a1 <-> a1inactive

Initiation
----------
Motors initiate their motion with the initiation rate constant k_INI if the switch is in the ON state. Switch flips between the OFF and the ON state with rate constants k_OPEN and k_CLOSE, respectively. Setting the OFF rate constant, k_CLOSE, to zero keeps the switch in the ON state. 
          
Motor progression
-----------------
Motors jump forward to the site on the biopolymer chain with the rate constant k_ELONG. Motors can occupy only one site at the time during progression if BUMP parameter is set to 1. BUMP equal zero allows many motor proteins to occupy a single site. In this case motors can pass each other and collisions between motors do not occur.

Motor pausing
-------------
Motors can pause their progression at every site on the biolpolymer. The rate constant for this reaction is given by stallF. Returning from the paused state is described bythe rate constant stallB.

Termination
-----------
Motor proteins leave the chain with the termination rate constant k_TER. In this reaction a product P is formed. Product degradation occurs with the rate constant k_DEG.

Program output
--------------
There are three output files created:
output*.par  - contains all simulation parameters, occupation statistics of the polymer, and the mean and the standard deviation of the product interarrival times,
output*.aN   - intervals between consecutive product releases (1 column),
output*.mrna - level of product P in time. 

