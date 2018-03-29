/*
  Written by Maciej Dobrzynski
  Email: m.dobrzynski@cwi.nl
  Centre for Mathematics and Computer Science (CWI), Amsterdam, The Netherlands.

  Gillespie simulation of the basic burst-generating model.
  The model defines the interrupted Poisson process.

  Reactions:
  Closed <-> Open
  Open -> Open + P
  P -> null
  
  Species:
  si[0..Nspec]
  si[0] - closed complex (init = 1)
  si[1] - open complex
  si[2] - P

  Reactions:
  cn[0..Nreact]
  cn[0] - closed -> open: k_on
  cn[1] - open -> closed: k_off
  cn[2] - open -> open + P: k_ini
  cn[Nreact-1] - P -> null (P decay): k_deg
*/

#include "main.h"

int main (int argc, char** argv) {

  long seed, n, i, ii, itmp;
  int simID;
  long Nreact, Nspec, Nstep, MUreact, Nout;
  double r1, r2, r2a0, Tcurr, Ttot, Tout, tau_next, sumtmp, Tequil, Tmeas;
  double a0, *an, *cn, *hn, *savg, *sstd;
  double k_on, k_off, k_ini, k_p_null;
  long **vni;
  long *si;

  double PrelT, PrelTprev, PrelTAVG, PrelTVAR;
  long PrelN, PrelTN;

  double NburstAVG, NburstVAR;
  long Nburst, NburstN;

  FILE *fspec, *fpar, *fprel, *fbl;
  char filename[128];

  double clock_diff;
  clock_t clock_start, clock_end;
  time_t time_start, time_end;

  seed = 5;

  if (argc != 9) {
    printf("wrong number of parameters!\n%s simID Ttot Tout Tequil k_ON k_OFF k_ini k_deg\n", argv[0]);
    printf("\nExplanation:\n");
    printf("simID  - number of the simulation added to the name of the output file,\n");
    printf("Ttot   - total simulation time,\n");
    printf("Tout   - output interval,\n");
    printf("Tequil - equilibration time,\n");
    printf("k_ON   - rate constant for opening of the switch,\n");
    printf("k_OFF  - rate constant for closing of the switch,\n");
    printf("k_ini  - rate constant for product P production,\n");
    printf("k_deg  - rate constant for product P degradation,\n");
    exit(1);
  }
  
  // command-line parameters
  simID = atoi(argv[1]); // simulation ID
  Ttot = atof(argv[2]); // Total sim time
  Tout = atof(argv[3]); // Output interval
  Tequil = atof(argv[4]); // Equilibration
  k_on = atof(argv[5]); // OFF -> ON rate
  k_off = atof(argv[6]); // ON -> OFF rate
  k_ini = atof(argv[7]); // ON -> P (initiation) rate
  k_p_null = atof(argv[8]); // P -> null rate

  // opening output files
  sprintf (filename, "output%d.spec", simID);
  fspec = fopen (filename, "w");

  sprintf (filename, "output%d.par", simID);
  fpar = fopen (filename, "w");

  sprintf (filename, "output%d.prel", simID);
  fprel = fopen (filename, "w");

  sprintf (filename, "output%d.bl", simID);
  fbl = fopen (filename, "w");

  // setting the number of reactions and species in the simulation
  Nreact = 4;
  Nspec  = 3;

  // writing the parameter output file
  fprintf(fpar, "Simulation id      : %d\n", simID);
  fprintf(fpar, "Number of reactions: %d\n", Nreact);
  fprintf(fpar, "Number of species  : %d\n", Nspec);
  fprintf(fpar, "Total sim. time    : %lf\n", Ttot);
  fprintf(fpar, "Output frequency   : %lf\n", Tout);
  fprintf(fpar, "Equilibration time : %lf\n", Tequil);

  // allocating memory for species, reactions, stoichiometric matrix
  cn = (double*)malloc(Nreact*sizeof(double)); // vector: reaction rates

  an = (double*)malloc(Nreact*sizeof(double)); // vector: propensity functions
  hn = (double*)malloc(Nreact*sizeof(double)); // vector: combinations of reactants

  vni = (long**)malloc(Nreact*sizeof(long*));  // stoichiometric matrix
  for (n = 0; n < Nreact; n++)
    vni[n] = (long*)malloc(Nspec*sizeof(long));

  si = (long*)malloc(Nspec*sizeof(long)); // vector: number of molecules of every species
  savg = (double*)malloc(Nspec*sizeof(double)); // vector: time average of number of molecules of every species
  sstd = (double*)malloc(Nspec*sizeof(double)); // vector: average squared

  // initialising species
  for (i = 0; i < Nspec; i++) {
    si[i] = 0;
    savg[i] = 0.;
    sstd[i] = 0.;
  }

  // initial condition
  si[1] = 1; // 1 open complex at Tini

  // cn (rates)
  // OFF -> ON
  cn[0] = k_on;

  // ON -> OFF
  cn[1] = k_off;

  // ON -> ON + P
  cn[2] = k_ini;

  // P -> null
  cn[3] = k_p_null;

  // setting vni
  for (n = 0; n < Nreact; n++)
    for (i = 0; i < Nspec; i++)
      vni[n][i] = 0;

  // OFF -> ON
  vni[0][0] = -1;
  vni[0][1] = +1;

  // ON -> OFF
  vni[1][0] = +1;
  vni[1][1] = -1;

  // ON -> ON + P
  vni[2][2] = +1;

  // P -> null
  vni[3][2] = -1;

  // Main loop
  Nstep = 0;
  Nout = 1;
  Tmeas = 0.;
  Tcurr = 0.;
  PrelN = 0;
  PrelTprev = 0.;
  PrelTAVG = 0.;
  PrelTVAR = 0.;
  PrelTN = 0;

  Nburst = 0;
  NburstAVG = 0.;
  NburstVAR = 0.;
  NburstN = 0;

  clock_start = clock();
  time_start = time(0);
  do {
    // setting hn
    // OFF <-> ON
    hn[0] = si[0];
    hn[1] = si[1];

    // ON -> ON + P
    hn[2] = si[1];

    // P -> null
    hn[3] = si[2];

    // computing propensity function(s)
    a0 = .0;
    for (n = 0; n < Nreact; n++) {
      an[n] = cn[n] * hn[n];
      a0 += an[n];
    }

    // drawing 2 uniform random numbers
    r1 = ran3(&seed);
    r2 = ran3(&seed);
    r2a0 = r2 * a0;
    
    // next time step
    tau_next = 1./a0 * log(1./r1);

    // next reaction
    MUreact = 0;
    sumtmp = an[0];
    while (sumtmp <= r2a0) {
      MUreact++;
      sumtmp += an[MUreact];
    };

    // output
    // computing time averages
    if (Tcurr > Tequil) {
      for (i = 0; i < Nspec; i++) {
	savg[i] += si[i]*tau_next;
	sstd[i] += si[i]*si[i]*tau_next;
      }

      Tmeas += tau_next;
    }
            
    // changing the state
    for (i = 0; i < Nspec; i++) {
      si[i] += vni[MUreact][i];
    }

    Tcurr += tau_next;
    Nstep++;

    // output of the state to a file at Tout intervals
    if (Tcurr > Tequil) {
      if (Tmeas > Nout*Tout) {
	fprintf(fspec, "%.6lf\t%d\t%d\n", Tcurr, si[1], si[Nspec-1]);
	Nout++;
      }
    }

    // interarrival times between P productions
    // output to a file
    if (MUreact == 2) {
      if (PrelN > 0) {
	PrelT = Tcurr - PrelTprev;
	fprintf(fprel, "%.4e\n", PrelT);
    	PrelTAVG += PrelT;
	PrelTVAR += PrelT*PrelT;
      }
      PrelN++;
      PrelTprev = Tcurr;
    }

    // length of a burst
    // (the amount of P produced during one opening)
    if (MUreact == 2) // P release
      Nburst++;
    else if (MUreact == 1) // closing of the complex
      {
	// write the number of initiation during the ON state to a file
	fprintf(fbl, "%d\n", Nburst);
	
	// accumulate averages
	NburstAVG += Nburst;
	NburstVAR += Nburst*Nburst;
	NburstN ++;
	Nburst = 0;
      }

  } while (Tcurr < Ttot);
  
  clock_end = clock();
  time_end = time(0);

  fprintf(fpar, "\nCPU time used (clock): %.4lf sec\n", clock_diff);
  fprintf(fpar, "CPU time used (time) : %.4lf sec\n", difftime(time_end, time_start));

  // output averages
  fprintf(fpar, "\nP interarrival times:\n");
  fprintf(fpar, "Trel AVG %.4e\n", PrelTAVG/PrelN);
  fprintf(fpar, "Trel STD %.4e\n", sqrt(PrelTVAR/PrelN - PrelTAVG*PrelTAVG/PrelN/PrelN));
  fprintf(fpar, "Trel N %d P releases\n\n", PrelN);

  fprintf(fpar, "Initiations per ON state:\n");
  fprintf(fpar, "Nburst AVG %.4e\n", NburstAVG/NburstN);
  fprintf(fpar, "Nburst STD %.4e\n", sqrt(NburstVAR/NburstN - NburstAVG*NburstAVG/NburstN/NburstN));
  fprintf(fpar, "Nburst N %d\n\n", NburstN);

  fprintf(fpar, "Species\tMean\tSTD\n");
  for (i = 0; i < Nspec; i++) {
    savg[i] /= Tmeas;
    sstd[i] /= Tmeas;
    sstd[i] = sqrt(sstd[i] - savg[i]*savg[i]);
    fprintf(fpar, "%d\t%.6lf\t%.6lf\n", i, savg[i], sstd[i]);

  }

  free(cn);
  free(hn);
  free(an);
  free(si);
  free(savg);
  free(vni);

  fclose(fspec);
  fclose(fpar);
  fclose(fprel);
  fclose(fbl);
  return(0);
}
