/*
  Written by Maciej Dobrzynski
  Email: m.dobrzynski@cwi.nl
  Centre for Mathematics and Computer Science (CWI), Amsterdam, The Netherlands.

  Model of motor protein sliding with pausing.
  Reactions:
  Switch:     OFF <-> ON
  Initiation: ON -> ON + a1
  Elongation: a1 -> a2 -> ... -> aN
  Product release: aN -> P
  Product degradation: P -> null  
  + pausing (at every site )
  a1 <-> a1inactive
  
  Notes:
  (1) If BUMPING is set to '1', only one motor protein can occupy a single site, i.e. elongation reaction is nonlinear: (1-si(n+1)).
  Otherwise, sites can be occupied by many motors (linear reactions)
  (2) A motor protein can occupy Npoly sites.

  Species:
  si[0..Nspec]
  si[0] - closed complex
  si[1] - open complex
  si[2] - a1
  ...
  si[Nspec-2] - aN
  si[Nspec-1] - mRNA

  Reactions:
  cn[0..Nreact]
  cn[0] - closed -> open: k_cl_open
  cn[1] - open -> closed: k_open_cl
  cn[2] - open -> a1 (injection): k_open_a1
  cn[3] - a1 -> a2: k_a_b
  ...
  cn[Nreact-3] - aN-1 -> aN: k_a_b
  cn[Nreact-2] - aN -> mRNA (mRNA production): k_aN_mRNA
  cn[Nreact-1] - mRNA -> null (mRNA decay): k_mRNA_null
*/

#include "main.h"

int main (int argc, char** argv) {

  long seed, n, i, k, ii, itmp;
  int simID, BUMPING, ipaused;
  long Nreact, Ndna, Npoly, Nspec, MUreact, Nout;
  double r1, r2, r2a0, Tcurr, Ttot, Tout, tau_next, sumtmp, Tequil, ttemp, Tmeas;
  double a0, *an, *cn, *hn, *savg, *cavg, *sstd, *cstd, occ_avg, occ_std;
  int *dna, *dna2; // vectors of dna occupation
  double k_cl_open, k_open_cl, k_open_a1, k_a_b, k_aN_rna, k_rna_null, k_stallF, k_stallB;
  long **vni;
  long *si;

  double aLastTon, aLastTprev, aLastTonAVG, aLastTonVAR;
  long aLastN, aLastPrev, occ_tmp;

  // output files
  FILE *fmrna, *fpar, *falast;
  char filename[128];

  // variables for measuring time
  double clock_diff;
  clock_t clock_start, clock_end;
  time_t time_start, time_end;

  seed = 11;

  if (argc != 16) {
    printf("Wrong number of parameters! Invoke:\n");
    printf("%s simID BUMP Ndna Npoly Ttot Tout Tequil k_OPEN k_CLOSE k_INI k_ELONG k_TER k_DEG stallF stallB\n", argv[0]);
    printf("\nExplanation:\n");
    printf("simID   - number of the simulation added to the name of the output file,\n");
    printf("BUMP    - a flag taking value 1 or 0 depending whether motor collisions are accounted for or not,\n");
    printf("Ndna    - length of the polymer chain,\n");
    printf("Npoly   - number of sites occupied by a single motor protein,\n");
    printf("Ttot    - total simulation time,\n");
    printf("Tequil  - equilibration time,\n");
    printf("k_OPEN  - rate constant for OFF -> ON reaction,\n");
    printf("k_CLOSE - rate constant for ON -> OFF reaction,\n");
    printf("k_INI   - initiation rate constant: ON -> ON + a1,\n");
    printf("k_ELONG - elongation rate constant: a1 -> a2,\n");
    printf("k_TER   - termination rate constant: an -> P,\n");
    printf("k_DEG   - product degradation rate constant: P -> null,\n");
    printf("stallF  - stalling rate constant: ai -> ai_stalled,\n");
    printf("stallB  - return from stalling rate constant: ai_stalled -> ai,\n");
    printf("\nNotes:\n");
    printf("(1) Units of Ttot, Tequil and Tout are in Time\n");
    exit(1);
  }
  
  // command-line parameters
  simID = atoi(argv[1]); // simulation ID
  BUMPING = atoi(argv[2]); // whether sliding is with bumping
  Ndna = atoi(argv[3]); // Number of dna sites
  Npoly = atoi(argv[4]); // number of DNA sites occupied by RNApol
  Ttot = atof(argv[5]); // Total sim time
  Tout = atof(argv[6]); // Output interval
  Tequil = atof(argv[7]); // Equilibration
  k_cl_open = atof(argv[8]); // closed -> open rate
  k_open_cl = atof(argv[9]); // open -> closed rate
  k_open_a1 = atof(argv[10]); // open -> a1 (injection) rate
  k_a_b = atof(argv[11]); // ai -> ai+1 rate
  k_aN_rna = atof(argv[12]); // aN -> mRNA rate
  k_rna_null = atof(argv[13]); // mRNA -> null rate
  k_stallF = atof(argv[14]); // rate to pause
  k_stallB = atof(argv[15]); // rate to go to active state after pausing

  if (Ndna + 2 < Npoly) {
    printf("Error: the number of DNA sites is smaller than\n");
    printf("the number of sites ocupied by RNA polymerase!\n");
    exit(1);
  }

  if ((k_stallF > 0 || k_stallB > 0) && (BUMPING == 0)) {
    printf("Error: stalling wthout bumping doesn't make sense!\n");
    exit(1);
  }

  if (Npoly < 1) {
    printf("Error: The size of the polymerase cannot be lower than 1!\n");
    exit(1);
  }

	// opening output files
  sprintf (filename, "output%d.mrna", simID);
  fmrna = fopen (filename, "w");

  sprintf (filename, "output%d.par", simID);
  fpar = fopen (filename, "w");

  sprintf (filename, "output%d.aN", simID);
  falast = fopen (filename, "w");

	// setting the number of reactions and species in the simulation
  Nreact = Ndna + 2 + 1 + 1 + 2*Ndna;
  Nspec  = Ndna + 3 + Ndna;

	// writing the parameter output file
  fprintf(fpar, "Simulation id      : %d\n", simID);
  if (BUMPING)
    fprintf(fpar, "Elongation with bumping\n");
  else
    fprintf(fpar, "Elongation without bumping\n");

  fprintf(fpar, "Number of DNA sites: %d\n", Ndna);
  fprintf(fpar, "Size of the RNAp   : %d\n", Npoly);
  fprintf(fpar, "Number of reactions: %d\n", Nreact);
  fprintf(fpar, "Number of species  : %d\n", Nspec);
  fprintf(fpar, "Total sim. time    : %.4e\n", Ttot);
  fprintf(fpar, "Output frequency   : %.4e\n", Tout);
  fprintf(fpar, "Equilibration time : %.4e\n", Tequil);
  fprintf(fpar, "\nReaction rate constants:\n");
  fprintf(fpar, "Opening of the complex : %lf\n", k_cl_open);
  fprintf(fpar, "Closing of the complex : %lf\n", k_open_cl);
  fprintf(fpar, "Initiation rate        : %lf\n", k_open_a1);
  fprintf(fpar, "Elongation rate        : %lf\n", k_a_b);
  fprintf(fpar, "mRNA release rate      : %lf\n", k_aN_rna);
  fprintf(fpar, "mRNA degradation       : %lf\n", k_rna_null);
  fprintf(fpar, "Stalling rate          : %lf\n", k_stallF);
  fprintf(fpar, "Return from paused     : %lf\n", k_stallB);

  fflush(fpar);

	// allocating memory for species, reactions, stoichiometric matrix
  cn = (double*)malloc(Nreact*sizeof(double)); // vector: reaction rates
  cavg = (double*)malloc(Nreact*sizeof(double)); // vector: reaction rates - mean
  cstd = (double*)malloc(Nreact*sizeof(double)); // vector: reaction rates - STD

  an = (double*)malloc(Nreact*sizeof(double)); // vector: propensity functions
  hn = (double*)malloc(Nreact*sizeof(double)); // vector: combinations of reactants
  vni = (long**)malloc(Nreact*sizeof(long*));  // stochiometric matrix

  dna = (int*)malloc((Ndna+Npoly)*sizeof(int)); // vector: dna occupancy (0 or 1)
  dna2 = (int*)malloc((Ndna+Npoly)*sizeof(int)); // vector: dna occupancy (differes between paused and active)

  for (n = 0; n < Nreact; n++)
    vni[n] = (long*)malloc(Nspec*sizeof(long));

  si = (long*)malloc(Nspec*sizeof(long)); // vector: number of molecules of every species
  savg = (double*)malloc(Nspec*sizeof(double)); // vector: time average of number of molecules of every species
  sstd = (double*)malloc(Nspec*sizeof(double)); // vector: time average of numbe

  // initialising species
  for (i = 0; i < Nspec; i++) {
    si[i] = 0;
    savg[i] = 0.;
    sstd[i] = 0.;
  }

  for (i = 0; i < Nreact; i++) {
    cavg[i] = 0.;
    cstd[i] = 0.;
  }

  for (i = 0; i < Ndna+Npoly; i++)
    dna[i] = 0;

  si[1] = 1; // 1 open complex at the beginning

  // cn (rates)
  // closed -> open
  cn[0] = k_cl_open;

  // open -> closed
  cn[1] = k_open_cl;

  // open -> a1; warning! changed later in the loop (to include nonlinearity)
  cn[2] = k_open_a1;

  // elongation
  for (n = 3; n < (3 + Ndna - 1); n++)
    cn[n] = k_a_b;

  // aN -> mRNA
  cn[3+Ndna-1] = k_aN_rna;

  // stalling - forward
  for (n = (3 + Ndna); n < (3 + 2*Ndna); n++)
    cn[n] = k_stallF;

  // stalling - backward
  for (n = (3 + 2*Ndna); n < (3 + 3*Ndna); n++)
    cn[n] = k_stallB;

  // last reaction: mRNA decay
  cn[Nreact-1] = k_rna_null;

  // setting vni
  for (n = 0; n < Nreact; n++)
    for (i = 0; i < Nspec; i++)
      vni[n][i] = 0;

  // closed -> open
  vni[0][0] = -1;
  vni[0][1] = +1;

  // open -> closed
  vni[1][0] = +1;
  vni[1][1] = -1;

  // open -> open + a1
  vni[2][2] = +1;

  // a1 -> a2...
  itmp = 2;
  for (n = 3; n < (3 + Ndna - 1); n++) {
    vni[n][itmp] = -1;
    vni[n][itmp+1] = +1;
    itmp++;
  }

  // last an -> mRNA
  vni[3 + Ndna - 1][1 + Ndna] = -1;
  vni[3 + Ndna - 1][Nspec-1] = +1;

  // a1 -> a1inact...
  itmp = 2;
  for (n = (3 + Ndna); n < (3 + 2*Ndna); n++) {
    vni[n][itmp] = -1;
    vni[n][itmp+Ndna] = +1;

    itmp++;
  }

  // a1inact -> a1...
  itmp = 2;
  for (n = (3 + 2*Ndna); n < (3 + 3*Ndna); n++) {
    vni[n][itmp] = +1;
    vni[n][itmp+Ndna] = -1;

    itmp++;
  }
  
  // mRNA -> null
  vni[Nreact-1][Nspec-1] = -1;

  // Main loop
  Nout = 0;
  Tmeas = 0.;
  Tcurr = 0.;
  aLastN = 0;
  aLastPrev = 0;
  aLastTprev = 0.;
  aLastTonAVG = 0.;
  aLastTonVAR = 0.;

  occ_avg = 0.;
  occ_std = 0.;

  clock_start = clock();
  time_start = time(0);
  do {
    // setting hn
    // closed <-> open cplx
    hn[0] = si[0];
    hn[1] = si[1];

    // open -> open + a1
    hn[2] = si[1];
    
    // sliding: s1 -> s2, s2 -> s3...
    itmp = 2;
    for (n = 3; n < (3 + Ndna); n++) {
      hn[n] = si[itmp];
      itmp++;
    }

    // stalling: si_active -> si_inactive
    itmp = 2;
    for (n = (3 + Ndna); n < (3 + 2*Ndna); n++) {
      hn[n] = si[itmp];
      itmp++;
    }

    // stalling: si_inactive -> si_active
    itmp = 2+Ndna;
    for (n = (3 + 2*Ndna); n < (3 + 3*Ndna); n++) {
      hn[n] = si[itmp];
      itmp++;
    }
    
    // mRNA -> null
    hn[Nreact-1] = si[Nspec-1];
    
    // inhibition (bumping)
    // mrna at i cannot proceed if there's mrna at i+1 
    // (both active and stalled)

    if (BUMPING) {
      cn[2] = k_open_a1 * (1 - dna[0+Npoly-1]); 
      for (ii = 0; ii < Ndna; ii++) {
	cn[ii+3] = k_a_b * (1 - dna[ii+Npoly]);
      }
    }

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

    // statistics for the average and std occupancy
    // gathered BEFORE changing the state
    if (Tcurr > Tequil) {
      for (i = 0; i < Nspec; i++) {
	savg[i] += si[i]*tau_next;
	sstd[i] += si[i]*si[i]*tau_next;
      }
      
      Tmeas += tau_next;
      while (Nout*Tout < Tmeas) {
	// output mRNA level
	fprintf(fmrna, "%.6lf\t%.6lf\t%ld\n", Nout*Tout, Tcurr, si[Nspec-1]);

	
	occ_tmp = 0;
	// number of occupied sites
	for (i = 2; i < 12; i++) {
	  occ_tmp += si[i];
	}

	occ_avg += (double) occ_tmp / 10;
	occ_std += (double) occ_tmp * occ_tmp / 100;
	

	Nout++;
      };
    }
        
    
    // changing the state
    for (i = 0; i < Nspec; i++) {
      si[i] += vni[MUreact][i];
    }

    // resetting dna vector
    for (i = 0; i < Ndna+Npoly; i++) {
      dna[i] = 0;
      dna2[i] = 0;
    }
    
    // filling dna vector for occupancy
    for (i = 0; i < Ndna; i++) {
      ipaused = 0;
      if (si[i+2] > 0)
	ipaused = 0;
      else if (si[i+2+Ndna] > 0)
	ipaused = 1;
      
      if (si[i+2] > 0 || si[i+2+Ndna] > 0) {
	for (k = 0; k < Npoly; k++) {
	  dna[i+k] = 1;
	  
	  if (ipaused)
	    dna2[i+k] = 2;
	  else
	    dna2[i+k] = 1;	  
	}
      }
    }

    Tcurr += tau_next;

    // time between mRNA production
    // makes sense only for the model with bumping
    if (MUreact == Ndna+2) {
      if (aLastN > 1) {
	aLastTon = Tcurr - aLastTprev;
	fprintf(falast, "%.4e\n", aLastTon);
	aLastTonAVG += aLastTon;
	aLastTonVAR += aLastTon*aLastTon;
      }
      aLastN++;
      aLastTprev = Tcurr;
    }
    
  } while (Tcurr < Ttot);
  
  clock_end = clock();
  time_end = time(0);

  fprintf(fpar, "\nCPU time used (clock): %.4lf sec\n", clock_diff);
  fprintf(fpar, "CPU time used (time) : %.4lf sec\n", difftime(time_end, time_start));

  // output averages
  fprintf(fpar, "\nDNA %d\n", Ndna);
  fprintf(fpar, "Spec\tMean\t\tSTD\n");
  for (i = 0; i < Nspec; i++) {
    savg[i] /= Tmeas;
    sstd[i] /= Tmeas;
    sstd[i] = sqrt(sstd[i] - savg[i]*savg[i]);

    if (i == Nspec-1)
      fprintf(fpar, "%d\t%.6lf\t%.6lf\tmRNA\n", i, savg[i], sstd[i]);
    else if (i == 1)
      fprintf(fpar, "%d\t%.6lf\t%.6lf\tO\n", i, savg[i], sstd[i]);
    else if (i < Nspec-1-Ndna && i > 1)
      fprintf(fpar, "%d\t%.6lf\t%.6lf\ta%d\n", i, savg[i], sstd[i], i-1);
  }

  
  occ_avg /= Nout;
  occ_std /= Nout;

  fprintf(fpar, "\nNo. of occupied sited (scaled)\n");
  fprintf(fpar, "occAVG\t%.6lf\n", occ_avg);
  fprintf(fpar, "occSTD\t%.6lf\n", sqrt(occ_std - occ_avg*occ_avg));
  

  aLastTonAVG /= aLastN;
  aLastTonVAR /= aLastN;

  fprintf(fpar, "\nmRNA releases\t%ld\n", aLastN);
  fprintf(fpar, "mRNA release AVG\t%.6lf\n", aLastTonAVG);
  fprintf(fpar, "mRNA release STD\t%.6lf\n\n", sqrt(aLastTonVAR - aLastTonAVG*aLastTonAVG));

  free(cn);
  free(hn);
  free(an);
  free(si);
  free(savg);
  free(dna);
  free(dna2);
  free(vni);

  fclose(fmrna);
  fclose(fpar);
  fclose(falast);
  return(0);
}
