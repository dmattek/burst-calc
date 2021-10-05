Supplementary code for M Dobrzyński & FJ Bruggeman, _Elongation dynamics shape bursty transcription and translation_, PNAS, 2009 ([link](https://doi.org/10.1073/pnas.0803507106)).

Analyze bursts in interarrival times.
Last modification: 30 March 2018

Authors: M Dobrzyński, Centre for Mathematics and Computer Science (CWI), Amsterdam, The Netherlands

Content of the directory:
burster         - C code for Gillespie simulation of the basic burst-generating model,
slide           - C code for Gillespie simulation of motor protein progression along the polymer chain,
matlab          - Matlab scripts
	mdloghist.m     - Matlab script to produce histograms with logarithmic binning,
	mdburstsz.m     - Matlab script to generate the burst size function.
	
R               - Rstudio project
	calcBurstSzFn.R - R function to generate the burst size function.
	analyzeBursts.Rmd - R notebook demonstrating burst analysis using example dataset
example-data-1	- Example dataset generated with `./burster.x 1 10000 1 0 0.1 0.1 10 1`
