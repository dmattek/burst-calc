Supplementary code for M Dobrzyński & FJ Bruggeman, _Elongation dynamics shape bursty transcription and translation_, PNAS, 2009 ([link](https://doi.org/10.1073/pnas.0803507106)).

Analyze bursts in interarrival times.
Last modification: 30 March 2018

Authors: M Dobrzyński, Centre for Mathematics and Computer Science (CWI), Amsterdam, The Netherlands

Content of the directory:

- `burster`         - C code for Gillespie simulation of the basic burst-generating model
- `slide`           - C code for Gillespie simulation of motor protein progression along the polymer chain
- `matlab`          - Matlab scripts
	- `mdloghist.m`     - a Matlab script to produce histograms with logarithmic binning
	- `mdburstsz.m`     - a Matlab script to generate the burst size function
- `R`               - an Rstudio project
	- `calcBurstSzFn.R` - an R function to generate the burst size function
	- `analyzeBursts.Rmd` - an R notebook that demonstrates burst analysis using the example dataset
- `example-data-1`	- an example dataset generated with `./burster.x 1 10000 1 0 0.1 0.1 10 1`
