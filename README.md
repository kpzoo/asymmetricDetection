# asymmetricDetection
Source code underlying the paper: 
Asymmetric limits on timely interventions from noisy epidemic data by Kris Parag, Ben Lambert, Christl Donnelly and Sandor Beregi.
Preprint: https://www.medrxiv.org/content/10.1101/2025.01.28.25321278v1

#### Matlab files

System Requirements:

Fig1.m and Fig2.m tested on Matlab v2023a and above and macOS 13
Slight dependence on the linespecer package (license included in main folder).


Instructions and installation

Run FigX.m where X is the figure in the manuscript to be reproduced.
Run times of all scripts are of the order of minutes or faster.

#### R files

To run the simulations use the following files in the R folder.

To test first lockdown deployment times:

EpiCont_parallel_covid.R
EpiCont_parallel_ebola.R
	
To test first lockdown relaxation times:
	
EpiCont_parallel_relax_covid.R
EpiCont_parallel_relax_ebola.R

To run simulations with surveillance noise use the settings delay = 1, ur = 1 when calling the function Epi_MPC_run_wd
and delay = 0, ur = 0 otherwise.


