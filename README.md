# ENACode
ENA Code from Kornbleuth et al. (2018, 2020, 2023)

ENA model to generate synthetic ENA maps from an MHD solution

Currently, ENA model configured for MHD-kinetic solutions derived from the Space Weather Modeling Framework (SWMF; Toth et al. 2005)

The model uses input plasma and neutral solutions from the global heliosphere solution to evaluate charge exchange events along radial trajectories to determine the probability of generating an ENA of a particular energy. Ions are propagated along streamlines from the termination shock, where the different characteristic properties of each ion species is determined. Inflowing neutrals will charge exchange with these different ion species, which can generate ENAs at particular energies. The total number of ENAs along a particular radial line-of-sight is integrated, with the survival probability that an ENA will survive transit to the observer without charge exchanging again taken into account. 

The original framework for this model is from Prested et al. (2008). The first version of this model following an overhaul is presented in Kornbleuth et al. (2018), which uses a multi-fluid treatment of neutrals and a simplified version of the ENA model. A more refined version is presented in Kornbleuth et al. (2020), while improved modeling of the ions at the termination shock is presented in Kornbleuth et al. (2023). 

Input plasma file is stored in secondary_files/input/opherwithneutralsplasmafinish.dat

Input neutral file is stored in secondary_files/input/opherwithneutralsnhfinish.dat

Input boundary (i.e. region) file is stored in secondary_files/input/boundaries.txt

List of files in top directory:

distribution_functions.pro: procedure with available distribution functions used in flux calculations

ena_prod_trans.pro: calculates charge-exchange x-section and applies constants to integrand

flux_and_counts.pro: sets up plotting procedure by reading in flux output from code

get_limits.pro: reads in regions of the heliosphere from input file (boundaries.txt)

load_opher.pro: sets up different arrays used as structures for plasma and neutral solutions separately

load_parameters.pro: parameter file where flags are set for the code and input files are specified - always start with this!

make_countrate_files_batch: batch script for running ENA code (see below)

make_countrate_files.pro: includes procedures for reading in input files, calculating ENA flux for each cell, integrating flux, and generating output - core of code!

make_plots_batch: batch script for plotting ENA maps (see below)

plotter_eps.pro: procedure used to plot ENA flux in Mollweide projections

ps_off.pro: procedure used in plotting ENA maps

ps_on.pro: procedure used in plotting ENA maps

read_stream_g.pro: calculates extinction along streamlines, propagates ion ratios from termination shock along streamlines, interpolates to ENA code grid

To run the ENA code from the top directory, in the IDL command line use "@make_countrate_files_batch". For plotting the resultant ENA maps, in the IDL command line use "@make_plots_batch"

Output ENA data files will be stored in secondary_files/output/

Output ENA maps following plotting will be stored in secondary_files/output/maps/

Output MHD/ENA data files for plotting in Tecplot (or other visualization software) will be stored in secondary_files/output/tec (view README in secondary_files/output/tec for processing information)

Need the Coyote IDL library for proper usage (http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD)
