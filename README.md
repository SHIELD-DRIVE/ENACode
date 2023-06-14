# ENACode
ENA Code from Kornbleuth et al. (2018, 2020, 2023)

ENA model to generate synthetic ENA maps from an MHD solution

Currently, ENA model configured for MHD-kinetic solutions derived from the Space Weather Modeling Framework (SWMF; Toth et al. 2005)

The model uses input plasma and neutral solutions from the global heliosphere solution to evaluate charge exchange events along radial trajectories to determine the probability of generating an ENA of a particular energy. Ions are propagated along streamlines from the termination shock, where the different characteristic properties of each ion species is determined. Inflowing neutrals will charge exchange with these different ion species, which can generate ENAs at particular energies. The total number of ENAs along a particular radial line-of-sight is integrated, with the survival probability that an ENA will survive transit to the observer without charge exchanging again taken into account. 

The original framework for this model is from Prested et al. (2008). The first version of this model following an overhaul is presented in Kornbleuth et al. (2018), which uses a multi-fluid treatment of neutrals and a simplified version of the ENA model. A more refined version is presented in Kornbleuth et al. (2020), while improved modeling of the ions at the termination shock is presented in Kornbleuth et al. (2023). 
