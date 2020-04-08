# Deltamethrin-Polymorph-Malaria
Code for Malaria Transmission Modeling in ‘A Faster-acting Deltamethrin Crystal Polymorph for Malaria Control’ by Yang et al.

Code Files:

DM.m Simulation code for malaria transmission model. Required by RunDMII20.m, RunDMII45.m, RunDMII70.m.

DMParameters.m Parameter file for simulation code. Required by RunDMII20.m, RunDMII45.m, RunDMII70.m.

mortality_ot.m Calculates the mortality of a population after contacting applied doses of DM Form II based on the mortality of a population exposed to the same amount of DM Form I. Required by RunDMII20.m, RunDMII45.m, RunDMII70.m, FeedersInfectedHuman_NoDM_DMI_DMII_20prevalence.m, FeedersInfectedHuman_NoDM_DMI_DMII_45prevalence.m.

brewermap.m The complete selection of ColorBrewer colorschemes (RGB colormaps). (c) 2014 Stephen Cobeldick Downloadable from https://github.com/spunt/bspmview/blob/master/supportfiles/brewermap.m Required by DMII20.m, RunDMII45.m, RunDMII70.m.

RunDMII20.m Calculates the effectiveness of DM Form I and II at 20% prevalence varying coverage (x-axis) and insecticide resistance (y-axis). Code for Figure S8A.

RunDMII45.m Calculates the effectiveness of DM Form I and II at 45% prevalence varying coverage (x-axis) and insecticide resistance (y-axis). Code for Figure 4A.

RunDMII70.m Calculates the effectiveness of DM Form I and II at 70% prevalence varying coverage (x-axis) and insecticide resistance (y-axis). Code for Figure S9.

FeedersInfectedHuman_NoDM_DMI_DMII_20prevalence.m Calculates effect of switching from Form I to Form II on the infectious mosquito population and human malaria infection dynamics at 20% prevalence.Code for Figure S8B.

FeedersInfectedHuman_NoDM_DMI_DMII_45prevalence.m Calculates effect of switching from Form I to Form II on the infectious mosquito population and human malaria infection dynamics at 45% prevalence.Code for Figure 4B.
