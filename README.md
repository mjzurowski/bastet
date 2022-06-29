This repository contains the files needed to run sensitivity tests for the SABRE experiment.

---

## Header files

The following files are needed for the fits/tests to run:

1. detector.h - contains the threshold, quenching factor, resolution, and loads background for the detector. At present, this only contains info for SABRE.
2. formfactors.h - contains the operators used to describe DM interactions. Expressions are taken from arxiv 1308.6288 and 1203.3542. TODO: instrument more operators
3. model.h - contains the collection of coupling constants that correspond to the models of interest.
4. veldist.h - loads the velocity distributions chosen from the velocity_distribution folder (done in config.h), and combines them with the appropriate form factors. See section 3.3 of arxiv 2005.10404 for more info as to how this is done.
5. config.h - this is file that should be edited for a users special case. Information on all of the variables that can be set and what they are can be found in this file (note that this is not saved in the headers folder).

Note that other velocity distributions can be included here if so desired, by saving them as a function of minimum velocity in a .dat or .txt file, and including them in veldist.h


---

## Model fitting

ModelFittingDM.c is able to compute fits for DM mass and cross section, as well as the mass splitting for an inelastic interaction, to the data given in data_dama.txt. To do this, input a set of normalised coupling constants into the model.h header, and call it by running DAMA_Fit() with the model number as the argument.
Note that the quenching factor and/or the velocity distribution used for the fit can be changed in the main code.

TODO: configure this to allow allow fitting of the coupling constants.

---

## Sensitivity plots

Sensitivity.c is able to compute the sensitivity of SABRE in a number of different ways. This is done by loading and running the Sensitivity.c script in ROOT, which will call either a model depedent or model independent method depending on the setting in config.h.

To test the influence of different energy bins, E_width, max_energy, and min_energy can all be changed by the user in the config file, though note that decreasing the energy bin width increases the run time/memory requirements of the code. The background model can also be changed in config.h. This defines functions defined in detector.h as a function of both energy and time. Again, note that more complex backgrounds will take longer to run.

The sensitivity computation/methodology is outlined in arxiv 2005.10404, but follows these basic steps:

1. Generate the fluctuations in data over the detector lifetime observed for both background only and signal+background models.
2. Fit the resulting time spectrum to a background+cosine model, and save the cosine amplitude.
3. Repeat steps 1 and 2 100 times, and use the amplitudes to construct gaussian distributions for both models.
4. Compute the probability of observing the background mean under the assumption of the signal+background distribution. This give the confidence level of signal exclusion, as a multiple of the signal+background uncertainty, so |μ_sb−μ_b|/σ_sb.
5. Where a p-value is preferred over a CL, this can be computed by 0.5*(1-TMath::Erf(CL/sqrt(2)))
          
