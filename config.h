#ifndef config_H
#define config_H
#endif 

// Defines a number of useful parameters and functions for sensitivity computation
// This is the file that should be edited by a user to find sensitivity for a given case.

#include "headers/detector.h"


// -----------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------- Statistical test settings ----------------------------------------------
// -----------------------------------------------------------------------------------------------------------------------

const int fluc_tests = 500;                                  // number of experimental lifetimes generated to construct gaussian
const int cycles = 1;                                        // cycles to produce error bars for model indep
bool exclusion = true;                                       // boolean setting for exclusion (true) vs. discovery (false) plots. 
const double benchmark_pval = 0.1;                           // p-value that dictates sensitivity level - (1-benchmark_pval)%
// for basic PID (for no PID, set both to 1)
const double pid_alpha = 1.0;                       // percentage of NR events (DM signal) passing selection
const double pid_beta = 1.0;                      // percentage of ER events (background) passing selection
const char output_file[] = "out_dep.root";              // name of output root file to be saved
const double set_threshval = 1.0;                                                // threshold vale for change from poissonian smearing to gaussian
bool bool_poiss = false ;                                                   // set true to use gaussian+poissonian smear, false for just gaussian
bool write_corres = false;                                   // boolean to write the correspondence func to txt file
bool read_corres = true;                                    // boolean to read correspondance from txt
const char corres_file[] = "sabresouth_correspondence_3yr.txt";              // name of output root file to be saved



// ---------------------------------------------------------------------------------------------------------------
// ---------------------------------------------- Detector settings ----------------------------------------------
// ---------------------------------------------------------------------------------------------------------------

double exposure = 50;                                            // exposure of the particular experiment we are performing a fit to (kg).
double Rb = 0.7228;                                         // Background value (if assumed constant) (cpd/kg/keV)
const double Nyears = 3;                                                 // Lifetime of experiment (only needed for model dependent testing)
double max_energy = 6.;                               // maximum energy value we want to test to
double min_energy = 1.;                               // minimum energy value of region of interest
double E_width = 0.25;                                 // half the width of the energy bin, so integral for a bin is performed from E-E_width to E+E_width
const double binperiod = 365.25/12;                         // time period of data taking - assume 1 month (days)

// Model setting
const int DM_model = 4;                                     // Select DM model from models.h. If==0, will use model independent method.

double background_model(double E, double t){                // define a background model in units of cpd/kg
  return pid_beta*Rb*2*E_width*Threshold(E);                                      //for assuming constant bkg of Rb in RoI
  //return 0;
}


// Quenching setting
// 0=DAMA quenching values, 1= Princeton quenching values, 2=University of Texas quenching values, 3=ANU
const int quenching = 0;

// Resolution function setting
// 0=DAMA, 1=COSINE, 2=SABRE PoP MC
int ResFactor = 0;

// ----------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------- Velocity distribution settings ----------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------------

// Load the required velocity distributions. Note that we need 4 separate files, a modulating and average component 
// for the "g" and "h" integrals (see arxiv 2005.10404 for more detail as to how these are implemented). Uncomment the
// desired distribution in config.h, or include new ones with the same structure (these must be placed in the velocity_distributions folder).

// Depending on which distribution will be used, slightly different constants will be used. 
// Simply comment out the unnecessary ones depending on user choice.
// Unchanging/common constants
const double vE = 232.8;       // solar velocity
const double vt = 30.;        // earth velocity
const double Pi = TMath::Pi();
const double sol_angle = TMath::Cos(TMath::Pi()/3.);       // Cos of Solar angle
const double sin_sol_angle = TMath::Sin(TMath::Pi()/3.);   // Sin of Solar angle


// - - - - - - - - - - - - - - - Maxwell Boltzmann dist. - - - - - - - - - - - - - - -
// constants:
const double vesc = 550.;     // escape velocity
const double v0 = 220.;       // dispersion velocity
const double RhoDM = 0.3;     // DM density in the galaxy
// plots:
TGraph *modh = new TGraph("velocity_distributions/modhMB.dat");
TGraph *modg = new TGraph("velocity_distributions/modgMB.dat");
TGraph *avh = new TGraph("velocity_distributions/avhMB.dat");
TGraph *avg = new TGraph("velocity_distributions/avgMB.dat");

// // - - - - - - - - - - - - - - - SHM+Stream dist. - - - - - - - - - - - - - - - 
// // constants:
// const double vesc = 520;
// const double v0 = 232.8;
// const double RhoDM = 0.5;
// // plots:
// TGraph *modh = new TGraph("velocity_distributions/modhStr.dat");
// TGraph *modg = new TGraph("velocity_distributions/modgStr.dat");
// TGraph *avh = new TGraph("velocity_distributions/avhStr.dat");
// TGraph *avg = new TGraph("velocity_distributions/avgStr.dat");