//--------------------------------------------------------------------------------------------------------------------------------------------------//
//                                                              SABRE Sensitivity                                                                   //
//                                                                                                                                                  //
//                                          Using RooFit to fit various DM models to DAMA/LIBRA data                                                //
//                                                                                                                                                  //
//                                                  Created by Madeleine Zurowski on 26/6/19.                                                       //
//                                                                                                                                                  //
//--------------------------------------------------------------------------------------------------------------------------------------------------//

#ifndef __CINT__
//#include "RooGlobalFunc.h"
#endif
#include "headers/models.h"


using namespace RooFit ;
using namespace RooStats;

//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                              D e f i n e   r e l e v a n t   c o n s t a n t s
//--------------------------------------------------------------------------------------------------------------------------------------------------

// Define the interaction rate for both targets based on theory for a particular detector
//---------------------------------------------------------------------------------------
double N_T(double A) {return (6.02E26)/A;}
const double c0_conversion = (c*1.E-8)*pow(10.,-28.)*86400*(197.327)*(197.327);
double infty = RooNumber::infinity();

// Set DM global variables, where m_x and sig_x can be cycled through.
// This will allow us to write the signal expressions only in terms of E and t.
// This is needed for model dependent cases.
const double DAMAoffset = 152.2;                    // DM offset observed by DAMA
double m_x;
double sig_x;
double d_x = 0;                               // mass splitting, only required for inelastic models

// year and Rb limits (units of years and cpd/kg/keV)
const double Rb_min = 0.5;
const double Rb_max = 5.5;
const double Rb_step = 0.1;
int Nbin_min = 6;
const int Nbin_max = 96;
const int Nbin_step = 6;

//SI model
const double min_mass = 0;
const double max_mass = 3.;
const double mstep = 0.05;
const double min_sig = -42;
const double max_sig = -38;
const double sigstep = 0.01;


//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                              D e f i n e   i n t e r a c t i o n   r a t e
//--------------------------------------------------------------------------------------------------------------------------------------------------
//  To compute the sensitivity we need the interaction rate of DM at our detector. For the model indep case we take those visible at D/L
//  Otherwise, use values computed using the form factors etc. constructed using my packages.

//----------------------------------
//             DAMA Rate
//----------------------------------

// construct the number of events seen in units of cpd/kg
// Should read modulation in from D/L data
double Rm_DAMA(double Er){
  TGraph *damadata = new TGraph("dama_data.txt") ;
  int start = 2*(Er-E_width-1);
  int finish = 2*(Er+E_width-1);
  double tot=0;
  for (int i = start; i < finish; i++) {
    // multiply the value in each bin by the bin width
    tot+=0.5*damadata->GetY()[i];
  }
  delete damadata;
  return tot*Threshold(Er);
  //return 0.00261*2*E_width;
}

// Compute R0
// using limits given by DAMA in their perspective paper
double R0_DAMA(double Er){
  int start = 2*(Er-E_width-1);
  int finish = 2*(Er+E_width-1);
  double tot=0;
  for (int i = start; i < finish; i++) {
    // multiply the value in each bin by the bin width
    double e = start+i;
    if(e<=2){tot+=0.8*1;}
    if(e>2 && e<3){tot+=0.24*1;}
    else{tot+=0.12*1;}
  }
  return tot*Threshold(Er);
  // return 0.0;
}

//----------------------------------
//    Model Dep Rate
//----------------------------------

double NadRdE(double E, double Eee, double m, double sig, double d){
    return N_T(ANa)*DerivQuench(ANa, Eee)*c0_squared(sig,m)*c0_conversion*fm(Quench(ANa, Eee), m, d, ANa, jx)*(RhoDM/m)*ResFunction(E,Eee, ResFactor)*Threshold(E);
}
double NadRdE0(double E, double Eee, double m, double sig, double d){
    return N_T(ANa)*DerivQuench(ANa, Eee)*c0_squared(sig,m)*c0_conversion*f0(Quench(ANa, Eee), m, d, ANa, jx)*(RhoDM/m)*ResFunction(E,Eee, ResFactor)*Threshold(E);
}
double IdRdE(double E, double Eee, double m, double sig, double d){
    return N_T(AI)*DerivQuench(AI, Eee)*c0_squared(sig,m)*c0_conversion*fm(Quench(AI, Eee), m, d, AI, jx)*(RhoDM/m)*ResFunction(E,Eee, ResFactor)*Threshold(E);
}
double IdRdE0(double E, double Eee, double m, double sig, double d){
    return N_T(AI)*DerivQuench(AI, Eee)*c0_squared(sig,m)*c0_conversion*f0(Quench(AI, Eee), m, d, AI, jx)*(RhoDM/m)*ResFunction(E,Eee, ResFactor)*Threshold(E);
}

// Number of modulation events seen across energy bin (cpd/kg)
double Rm(double E){
    TF2 *NaMod = new TF2("NaMod", "NadRdE(x, y, [0], [1], [2])", 1E-1, 2*max_energy, 1E-1, infty);
    TF2 *IMod = new TF2("IMod", "IdRdE(x, y, [0], [1], [2])", 1E-1, 2*max_energy, 1E-1, infty);
    NaMod->SetParameter(0, m_x);
    NaMod->SetParameter(1, sig_x);
    NaMod->SetParameter(2, d_x);
    IMod->SetParameter(0, m_x);
    IMod->SetParameter(1, sig_x);
    IMod->SetParameter(2, d_x);
    double na = (ANa/(ANa+AI))*NaMod->TF2::Integral(E-E_width, E+E_width, 1E-1, 2*max_energy, 1E-7);
    double i = (AI/(ANa+AI))*IMod->TF2::Integral(E-E_width, E+E_width, 1E-1, 2*max_energy, 1E-7) ;
    if(i!=i){
      cout<<"integral zero"<<endl;
      i=0;
    }
    if(na!=na){
      cout<<"na integral zero"<<endl;
      na = 0;
    }
    delete NaMod;
    delete IMod;
    return na+i;
}

double R0(double E){
    TF2 *NaAv = new TF2("NaAv", "NadRdE0(x, y, [0], [1], [2])", 1E-1, 2*max_energy, 1E-1, infty);
    TF2 *IAv = new TF2("IAv", "IdRdE0(x, y, [0], [1], [2])", 1E-1, 2*max_energy, 1E-1, infty);
    NaAv->SetParameter(0, m_x);
    NaAv->SetParameter(1, sig_x);
    NaAv->SetParameter(2, d_x);
    IAv->SetParameter(0, m_x);
    IAv->SetParameter(1, sig_x);
    IAv->SetParameter(2, d_x);
    double na = (ANa/(ANa+AI))*NaAv->TF2::Integral(E-E_width, E+E_width, 1E-1,max_energy, 1E-7);
    double i = (AI/(ANa+AI))*IAv->TF2::Integral(E-E_width, E+E_width, 1E-1, max_energy, 1E-7) ;
    if(i!=i){
      cout<<"integral zero"<<endl;
      i=0;
    }
    if(na!=na){
      cout<<"na integral zero"<<endl;
      na = 0;
    }
    delete NaAv;
    delete IAv;
    return na+i;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                              D e f i n e   s t a t i s t i c a l   f u n c t i o n s
//--------------------------------------------------------------------------------------------------------------------------------------------------

// function for the expected DM interaction rate
double cosine_function(double x, double a, double b, double x0){
  //return a +b*TMath::Cos(2.*TMath::Pi()*(x+91.2-x0)/365.25);
  return a +b*TMath::Cos(2.*TMath::Pi()*(x-x0)/365.25);
}

// compute/generate the rate of fluctuations Rf for some signal and background case
// Takes as an input the number of bin periods the detector has run for, and functions for the signal and background rates in cpd/kg
// event(Er,t) will depend on some chosen signal/background model, so also has units of cpd/kg
// These must be functions of energy, and time, in that order.
double Rf(int Nbins, double Er, double RM, double RO, bool bkg = true){
    // compute limits of the plot
    double ymax = exposure*binperiod*2*E_width;//*binned_background(Er); // assume that the background or signal (even combined) won't exceed 10
    double xmin = -2*pow(10*ymax,0.5);
    double xmax = 2*pow(10*ymax,0.5);
    // construct roovariable
    RooRealVar t("t","t",0,Nbins*binperiod);
    RooRealVar yval("y","y",0,ymax*10+0.5*xmax) ;
    RooDataSet *dh = new RooDataSet("dh", "dh", RooArgSet(t,yval), StoreError(RooArgSet(t,yval)));
    // construct fitted signal function (cosine)
    RooRealVar *cons = new RooRealVar("cons", "cons", 0, ymax*10+xmax);
    RooRealVar *amp = new RooRealVar("amp", "amp",7*xmin, 8*xmax);
    //RooRealVar offset("offset", "offset",0,370);        // for fitting of offset
    RooRealVar offset("offset", "offset",DAMAoffset);     // for setting to DM mod
    RooAbsReal *fit_model = bindFunction("fit_signal", cosine_function,t, *cons, *amp, offset) ;
    RooRealVar coeff = RooRealVar("coeff","coeff",1);
  

    TRandom2 *a = new TRandom2();
    a->SetSeed();
    double X,Y,XE,YE;
    int n=0;
    for(n=0; n<Nbins; n++){
        X = n*binperiod;
        XE = 0.5*binperiod;
        //generate fluctuations using a Poissonian distribution
        if(bkg==true){Y = (a->Poisson(exposure*binperiod*(background_model(Er,n*binperiod)+pid_alpha*cosine_function(n*binperiod,RO,RM,DAMAoffset))));}
        else{Y = (a->Poisson(exposure*binperiod*cosine_function(n*binperiod,RO,RM,DAMAoffset)));}
        YE = pow(Y,0.5);

        t = X;
        t.setError(XE);
        yval = Y;
        yval.setError(YE);
        dh->add(RooArgSet(t,yval));
    }
    delete a;

    // compute fit
    fit_model->chi2FitTo(*dh,YVar(yval));
    double output_amp = coeff.getVal()*amp->getVal();

    delete dh;
    delete cons;
    delete amp;
    delete fit_model;
    return output_amp;
}

// construct the pdfs to compare for the sensitivity
// For a CL exlcusion plot, CL=true, and for m-sig CL=false
double test_sensitivity(int Nbins, double(*mod)(double), double(*cons)(double), bool CL){
  int points = (max_energy-min_energy)/(2.*E_width);
  double cl2 = 0;
  double Er;

  for (int i =0; i <points; i++) {
    Er = min_energy+E_width +2*E_width*i;
    // compute the expected signal in this bin (cpd/kg)
    double RM = mod(Er);
    double RO = cons(Er);

    double ymax = exposure*binperiod*(background_model(Er,0)+RO+RM);
    //double max = ymax*2;
    double mean_max = pow(ymax,0.5);

    // select sensible x values based on approximate modulation in that bin
    double xmin = -2*mean_max;
    double xmax = 3*mean_max;
    RooRealVar x("x","x",xmin,xmax);
    double min_std = 2*(xmax-xmin)/100;
    TH1F *amp_collection_background = new TH1F("amp_coll_b","amp_b",100,xmin,xmax);
    TH1F *amp_collection_signal_background = new TH1F("amp_coll_sb","amp_sb",100,xmin,xmax);

    // fill histograms
    //double test_b, test_sb,progress;
    for(int i=0;i<fluc_tests;i++){
      amp_collection_background->Fill(Rf(Nbins, Er, 0,0));
      amp_collection_signal_background->Fill(Rf(Nbins, Er, RM,RO));
    }

    //background
    RooRealVar *sigma_b = new RooRealVar("sigma_b","sigma_b",min_std,mean_max);
    RooRealVar *mean_b = new RooRealVar("mean_b", "mean_b",-mean_max,3*mean_max);
    RooGaussian *gauss_b = new RooGaussian("gauss_b","gauss_b",x,*mean_b,*sigma_b) ;
    RooDataHist *dh_b = new RooDataHist("dh_b", "dh_b", x, Import(*amp_collection_background));
    delete amp_collection_background;

    // //signal + background
    RooRealVar *sigma_sb = new RooRealVar("sigma_sb","sigma_sb",min_std,mean_max);
    RooRealVar *mean_sb = new RooRealVar("mean_sb", "mean_sb",-mean_max,3*mean_max);
    RooGaussian *gauss_sb = new RooGaussian("gauss_sb","gauss_sb",x,*mean_sb,*sigma_sb) ;
    RooDataHist *dh_sb = new RooDataHist("dh_sb", "dh_sb", x, Import(*amp_collection_signal_background));
    delete amp_collection_signal_background;
    // //----------------------------------------------------------------------------
    // fit it
    // background only
    gauss_b->fitTo(*dh_b);
    // signal+background
    gauss_sb->fitTo(*dh_sb);

    // for exclusion plots use sigma_sb, for discovery sigma_b.
    double diff_test;
    if(exclusion==true){
      diff_test = pow(fabs(mean_sb->getVal() - mean_b->getVal())/sigma_sb->getVal(),2.);
    }
    else{diff_test = pow(fabs(mean_sb->getVal() - mean_b->getVal())/sigma_b->getVal(),2.);}
    cl2+=diff_test;
    delete dh_b;
    delete dh_sb;
    delete gauss_b;
    delete gauss_sb;
    delete mean_sb;
    delete mean_b;
    delete sigma_sb;
    delete sigma_b;
  }

  if(CL==true){
    // return  the CL as a multiple of signal+background uncertainty
    return pow(cl2,0.5);
  }

  else{
    // compute p value
    return 0.5*(1-TMath::Erf(pow(cl2/2,0.5)));
  }

}


//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                             C a l c u l a t e   s e n s i t i v i t y
//--------------------------------------------------------------------------------------------------------------------------------------------------
//                  Confidence level plots as function of detector
//

// Function to obtain a plot of the model indepedent exclusion power
// I.e., confidence level of excluding DAMA assuming null results
void model_indep(){
  cout<<"Model independent calculation. DM signal based on DAMA data."<<endl;
  cout<<"Output sensitivity plot will be saved as "<<output_file<<endl;
  cout<<"Assuming active mass of "<<exposure<<" kg and background "<<background_model(min_energy+E_width, 0)/(2*E_width)<<" cpd/kg/keV"<<endl;
  int run_bins = 30;
  int count = 0;
  double Nbin;
  double rb;
  double test_cl;
  double ME;
  double progress;
  double r0;
  double rm;

  int n_points = int((Nbin_max-Nbin_min)/Nbin_step);
  
  // arrays for making error bars
  double cl_array[cycles][n_points];
  double Rb_array[cycles][n_points];
  double list_cl[cycles];


  // arrays for plotting
  double median_cl[n_points];
  double min_cl[n_points];
  double max_cl[n_points];
  double nbin_plot[n_points];


  //changing n bin
  for(int cycle=0; cycle<cycles; cycle++){
    count = 0;
    for(Nbin = Nbin_min; Nbin < Nbin_max; Nbin = Nbin+Nbin_step){
      nbin_plot[count] = Nbin;
      cl_array[cycle][count] =test_sensitivity(Nbin,Rm_DAMA,R0_DAMA,true);
      count++;
      // counter
      int barWidth = 70;
      progress = (Nbin-Nbin_min)/(Nbin_max-Nbin_min-Nbin_step);
      double pos = barWidth*progress;
      std::cout<<"Progress for cycle "<<cycle+1<<" of "<<cycles<<": [";
      for (int i = 0; i <= barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else std::cout << " ";
      }
      std::cout << "] " << int(progress * 100.0) << " %\r";
      std::cout.flush();
    }
    std::cout << std::endl;
  }

  // TGraph for shading
  TGraph *errshade = new TGraph(2*n_points);
  
  for(int i=0;i<count;i++){
    for(int j = 0; j<cycles; j++){
      list_cl[j] = cl_array[j][i];
    }
    
    median_cl[i]= TMath::Mean(cycles, list_cl);
    min_cl[i]= median_cl[i] - TMath::StdDev(cycles, list_cl);
    max_cl[i]= median_cl[i] + TMath::StdDev(cycles, list_cl);

    errshade->SetPoint(i,nbin_plot[i],max_cl[i]);
    errshade->SetPoint(2*count-i-1,nbin_plot[i],min_cl[i]);
  }

  // P l o t
  // -------
  gStyle->SetPalette(kSolar);
  TGraph* cl_plot = new TGraph(n_points,nbin_plot,median_cl);
  TGraph* cl_min_err = new TGraph(n_points,nbin_plot,min_cl);
  TGraph* cl_max_err = new TGraph(n_points,nbin_plot,max_cl);
  TCanvas *c1 = new TCanvas("c1","SABRE Sensitivity");
  //c1->SetLogy();
  //c1->SetLogx();
  c1->SetGrid();
  TH1F *cl_plotframe = c1->DrawFrame(0, 0,100,10);
  //TH1F *cl_plotframe = c1->DrawFrame(0, 0,150,4);
  cl_plotframe->SetXTitle("Months");
  cl_plotframe->SetYTitle("#sigma C.L.");
  cl_plotframe->SetTitle("Sensitivity for 50 kg, 36 bins");


  if(cycles!=1){
    errshade->SetFillStyle(3013);
    errshade->SetFillColor(2);
    errshade->Draw("f");
    cl_min_err->SetLineColor(2);
    cl_min_err->Draw("l");
    cl_max_err->SetLineColor(2);
    cl_max_err->Draw("l");
  }

  // include the sensitivity
  cl_plot->Draw("l");

  // // S a v e  f i l e
  // // -----------------
  TFile *f = new TFile(output_file,"RECREATE");
  cl_plot->Write("MedianSensitivityGraph");
  c1->Write("Sense_Plot");
    if(cycles!=1){
    cl_min_err->Write("MinSensitivityError");
    cl_max_err->Write("MaxSensitivityError");
    errshade->Write("ErrorShading");
  }
  f->Close();
}

// Function to plot traditional m-sigma exclusion plots
// Requires assumption of some model
void model_dep(int model_setting){
  cout<<"Model dependent calculation. Model selected:"<<endl;
  // S e t u p   m o d e l
  // ---------------------
  // Select a case to investigate, that is, a set of coupling constants from models.h
  model(model_setting);
  cout<<"Output sensitivity plot will be saved as "<<output_file<<endl;
  cout<<"Assuming active mass of "<<exposure<<" kg and background "<<background_model(min_energy+E_width, 0)/(2*E_width)<<" cpd/kg/keV"<<endl;


  // construct the expected uncertainty as a function of rate over the specified time period
    TGraph* sigma_r0 = new TGraph();
    int i = 0;
    ofstream myfile;
    if(read_corres==false){
        myfile.open (corres_file);
        for(double r0=0.05;r0<10.0;r0=r0+0.01){
            // compute the expected signal in this bin (cpd/kg)
            double ymax = exposure*binperiod*r0;
            //double max = ymax*2;
            double mean_max = pow(ymax,0.5);

            // select sensible x values based on approximate modulation in that bin
            double xmin = -2*mean_max;
            double xmax = 3*mean_max;
            RooRealVar x("x","x",xmin,xmax);
            double min_std = 2*(xmax-xmin)/100;
            TH1F *amp_collection_signal_background = new TH1F("amp_coll_sb","amp_sb",100,xmin,xmax);

            // fill histograms
            //double test_b, test_sb,progress;
            for(int i=0;i<2000;i++){
                amp_collection_signal_background->Fill(Rf(Nyears*12, 4, 0.0, r0,false)); 
            }

            // //signal + background
            RooRealVar *sigma_sb = new RooRealVar("sigma_sb","sigma_sb",min_std,mean_max);
            RooRealVar *mean_sb = new RooRealVar("mean_sb", "mean_sb",-mean_max,3*mean_max);
            RooGaussian *gauss_sb = new RooGaussian("gauss_sb","gauss_sb",x,*mean_sb,*sigma_sb) ;
            RooDataHist *dh_sb = new RooDataHist("dh_sb", "dh_sb", x, Import(*amp_collection_signal_background));
            delete amp_collection_signal_background;
            // //----------------------------------------------------------------------------
            // fit it and save into graph
            gauss_sb->fitTo(*dh_sb);
            sigma_r0->SetPoint(i, r0, sigma_sb->getVal());
            if(write_corres==true){
                myfile<<r0<<"   "<<sigma_sb->getVal()<<"\n";
            }
            delete dh_sb;
            delete gauss_sb;
            delete mean_sb;
            delete sigma_sb;
            i++;
        }
        myfile.close();
    }
    else{
        //sigma_r0 = new TGraph("text.csv");
        //sigma_r0 = new TGraph("anais_8yr_correspondence_nomod.txt");
        sigma_r0 = new TGraph(corres_file);
    }
    
    cout<<"Correspondance computed, now testing sensitivity"<<endl;

    // C o n s t r u c t   g r a p h
    // -----------------------------
    int n_points = int((max_mass-min_mass)/mstep)+1;
    double mass_array[n_points];
    double sig_array[n_points];
    double list_sig[cycles];

    double sigma_b = sigma_r0->Eval(0);

    // arrays for plotting
    double median_sig[n_points];
    double min_sig_err[n_points];
    double max_sig_err[n_points];

    int count = 0;
    double progress;
    double delta = model_delta;
    double Mstep=mstep;
    double last_cl;
    for(int cycle = 0; cycle<cycles;cycle++){
        count = 0;
        for(double m = min_mass; m <= max_mass+mstep; m = m+Mstep){
            double teststat;
            double sigma;
            for (sigma = min_sig; sigma<=max_sig; sigma=sigma+sigstep) {
                m_x = pow(10,m+9);
                sig_x = pow(10,sigma);
                double test_stat = 0;
                for(double E=min_energy+E_width; E<max_energy;E +=2*E_width){
                    double rate = pid_alpha*R0(E)+background_model(E,0); 
                    double mod_rate = Rm(E)*exposure*binperiod;
                    double sigma_sb = sigma_r0->Eval(rate);
                    test_stat += pow(fabs(mod_rate)/sigma_sb,2);
                }
                double pval = 0.5*(1-TMath::Erf(pow(test_stat/2,0.5)));
                if(pval<=benchmark_pval){
                    break;
                }
            }
            mass_array[count] = pow(10,m);
            sig_array[count] = pow(10,sigma)/4; //factor 4 to transform from isospin to nucleon xsec
            count++;
            // counter
            int barWidth = 70;
            progress = (m-min_mass)/(max_mass-min_mass);
            double pos = barWidth*progress;
            std::cout<<"Progress for cycle "<<cycle+1<<" of "<<cycles<<": [";
            for (int i = 0; i <= barWidth; ++i) {
                if (i < pos) std::cout << "=";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
        }
    std::cout << std::endl;
    }

  
  // P l o t
  // -------
  gStyle->SetPalette(kSolar);
  TGraph* sense = new TGraph(n_points,mass_array,sig_array);
  TCanvas *c1 = new TCanvas("c1","SABRE Sensitivity");
  c1->SetLogy();
  c1->SetLogx();
  c1->SetGrid();
  TH1F *senseframe = c1->DrawFrame(1, pow(10,min_sig),pow(10,max_mass),pow(10,max_sig));
  senseframe->SetXTitle("DM Mass (GeV)");
  senseframe->SetYTitle("#sigma (cm^{2})");
  senseframe->SetTitle( Form("SABRE Sensitivity after %f years",Nyears));

  // include the sensitivity
  sense->Draw("l");

  // S a v e  f i l e
  //-----------------
  TFile *f = new TFile(output_file,"RECREATE");
  sense->Write("sense");
  c1->Write("Sense_Plot");    
}

void Sensitivity(){
  threshval=set_threshval;
  poiss = bool_poiss;
  setquenching(quenching);
  gErrorIgnoreLevel = kWarning;
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  RooMsgService::instance().setSilentMode(kTRUE);
  if(DM_model==0){model_indep();}
  else{model_dep(DM_model);}
}
