//--------------------------------------------------------------------------------------------------------------------------------------------------//
//                                                              DM Model Fitting                                                                    //
//                                                                                                                                                  //
//                                          Using RooFit to fit various DM models to DAMA/LIBRA data                                                //
//                                                                                                                                                  //
//                                                  Created by Madeleine Zurowski on 26/6/19.                                                       //
//                                                                                                                                                  //
//--------------------------------------------------------------------------------------------------------------------------------------------------//

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "headers/models.h"
#include "headers/detector.h"

using namespace RooFit ;
using namespace RooStats;

//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                              D e f i n e   i n t e r a c t i o n   r a t e
//--------------------------------------------------------------------------------------------------------------------------------------------------

// Define the interaction rate for both targets based on theory for a particular detector
//---------------------------------------------------------------------------------------
const double c0_conversion = (c*1.E-8)*pow(10.,-28.)*86400*(197.327)*(197.327);
double N_T(double A) {return (6.02E26)/A;}
double ev_to_kg = (1.783E-36);
double amu_to_kg = 1.66E-27;                                // convert from amu to kg
double DAMAexposure = 1;//411137.;                              // total exposure of the DAMA experiment (kg x days)
// max_energy = 16.;                                     // maximum energy value we want to test to
// min_energy = 1.;
int npoints = 15;//(max_energy-1)*2;                             // number of DAMA/LIBRA points/bins in region of interest
int Nyear = 2;                                              // number of years of data taking
double SABREexposure = 50*Nyear*365.25;                     // total exposure of SABRE experiment, assuming 50kg of active target
double infty = TMath::Infinity();
double zero = 0.;


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

TF2 *NaAv = new TF2("NaAv", "NadRdE0([3], x, [0], [1], [2])", 0., infty);
TF2 *NaMod = new TF2("NaMod", "NadRdE([3], x, [0], [1], [2])", 0., infty);
TF2 *IAv = new TF2("IAv", "IdRdE0([3], x, [0], [1], [2])", 0., infty);
TF2 *IMod = new TF2("IMod", "IdRdE([3], x, [0], [1], [2])", 0., infty);

double av_rate(double E, double m, double sig, double d){
    NaAv->SetParameter(0, m*1E9);
    NaAv->SetParameter(1, sig);
    NaAv->SetParameter(2, d);
    NaAv->SetParameter(3, E);
    IAv->SetParameter(0, m*1E9);
    IAv->SetParameter(1, sig);
    IAv->SetParameter(2, d);
    IAv->SetParameter(3, E);
    return (ANa/(ANa+AI))*NaAv->TF2::Integral(zero, infty, 1E-7) + (AI/(ANa+AI))*IAv->TF2::Integral(zero, infty, 1E-7) ;
}

double mod_rate(double E, double m, double sig, double d){
    NaMod->SetParameter(0, m*1E9);
    NaMod->SetParameter(1, sig);
    NaMod->SetParameter(2, d);
    NaMod->SetParameter(3, E);
    IMod->SetParameter(0, m*1E9);
    IMod->SetParameter(1, sig);
    IMod->SetParameter(2, d);
    IMod->SetParameter(3, E);
    return (ANa/(ANa+AI))*NaMod->TF2::Integral(zero, infty, 1E-7) + (AI/(ANa+AI))*IMod->TF2::Integral(zero, infty, 1E-7) ;
}

// Target specific rate for plotting
double Na_av_rate(double E, double m, double sig, double d){
    NaAv->SetParameter(0, m*1E9);
    NaAv->SetParameter(1, sig);
    NaAv->SetParameter(2, d);
    NaAv->SetParameter(3, E);
    return (ANa/(ANa+AI))*NaAv->TF2::Integral(zero, infty, 1E-7);
}

double Na_mod_rate(double E, double m, double sig, double d){
    NaMod->SetParameter(0, m*1E9);
    NaMod->SetParameter(1, sig);
    NaMod->SetParameter(2, d);
    NaMod->SetParameter(3, E);
    return (ANa/(ANa+AI))*NaMod->TF2::Integral(zero,infty, 1E-7) ;
}

double I_av_rate(double E, double m, double sig, double d){
    IAv->SetParameter(0, m*1E9);
    IAv->SetParameter(1, sig);
    IAv->SetParameter(2, d);
    IAv->SetParameter(3, E);
    return (AI/(ANa+AI))*IAv->TF2::Integral(zero, infty, 1E-7) ;
}

double I_mod_rate(double E, double m, double sig, double d){
    IMod->SetParameter(0, m*1E9);
    IMod->SetParameter(1, sig);
    IMod->SetParameter(2, d);
    IMod->SetParameter(3, E);
    return (AI/(ANa+AI))*IMod->TF2::Integral(zero, infty, 1E-7) ;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                                   P e r f o r m   t h e   f i t
//--------------------------------------------------------------------------------------------------------------------------------------------------

void DAMA_Fit(int model_choice)
{

    // S e t u p   m o d e l
    // ---------------------
    // Select a case to investigate, that is, a set of coupling constants from models.h
    model(model_choice);
    setquenching(0);
    loadbackground();
    //exposure = DAMAexposure;                                                    // for fitting, we want the DAMA exposure
    //
    // Declare variables
    RooRealVar E("E","Observed recoil energy",1,16) ;                  // Er passed in units of keV
    RooRealVar mass("mass","Dark matter mass",30,150) ;               // Mass passed in units of GeV
    RooRealVar sigma("sigma","Dark matter cross section",1E-41,5E-38) ;   // Sigma passed in units of cm^2
    //RooRealVar sigma("sigma","Dark matter cross section",1E-30,5E-27) ;   // Sigma passed in units of cm^2
    RooRealVar delta("delta","Dark matter mass splitting",0E3,0E3) ;     // delta passed in units of eV

    // Build average and modulating rate p.d.f in terms of Er, mass, sigma, and delta
    RooAbsReal* AvRate = bindFunction("AvRate", av_rate, E, mass, sigma, delta) ;
    RooAbsReal* ModRate = bindFunction("ModRate", mod_rate, E, mass, sigma, delta) ;
    RooAbsReal* IModRate = bindFunction("IModRate", I_mod_rate, E, mass, sigma, delta) ;
    RooAbsReal* NaModRate = bindFunction("NaModRate", Na_mod_rate, E, mass, sigma, delta) ;
    RooAbsReal* IAvRate = bindFunction("IAvRate", I_av_rate, E, mass, sigma, delta) ;
    RooAbsReal* NaAvRate = bindFunction("NaAvRate", Na_av_rate, E, mass, sigma, delta) ;

    // Construct plot frame in energy
    RooPlot* Eframe1 = E.frame(Title("Modulating Interaction Rate")) ;
    RooPlot* Eframe2 = E.frame(Title("AvRate")) ;



    // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
    // ---------------------------------------------------------------------------
    //
    // mass.setVal(10) ;
    // sigma.setVal(sig_lt) ;
    // delta.setVal(model_delta) ;
    // ModRate->plotOn(Eframe1) ;
    // AvRate->plotOn(Eframe2,LineColor(kRed)) ;

    // I m p o r t   d a t a
    // ---------------------

    // Create a dataset of the DAMA data imported from a .txt file
    TGraph2D *damadata = new TGraph2D("dama_data_err.txt") ;
    RooRealVar y("y","y",0,1) ;
    RooDataSet dxy("dxy","dxy",RooArgSet(E,y),StoreError(RooArgSet(E,y))) ;

    for (int i=0 ; i<=npoints-1 ; i++) {
        double X,Y,Z;
        X = damadata->GetX()[i] ;
        Y = damadata->GetY()[i];
        Z = damadata->GetZ()[i];
        // Set X value and error
        E = X;
        if (i==14) {
          E.setError(4);
        }
        else{
          E.setError(0.25) ;
        }

        // Set Y value and error
        y = Y*Threshold(X);
        y.setError(Z*Threshold(X)) ;

        dxy.add(RooArgSet(E,y)) ;
    }
    dxy.plotOnXY(Eframe1,YVar(y)) ;

    // F i t   m o d e l   t o   d a t a
    // ---------------------------------

    // Fit rate to data
    ModRate->chi2FitTo(dxy,YVar(y)) ;

    // Print values of mass, sigma and delta (that now reflect fitted values and errors)
    mass.Print() ;
    sigma.Print() ;
    delta.Print() ;

    // Calculate and print Chi2 value
    RooAbsReal* Chi2 = ModRate->createChi2(dxy, YVar(y)) ;
    cout<< "chi2 value is " << Chi2->getVal() <<endl;



    // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
    // -------------------------------------------------------

    ModRate->plotOn(Eframe1) ;

    //IAvRate->plotOn(Eframe2,LineColor(kGreen),LineStyle(kDashed)) ;
    //NaAvRate->plotOn(Eframe2,LineColor(kRed),LineStyle(kDashed)) ;

    RooHist* hpull = Eframe1->pullHist() ;

    IModRate->plotOn(Eframe1,LineColor(kGreen),LineStyle(kDashed)) ;
    NaModRate->plotOn(Eframe1,LineColor(kRed),LineStyle(kDashed)) ;
    AvRate->plotOn(Eframe2,LineColor(kBlue),Precision(1E-1)) ;

    RooPlot* Eframe3 = E.frame(Title(" ")) ;
    Eframe3->addPlotable(hpull,"P") ;

    gStyle->SetOptStat(0);
    TCanvas *canv = new TCanvas("canv","c1",800,600);
    TPad *pad1 = new TPad("pad1","pad1",0,0.28,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
    pad1->SetBottomMargin(0.03);
    pad1->SetBorderMode(0);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.4);
    pad2->SetBorderMode(0);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    Eframe1->GetYaxis()->SetTitleOffset(0.8);
    Eframe1->GetYaxis()->SetTitleSize(0.04);
    //Eframe1->GetXaxis()->SetLabelOffset(0.045);
    //Eframe1->GetXaxis()->SetLabelSize(0.05);
    //Eframe1->GetXaxis()->SetTitle("");
    Eframe1->GetYaxis()->SetTitle("cpd/kg/keV_{ee}");
    Eframe1->GetYaxis()->SetRangeUser(0,0.03);
    Eframe1->GetXaxis()->SetLabelOffset(999);
    Eframe1->GetXaxis()->SetLabelSize(0);
    Eframe1->Draw();

    auto* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetHeader(Form("Fits for model %d",model_choice), "C");
    legend->AddEntry((TObject*)0, Form("mass = %.2f GeV",mass.getVal()), "");
    legend->AddEntry((TObject*)0, Form("#sigma = %0.2e cm^{2}",sigma.getVal()), "");
    legend->AddEntry((TObject*)0, Form("#delta = %.2f keV",delta.getVal()*1E-3), "");
    gStyle->SetLegendTextSize(0.05);
    //legend->AddEntry(IModRate, "I interaction rate",LineColor(kGreen),LineStyle(kDashed));
    //legend->AddEntry(NaModRate, "Na interaction rate");
    //legend->AddEntry(ModRate, "Total interaction rate");
    legend->Draw();

    pad2->cd();
    Eframe3->GetYaxis()->SetTitleOffset(0.55);
    Eframe3->GetYaxis()->SetTitleSize(0.13);
    Eframe3->GetYaxis()->SetLabelSize(0.1);
    Eframe3->GetXaxis()->SetLabelOffset(0.06);
    Eframe3->GetXaxis()->SetLabelSize(0.11);
    Eframe3->GetXaxis()->SetTitleOffset(1.2);
    Eframe3->GetXaxis()->SetTitleSize(0.1);
    Eframe3->GetXaxis()->SetTitle("Recoil energy (keV_{ee})");
    //Eframe3->GetYaxis()->SetTitle("Pull");
    //Eframe3->GetXaxis()->SetRangeUser(2.91,3.2);
    Eframe3->Draw();

    TLine *line3 = new TLine(1,0,max_energy,0);
    line3->SetLineColor(4);
    line3->Draw();

}
