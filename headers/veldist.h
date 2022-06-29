#ifndef veldist_H
#define veldist_H


#include "formfactors.h"
#include "../config.h"
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------- Velocity distribution --------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------------------
// To save on computation time, we follow the notation from arxiv 1607.04418 and split all velocity distributions for DM into two, labelled g and h.
// These correspond to \int v*f(v)dv and \int v^3*f(v)dv respectively, and are ultimately multiplied by the form factors with the same label.

const double c = 299792458.;
double c0_squared(double sig, double Mw){   // used to convert from c0^2 to proton cross sections
  double conversion = 2.5682E45; //should note where this number has come from
  return conversion*TMath::Pi()*sig/(mu(Mw,1.)*mu(Mw,1.));
}

double vmin(double Er,double M,double A,double del){
  // if(Er==0){return pow(2*del/mu(M,A),0.5);}
  // else{
    double speed=(c*1.E-3)*fabs(del +(MT(A)*Er/mu(M,A)))/pow(2.*MT(A)*Er,0.5);
    return speed;
  //}
}



/*
// SHM+Stream
double vesc = 520;
double v0 = 232.8;
double RhoDM = 0.5;
*/

//----------------------------------------------------Cross section----------------------------------------------------
//
// Cross sections for each target are constructed in formfactors.h, and loaded here as FNa and FI
// dsigdE is the general expression for cross section defined in my moleskine/Cirelli paper
// These terms are split up by velocity dependence to be multiplied by the appropriate vel integral.
//


// For newly instrumented EFT FFs
//double dsigdE(double A, double Mw){return MT(A)/(2*pow(mu(Mw,1),2));}      // Sigma normalisation
double dsigdE(double A, double Mw){return MT(A)/(2.*TMath::Pi());}      // c0_squared normalisation

// new g cross section terms
double gFNa(double jx, double Er, double Mw, double minv){return gNa_total_form(jx, Er, Mw, minv);}
double gFI(double jx, double Er, double Mw, double minv){return gI_total_form(jx, Er, Mw, minv);}
double gsigma1Na(double Er,double Mw,double jx, double minv){return gFNa(jx,Er,Mw, minv)*dsigdE(ANa,Mw);}
double gsigma1I(double Er,double Mw,double jx, double minv){return gFI(jx,Er,Mw, minv)*dsigdE(AI,Mw);}

// h cross section terms
double hFNa(double jx, double Er, double Mw){return hNa_total_form(jx, Er, Mw);}
double hFI(double jx, double Er, double Mw){return hI_total_form(jx, Er, Mw);}
double hsigma1Na(double Er,double Mw,double jx){return hFNa(jx,Er,Mw)*dsigdE(ANa,Mw);}
double hsigma1I(double Er,double Mw,double jx){return hFI(jx,Er,Mw)*dsigdE(AI,Mw);}


// To avoid too many long functions later, write a function able to distinguish which target is being used, and so select the correct cross section
// Note that this is only instrumented for Na and I targets, and will otherwise return 0.

double g_crosssection(double Er, double Mw, double jx, double a, double d){
    if(a==ANa){
        //double v = vmin(Er, Mw, a, d)/((1E-3)*c);
        double v = vmin(Er, Mw, a, d)/((1.E-3)*c);
        return gsigma1Na(Er, Mw, jx, v);
    }
    else if (a==AI){
        //double v = vmin(Er, Mw, a, d)/((1E-3)*c);
        double v = vmin(Er, Mw, a, d)/((1.E-3)*c);
        return gsigma1I(Er, Mw, jx,v);
    }
    else {return 0.;}

}

double h_crosssection(double Er, double Mw, double jx, double a){
    if(a==ANa){
        return hsigma1Na(Er, Mw, jx);
    }
    else if (a==AI){
        return hsigma1I(Er, Mw, jx);
    }
    else {return 0.;}
}

//----------------------------------------------------Velocity distribution----------------------------------------------------------
// graphs are loaded in config.h

//  - - - - - - - - - - - - - - - -  Modulation Component - - - - - - - - - - - - - - - -
// Define the h velocity distribution
double h_fm(double Er, double Mw, double d, double a){
    double min_vel = vmin(Er, Mw, a, d);
    if(min_vel<=(vesc+vE)){
        double h = modh->Eval(min_vel,0,"S");
        return h/((1.E-3)*c);
    }
    else{
        return 0.0;
    }
}

// Define the g velocity distribution
double g_fm(double Er, double Mw, double d, double a){
    double min_vel = vmin(Er, Mw, a, d);
    if(min_vel<=(vesc+vE)){
        double g = modg->Eval(min_vel,0,"S");
        return g*(c*1.E-3);
    }
    else{
        return 0.0;
    }
}

// Add the two components together and multiply by the cross section
double fm(double Er, double Mw, double d, double a, double jx){
    return h_fm(Er,Mw,d,a)*h_crosssection(Er,Mw,jx,a)+g_fm(Er,Mw,d,a)*g_crosssection(Er,Mw,jx,a,d);
}

//  - - - - - - - - - - - - - - - - Average Component - - - - - - - - - - - - - - - -
// Define the h velocity distribution
double h_f0(double Er, double Mw, double d, double a){
    double min_vel = vmin(Er, Mw, a, d);
    if(min_vel<=(vesc+vE)){
        double h = avh->Eval(min_vel,0,"S");
        return h/((1.E-3)*c);
    }
    else{
        return 0.0;
    }
}

// Define the g velocity distribution
double g_f0(double Er, double Mw, double d, double a){
    double min_vel = vmin(Er, Mw, a, d);
    if(min_vel<=(vesc+vE)){
        double g = avg->Eval(min_vel,0,"S");
        return g*(c*1.E-3);
    }
    else{
        return 0.0;
    }
}

// Add the two components together and multiply by the cross section
double f0(double Er, double Mw, double d, double a, double jx){
    return h_f0(Er,Mw,d,a)*h_crosssection(Er,Mw,jx,a)+g_f0(Er,Mw,d,a)*g_crosssection(Er,Mw,jx,a,d);
}


#endif
