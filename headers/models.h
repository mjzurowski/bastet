//
// Header file to house possible DM models fits can be performed to.
// Choice is made by calling model function within the fitting file.
// This will then set the coupling constants to chosen values.
// Mass should be input in units of GeV, cross section in cm^2, and delta in eV.
//

#ifndef models_h
#define models_h

#include "veldist.h"

double jx;                                                  // DM spin - depends on model under investigation
double model_delta;
double m_lt;                        // computed lowest tension mass
double sig_lt;                       // copmuted lowest tension mass
double m_bf;
double sig_bf;
double d_bf;
double m_bf_str;
double sig_bf_str;
double d_bf_str;

void model(int i ){
    if(i==1){ // Choose case 1
        model_delta = 22.84E3;
        m_lt = 11.08;
        sig_lt = 3.93E-27;
        m_bf = 13.87;
        sig_bf = 7.53E-29;//1.26E-26;
        d_bf = 20.17E3;
        m_bf_str = 14.72;
        sig_bf_str = 4.89E-29;
        d_bf_str = 19.81E3;
        jx = 0.;
        c0_1=0.;
        c1_1=0.;
        c0_4=0.;
        c1_4=0.;
        c0_5=0.;
        c1_5=0.;
        c0_6=0.;
        c1_6=0.;
        c0_7=2*0.6827;
        c1_7=2*0.7307;
        c0_10=0.;
        c1_10=0.;
        c_Helm=0.;
    }
    else if(i==2){ // Choose case 2
        model_delta = 23.7E3;
        m_lt = 11.64;
        sig_lt = 4.68E-28;
        m_bf = 13.47;
        sig_bf = 2.09E-29;
        d_bf = 20.82E3;
        m_bf_str = 14.29;
        sig_bf_str = 1.36E-29;
        d_bf_str = 20.67E3;
        jx = 0.5;
        c0_1=0.;
        c1_1=0.;
        c0_4=2*-0.0014;
        c1_4=2*-0.0015;
        c0_5=2*-0.032;
        c1_5=2*-0.0166;
        c0_6=2*0.692;
        c1_6=2*0.7217;
        c0_7=0.;
        c1_7=0.;
        c0_10=0.;
        c1_10=0.;
        c_Helm=0.;
    }
    else if(i==3){ // Choose case 3
        model_delta = 23.4E3;
        m_lt = 11.36;
        sig_lt = 5.71E-32;
        m_bf = 13.17;//12.453;
        sig_bf = 2.45E-33;//9.926E-32;
        d_bf = 20.42E3;
        m_bf_str = 13.96;
        sig_bf_str = 1.26E-33;
        d_bf_str = 19.70E3;
        jx = 1.;
        c0_1=0.;
        c1_1=0.;
        c0_4=2*0.0717;
        c1_4=2*0.0753;
        c0_5=2*0.1892;
        c1_5=2*0.9764;
        c0_6=0.;
        c1_6=0.;
        c0_7=0.;
        c1_7=0.;
        c0_10=0.;
        c1_10=0.;
        c_Helm=0.;
    }
    else if(i==4){ // Choose best elastic, SI fit for D/L phase 2
      cout<< "High mass elastic O1"<<endl;
        model_delta = 0;
        m_lt = 11.14;
        sig_lt = 9.39E-40;
        jx = 0;
        c0_1=1.;//0.34;//2*0.12;//
        c1_1=0.;//1.66;//2*0.88;//
        c0_4=0.;
        c1_4=0.;
        c0_5=0.;
        c1_5=0.;
        c0_6=0.;
        c1_6=0.;
        c0_7=0.;
        c1_7=0.;
        c0_10=0.;
        c1_10=0.;
        c_Helm=1.;
    }
    else if(i==5){ // Choose best elastic, SI fit for D/L phase 2
      cout<< "Low mass elastic O7"<<endl;
        model_delta = 0;
        m_lt = 13.41;//49.24;//
        sig_lt = 4.75E-30;//1.35E-30;//
        jx = 0;
        c0_1=0;
        c1_1=0;
        c0_4=0.;
        c1_4=0.;
        c0_5=0.;
        c1_5=0.;
        c0_6=0.;
        c1_6=0.;
        c0_7=-3.32;//0.35;//
        c1_7=5.32;//1.65;//
        c0_10=0.;
        c1_10=0.;
        c_Helm=1.;
    }
    else if(i==6){ // Choose best elastic, SI fit for D/L phase 2
      cout<< "Low mass elastic O4"<<endl;
        model_delta = 0;
        m_lt = 11.22;//44.71;//
        sig_lt = 2.95E-36;//x5.96E-36;//
        jx = 0.5;
        c0_1=0;
        c1_1=0;
        c0_4=2.71;//-7.34;//
        c1_4=-0.71;//9.34;//
        c0_5=0.;
        c1_5=0.;
        c0_6=0.;
        c1_6=0.;
        c0_7=0;
        c1_7=0;
        c0_10=0.;
        c1_10=0.;
        c_Helm=1.;
    }
    else if(i==7){ // Choose best elastic, SI fit for D/L phase 2
        cout<< "Low mass elastic O5"<<endl;
        model_delta = 0;
        m_lt = 8.34;//96.13;//
        sig_lt = 1.62E-29;//3.63E-34;//
        jx = 0.5;
        c0_1=0;
        c1_1=0;
        c0_4=0;
        c1_4=0;
        c0_5=0.39;//-4.74;//
        c1_5=1.61;//6.74;//
        c0_6=0.;
        c1_6=0.;
        c0_7=0;
        c1_7=0;
        c0_10=0.;
        c1_10=0.;
        c_Helm=1.;
    }
    else if(i==8){ // Choose best elastic, SI fit for D/L phase 2
        cout<< "Low mass elastic O6"<<endl;
        model_delta = 0;
        m_lt = 8.09;
        sig_lt = 5.05E-28;
        jx = 0.5;
        c0_1=0;
        c1_1=0;
        c0_4=0;
        c1_4=0;
        c0_5=0.;
        c1_5=0.;
        c0_6=-6.20;
        c1_6=8.20;
        c0_7=0;
        c1_7=0;
        c0_10=0.;
        c1_10=0.;
        c_Helm=1.;
    }
    else if(i==9){ // Choose best elastic, SI fit for D/L phase 2
      cout<< "pSIDM model"<<endl;
        model_delta = 18.3E3;
        m_lt = 12.1;
        sig_lt = 7.95E-35;
        jx = 0.5;
        c0_1=0;
        c1_1=0;
        c0_4=0.972;
        c1_4=1.028;
        c0_5=0.;
        c1_5=0.;
        c0_6=0.;
        c1_6=0.;
        c0_7=0;
        c1_7=0;
        c0_10=0.;
        c1_10=0.;
        c_Helm=1.;
    }
    else{
        cout<< "Invalid model choice, please try again"<<endl;
    }
}



#endif /* models_h */
