#ifndef formfactors_H
#define formfactors_H

//--------------------------------------------------------Form Factors--------------------------------------------------------
// A setup for the various form factors used in DM scattering. Based on expressions from 1308.6288 version
// Note that following arxiv 1607.04418 these can be separated into velocity independent and velocity dependent
// form factors. Following the notation developed in that paper, independent FFs will be labelled g and dependent h.
//
// Atomic mass numbers for targets of interest
const double ANa = 22.99;
const double AI = 126.90;
const double mp =  0.9314941*1E9; //mass of a nulceon in eV.

// Terms in common/expressions used across all operators
double Q(double Er,double A){return pow(2.*A*mp * Er , 0.5);}
double B(double A){return pow(41.467/(45.*pow(A,-1./3.)-25.*pow(A,-2./3.)),0.5);}
double Y(double Er, double A){return pow((1/(197.327 *1E6))*Q(Er,A)*B(A)/2.,2.);}
double C(double jx){return (4.*jx/3.)*(jx+1.);}
double MT(double A){return mp*A;}               // Mass of target in eV (based on A value)
double mu(double Mw,double A){return ((MT(A)*Mw)/(MT(A)+Mw));} // Reduced mass of the system in units of eV
const double Q_conv = 1/(197.327 *1E6); // converts q into units of inverse fm from eV.


//--------------------------------------------------------Coupling coefficients--------------------------------------------------------
// Can input here manually, or load a header (either here, or it might be better to include in Fit_to_DAMA?) with the desired values.
double c0_1;//
double c1_1;//
double c0_4;//
double c1_4;//
double c0_5;//
double c1_5;//
double c0_6;//
double c1_6;//
double c0_7;//
double c1_7;//
double c0_10;//
double c1_10;//
double c_Helm;//

//--------------------------------------------------------Target dependent terms--------------------------------------------------------

// For the np basis:

// 23Na
double NaFS1pp(double y){return exp(-2.*y)*(0.273-0.824*y+1.19*pow(y,2.)-0.477*pow(y,3.)+0.0593*pow(y,4.));}
double NaFS1np(double y){return exp(-2.*y)*(0.0219-0.0578*y+0.036*pow(y,2.)-0.003*pow(y,3.)+0.000363*pow(y,4.));}
double NaFS1nn(double y){return exp(-2.*y)*(0.00176-0.00396*y+0.00228*pow(y,2.)-0.0000195*pow(y,3.));}
double NaFS2pp(double y){return exp(-2.*y)*(0.136-0.267*y+0.458*pow(y,2.)-0.112*pow(y,3.)+0.00828*pow(y,4.));}
double NaFS2np(double y){return exp(-2.*y)*(0.011-0.03*y+0.0217*pow(y,2.)-0.00897*pow(y,3.)+0.000592*pow(y,4.));}
double NaFS2nn(double y){return exp(-2.*y)*(0.000882-0.0031*y+0.00399*pow(y,2.)-0.00203*pow(y,3.)+0.000409*pow(y,4.));}
double NaFMpp(double y){return exp(-2.*y)*(120.-180.*y+87.*pow(y,2.)-17.*pow(y,3.)+1.2*pow(y,4.));}
double NaFMnp(double y){return exp(-2.*y)*(130.-200.*y+100.*pow(y,2.)-20.*pow(y,3.)+1.5*pow(y,4.));}
double NaFMnn(double y){return exp(-2.*y)*(140.-220.*y+120.*pow(y,2.)-25.*pow(y,3.)+1.8*pow(y,4.));}
double NaFDpp(double y){return exp(-2.*y)*(0.231-0.185*y+0.0502*pow(y,2.));}
double NaFDnp(double y){return exp(-2.*y)*(0.0812-0.065*y+0.0138*pow(y,2.));}
double NaFDnn(double y){return exp(-2.*y)*(0.0286-0.0228*y+0.00462*pow(y,2.));}
double NaFS1Dpp(double y){return exp(-2.*y)*(-0.25+0.48*y-0.29*pow(y,2.)+0.049*pow(y,3.));}
double NaFS1Dpn(double y){return exp(-2.*y)*(-0.088+0.17*y-0.081*pow(y,2.)+0.011*pow(y,3.));}
double NaFS1Dnp(double y){return exp(-2.*y)*(-0.02+0.031*y-0.0076*pow(y,2.)+0.00027*pow(y,3.));}
double NaFS1Dnn(double y){return exp(-2.*y)*(-0.0071+0.011*y-0.003*pow(y,2.));}

// 127I
double IFS1pp(double y){return exp(-2.*y)*(0.26-1.6*y+5.3*pow(y,2.)-8.9*pow(y,3.)+8.7*pow(y,4.)-4.9*pow(y,5.)+1.5*pow(y,6.)-0.25*pow(y,7.)+0.016*pow(y,8.));}
double IFS1np(double y){return exp(-2.*y)*(0.065-0.46*y+1.3*pow(y,2.)-1.8*pow(y,3.)+1.4*pow(y,4.)-0.65*pow(y,5.)+0.17*pow(y,6.)-0.026*pow(y,7.)+0.002*pow(y,8.));}
double IFS1nn(double y){return exp(-2.*y)*(0.016-0.13*y+0.37*pow(y,2.)-0.48*pow(y,3.)+0.34*pow(y,4.)-0.14*pow(y,5.)+0.033*pow(y,6.)-0.0048*pow(y,7.)+0.00041*pow(y,8.));}
double IFS2pp(double y){return exp(-2.*y)*(0.13-0.49*y+1.8*pow(y,2.)-2.8*pow(y,3.)+2.7*pow(y,4.)-1.6*pow(y,5.)+0.53*pow(y,6.)-0.092*pow(y,7.)+0.0067*pow(y,8.));}
double IFS2np(double y){return exp(-2.*y)*(0.032-0.13*y+0.26*pow(y,2.)-0.3*pow(y,3.)+0.21*pow(y,4.)-0.098*pow(y,5.)+0.027*pow(y,6.)-0.0042*pow(y,7.)+0.00033*pow(y,8.));}
double IFS2nn(double y){return exp(-2.*y)*(0.008-0.032*y+0.053*pow(y,2.)-0.046*pow(y,3.)+0.025*pow(y,4.)-0.0086*pow(y,5.)+0.0019*pow(y,6.)-0.00026*pow(y,7.));}
double IFMpp(double y){return exp(-2.*y)*(2800.-10000.*y+14000.*pow(y,2.)-9800.*pow(y,3.)+3800.*pow(y,4.)-840.*pow(y,5.)+100.*pow(y,6.)-6.3*pow(y,7.)+0.15*pow(y,8.));}
double IFMnp(double y){return exp(-2.*y)*(3900.-15000.*y+23000.*pow(y,2.)-18000.*pow(y,3.)+7900.*pow(y,4.)-2000.*pow(y,5.)+290.*pow(y,6.)-23.*pow(y,7.)+0.75*pow(y,8.)-0.0048*pow(y,9.));}
double IFMnn(double y){return exp(-2.*y)*(5500.-23000.*y+38000.*pow(y,2.)-32000.*pow(y,3.)+16000.*pow(y,4.)-4600.*pow(y,5.)+790.*pow(y,6.)-75.*pow(y,7.)+3.3*pow(y,8.)-0.041*pow(y,9.)+0.00015*pow(y,10.));}
double IFDpp(double y){return exp(-2.*y)*(0.54-1.3*y+1.6*pow(y,2.)-1.2*pow(y,3.)+0.51*pow(y,4.)-0.11*pow(y,5.)+0.0097*pow(y,6.));}
double IFDnp(double y){return exp(-2.*y)*(0.23-0.65*y+0.79*pow(y,2.)-0.54*pow(y,3.)+0.2*pow(y,4.)-0.04*pow(y,5.)+0.0039*pow(y,6.)-0.00014*pow(y,7.));}
double IFDnn(double y){return exp(-2.*y)*(0.1-0.32*y+0.4*pow(y,2.)-0.25*pow(y,3.)+0.084*pow(y,4.)-0.016*pow(y,5.)+0.0018*pow(y,6.)-0.00011*pow(y,7.));}
double IFS1Dpp(double y){return exp(-2.*y)*(-0.37+1.6*y-3.1*pow(y,2.)+3.4*pow(y,3.)-2.1*pow(y,4.)+0.73*pow(y,5.)-0.13*pow(y,6.)+0.0086*pow(y,7.));}
double IFS1Dpn(double y){return exp(-2.*y)*(-0.16+0.75*y-1.5*pow(y,2.)+1.5*pow(y,3.)-0.82*pow(y,4.)+0.26*pow(y,5.)-0.047*pow(y,6.)+0.0043*pow(y,7.)-0.00015*pow(y,8.));}
double IFS1Dnp(double y){return exp(-2.*y)*(-0.093+0.48*y-0.85*pow(y,2.)+0.79*pow(y,3.)-0.43*pow(y,4.)+0.13*pow(y,5.)-0.021*pow(y,6.)+0.0017*pow(y,7.));}
double IFS1Dnn(double y){return exp(-2.*y)*(-0.04+0.22*y-0.43*pow(y,2.)+0.38*pow(y,3.)-0.19*pow(y,4.)+0.053*pow(y,5.)-0.009*pow(y,6.)+0.00088*pow(y,7.));}


// For the isospin basis:
//23Na
double NaF00_S1(double y){return exp(-2.*y)*(0.0795907-0.235885*y+0.314987*y*y-0.120719*pow(y,3.)+0.0146522*pow(y,4.));}
double NaF11_S1(double y){return exp(-2.*y)*(0.0576531-0.178131*y+0.278909*y*y-0.117715*pow(y,3.)+0.0150154*pow(y,4.));}
double NaF01_S1(double y){return exp(-2.*y)*(0.0677396-0.205029*y+0.295762*y*y-0.119227*pow(y,3.)+0.0148326*pow(y,4.));}
double NaF10_S1(double y){return exp(-2.*y)*(0.0677396-0.205029*y+0.295762*y*y-0.119227*pow(y,3.)+0.0148326*pow(y,4.));}

double NaF00_S2(double y){return exp(-2.*y)*(0.0397953-0.0824771*y+0.126256*y*y-0.0330308*pow(y,3.)+0.00246945*pow(y,4.));}
double NaF11_S2(double y){return exp(-2.*y)*(0.0288265-0.0524812*y+0.104537*y*y-0.0240558*pow(y,3.)+0.00187765*pow(y,4.));}
double NaF01_S2(double y){return exp(-2.*y)*(0.0338698-0.0659295*y+0.113402*y*y-0.027527*pow(y,3.)+0.00196889*pow(y,4.));}
double NaF10_S2(double y){return exp(-2.*y)*(0.0338698-0.0659295*y+0.113402*y*y-0.027527*pow(y,3.)+0.00196889*pow(y,4.));}

double NaF00_S1D(double y){return exp(-2.*y)*(-0.0916196+0.172416*y-0.095927*y*y+0.0149661*pow(y,3.));}
double NaF11_S1D(double y){return exp(-2.*y)*(-0.0374013+0.0727401*y-0.051533*y*y+0.00974633*pow(y,3.));}
double NaF01_S1D(double y){return exp(-2.*y)*(-0.0439447+0.0826981*y-0.053835*y*y+0.0096358*pow(y,3.));}
double NaF10_S1D(double y){return exp(-2.*y)*(-0.0779774+0.151655*y-0.0906034*y*y+0.0151226*pow(y,3.));}

double NaF00_D(double y){return exp(-2.*y)*(0.105467-0.0843733*y+0.020637*y*y);}
double NaF11_D(double y){return exp(-2.*y)*(0.0242634-0.0194107*y+0.00679182*y*y);}
double NaF01_D(double y){return exp(-2.*y)*(0.0505863-0.040469*y+0.0114025*y*y);}
double NaF10_D(double y){return exp(-2.*y)*(0.0505863-0.040469*y+0.0114025*y*y);}

double NaF00_M(double y){return exp(-2.*y)*(132.25-199.333*y+102.389*y*y-20.6679*pow(y,3.)+1.51791*pow(y,4.));}
double NaF11_M(double y){return exp(-2.*y)*(0.25-0.666668*y+0.574725*y*y-0.170869*pow(y,3.)+0.0164309*pow(y,4.));}
double NaF01_M(double y){return exp(-2.*y)*(-5.75+12.0*y-7.86797*y*y+1.87811*pow(y,3.)-0.142785*pow(y,4.));}
double NaF10_M(double y){return exp(-2.*y)*(-5.75+12.0*y-7.86797*y*y+1.87811*pow(y,3.)-0.142785*pow(y,4.));}


//127I
double IF00_S1(double y){return exp(-2.*y)*(0.101358-0.658425*y+2.05704*y*y-3.24061*pow(y,3.)+2.96696*pow(y,4.)-1.58699*pow(y,5.)+0.479893*pow(y,6.)-0.0757009*pow(y,7.)+0.00505039*pow(y,8.)-0.0000350195*pow(y,9.)+9.19758*pow(10.,-8.)*pow(y,10.));}
double IF11_S1(double y){return exp(-2.*y)*(0.0367843-0.202443*y+0.77617*y*y-1.44636*pow(y,3.)+1.55152*pow(y,4.)-0.935991*pow(y,5.)+0.305841*pow(y,6.)-0.0495074*pow(y,7.)+0.00303135*pow(y,8.)-0.0000254083*pow(y,9.)+9.19753*pow(10.,-8.)*pow(y,10.));}
double IF01_S1(double y){return exp(-2.*y)*(0.0610604-0.366349*y+1.23145*y*y-2.10201*pow(y,3.)+2.09138*pow(y,4.)-1.19396*pow(y,5.)+0.376508*pow(y,6.)-0.0602037*pow(y,7.)+0.00383369*pow(y,8.)-4.80553*pow(10.,-6.)*pow(y,9.)+9.19755*pow(10.,-8.)*pow(y,10.));}
double IF10_S1(double y){return exp(-2.*y)*(0.0610604-0.366349*y+1.23145*y*y-2.10201*pow(y,3.)+2.09138*pow(y,4.)-1.19396*pow(y,5.)+0.376508*pow(y,6.)-0.0602037*pow(y,7.)+0.00383369*pow(y,8.)-4.80553*pow(10.,-6.)*pow(y,9.)+9.19755*pow(10.,-8.)*pow(y,10.));}

double IF00_S2(double y){return exp(-2.*y)*(0.0506789-0.192735*y+0.604006*y*y-0.85886*pow(y,3.)+0.789315*pow(y,4.)-0.446953*pow(y,5.)+0.146138*pow(y,6.)-0.0251915*pow(y,7.)+0.00183763*pow(y,8.)-5.09683*pow(10.,-6.)*pow(y,9.)+4.21877*pow(10.,-9.)*pow(y,10.));}
double IF11_S2(double y){return exp(-2.*y)*(0.0183921-0.0675558*y+0.339556*y*y-0.562028*pow(y,3.)+0.57516*pow(y,4.)-0.348172*pow(y,5.)+0.119013*pow(y,6.)-0.0209647*pow(y,7.)+0.00150736*pow(y,8.)-4.62869*pow(10.,-6.)*pow(y,9.)+4.21876*pow(10.,-9.)*pow(y,10.));}
double IF01_S2(double y){return exp(-2.*y)*(0.0305302-0.114124*y+0.445065*y*y-0.687362*pow(y,3.)+0.669869*pow(y,4.)-0.393542*pow(y,5.)+0.131624*pow(y,6.)-0.0229471*pow(y,7.)+0.00166181*pow(y,8.)-2.34066*pow(10.,-7.)*pow(y,9.)+4.21876*pow(10.,-9.)*pow(y,10.));}
double IF10_S2(double y){return exp(-2.*y)*(0.0305302-0.114124*y+0.445065*y*y-0.687362*pow(y,3.)+0.669869*pow(y,4.)-0.393542*pow(y,5.)+0.131624*pow(y,6.)-0.0229471*pow(y,7.)+0.00166181*pow(y,8.)-2.34066*pow(10.,-7.)*pow(y,9.)+4.21876*pow(10.,-9.)*pow(y,10.));}

double IF00_S1D(double y){return exp(-2.*y)*(-0.166894+0.762496*y-1.46788*y*y+1.50582*pow(y,3.)-0.884208*pow(y,4.)+0.293623*pow(y,5.)-0.0512228*pow(y,6.)+0.00385863*pow(y,7.)-0.0000619248*pow(y,8.)+2.40277*pow(10.,-7.)*pow(y,9.));}
double IF11_S1D(double y){return exp(-2.*y)*(-0.0398501+0.14534*y-0.310762*y*y+0.376832*pow(y,3.)-0.260324*pow(y,4.)+0.0974939*pow(y,5.)-0.0171823*pow(y,6.)+0.000871078*pow(y,7.)+0.0000391843*pow(y,8.)+2.40276*pow(10.,-7.)*pow(y,9.));}
double IF01_S1D(double y){return exp(-2.*y)*(-0.0661495+0.274086*y-0.522867*y*y+0.583043*pow(y,3.)-0.379992*pow(y,4.)+0.135416*pow(y,5.)-0.0232947*pow(y,6.)+0.00129089*pow(y,7.)+0.0000358032*pow(y,8.)-2.40276*pow(10.,-7.)*pow(y,9.));}
double IF10_S1D(double y){return exp(-2.*y)*(-0.100541+0.409451*y-0.830535*y*y+0.917057*pow(y,3.)-0.576942*pow(y,4.)+0.202249*pow(y,5.)-0.0360983*pow(y,6.)+0.00255991*pow(y,7.)-0.0000130629*pow(y,8.)-2.40276*pow(10.,-7.)*pow(y,9.));}

double IF00_D(double y){return exp(-2.*y)*(0.274804-0.725884*y+0.906019*y*y-0.627482*pow(y,3.)+0.249047*pow(y,4.)-0.0520121*pow(y,5.)+0.00482199*pow(y,6.)-0.0000966194*pow(y,7.)+6.31504*pow(10.,-7.)*pow(y,8.));}
double IF11_D(double y){return exp(-2.*y)*(0.0431713-0.077312*y+0.115019*y*y-0.0918722*pow(y,3.)+0.0483834*pow(y,4.)-0.0119796*pow(y,5.)+0.000939474*pow(y,6.)-0.0000430693*pow(y,7.)+6.315*pow(10.,-7.)*pow(y,8.));}
double IF01_D(double y){return exp(-2.*y)*(0.108921-0.241383*y+0.31163*y*y-0.236906*pow(y,3.)+0.106749*pow(y,4.)-0.0237535*pow(y,5.)+0.00196234*pow(y,6.)-0.0000267748*pow(y,7.)+6.31502*pow(10.,-7.)*pow(y,8.));}
double IF10_D(double y){return exp(-2.*y)*(0.108921-0.241383*y+0.31163*y*y-0.236906*pow(y,3.)+0.106749*pow(y,4.)-0.0237535*pow(y,5.)+0.00196234*pow(y,6.)-0.0000267748*pow(y,7.)+6.31502*pow(10.,-7.)*pow(y,8.));}

double IF00_M(double y){return exp(-2.*y)*(4032.23-15747.9*y+24396.4*y*y-19508.4*pow(y,3.)+8873.75*pow(y,4.)-2374.58*pow(y,5.)+370.268*pow(y,6.)-31.609*pow(y,7.)+1.24595*pow(y,8.)-0.0127146*pow(y,9.)+0.0000381874*pow(y,10.));}
double IF11_M(double y){return exp(-2.*y)*(110.247-615.985*y+1356.88*y*y-1529.79*pow(y,3.)+971.522*pow(y,4.)-359.49*pow(y,5.)+77.0056*pow(y,6.)-9.02569*pow(y,7.)+0.495715*pow(y,8.)-0.00789838*pow(y,9.)+0.0000381874*pow(y,10.));}
double IF01_M(double y){return exp(-2.*y)*(-666.74+3164.62*y-5884.09*y*y+5603.13*pow(y,3.)-3011.45*pow(y,4.)+945.876*pow(y,5.)-172.298*pow(y,6.)+17.1571*pow(y,7.)-0.794693*pow(y,8.)+0.0103065*pow(y,9.)-0.0000381874*pow(y,10.));}
double IF10_M(double y){return exp(-2.*y)*(-666.74+3164.62*y-5884.09*y*y+5603.13*pow(y,3.)-3011.45*pow(y,4.)+945.876*pow(y,5.)-172.298*pow(y,6.)+17.1571*pow(y,7.)-0.794693*pow(y,8.)+0.0103065*pow(y,9.)-0.0000381874*pow(y,10.));}


//--------------------------------------------------------Operator O1--------------------------------------------------------
double FNN_11(double FM){return FM;}

double NaFnn_11(double Er){return FNN_11(NaFMnn(Y(Er,ANa)));}
double NaFpp_11(double Er){return FNN_11(NaFMpp(Y(Er,ANa)));}
double NaFnp_11(double Er){return FNN_11(NaFMnp(Y(Er,ANa)));}

double IFnn_11(double Er){return FNN_11(IFMnn(Y(Er,AI)));}
double IFpp_11(double Er){return FNN_11(IFMpp(Y(Er,AI)));}
double IFnp_11(double Er){return FNN_11(IFMnp(Y(Er,AI)));}

// 23Na
// double NaF00_11(double Er){return 0.25*(NaFnn_11(Er)+NaFpp_11(Er)+2*NaFnp_11(Er));}
// double NaF11_11(double Er){return 0.25*(NaFnn_11(Er)+NaFpp_11(Er)-2*NaFnp_11(Er));}
// double NaF01_11(double Er){return 0.25*(-NaFnn_11(Er)+NaFpp_11(Er));}
// double NaF10_11(double Er){return 0.25*(-NaFnn_11(Er)+NaFpp_11(Er));}

// from mathematica
double NaF00_11(double Er){return FNN_11(NaF00_M(Y(Er,ANa)));}
double NaF11_11(double Er){return FNN_11(NaF11_M(Y(Er,ANa)));}
double NaF01_11(double Er){return FNN_11(NaF01_M(Y(Er,ANa)));}
double NaF10_11(double Er){return FNN_11(NaF10_M(Y(Er,ANa)));}

// 23Na total F11
double NaF11(double Er){
    double term1 = c0_1*c0_1*NaF00_11(Er);
    double term2 = c1_1*c1_1*NaF11_11(Er);
    double term3 = c0_1*c1_1*NaF01_11(Er);
    double term4 = c1_1*c0_1*NaF10_11(Er);
    return term1+term2+term3+term4;
}

// 127I
// double IF00_11(double Er){return 0.25*(IFnn_11(Er)+IFpp_11(Er)+2*IFnp_11(Er));}
// double IF11_11(double Er){return 0.25*(IFnn_11(Er)+IFpp_11(Er)-2*IFnp_11(Er));}
// double IF01_11(double Er){return 0.25*(-IFnn_11(Er)+IFpp_11(Er));}
// double IF10_11(double Er){return 0.25*(-IFnn_11(Er)+IFpp_11(Er));}

// from mathematica
double IF00_11(double Er){return FNN_11(IF00_M(Y(Er,AI)));}
double IF11_11(double Er){return FNN_11(IF11_M(Y(Er,AI)));}
double IF01_11(double Er){return FNN_11(IF01_M(Y(Er,AI)));}
double IF10_11(double Er){return FNN_11(IF10_M(Y(Er,AI)));}

// 127I total F11
double IF11(double Er){
    double term1 = c0_1*c0_1*IF00_11(Er);
    double term2 = c1_1*c1_1*IF11_11(Er);
    double term3 = c0_1*c1_1*IF01_11(Er);
    double term4 = c1_1*c0_1*IF10_11(Er);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O4--------------------------------------------------------
double FNN_44(double jx, double FS1, double FS2, double Er, double A){return C(jx)*(FS1+FS2)/16.;}

double NaFnn_44(double jx, double Er){return FNN_44(jx,NaFS1nn(Y(Er,ANa)),NaFS2nn(Y(Er,ANa)),Er,ANa);}
double NaFpp_44(double jx, double Er){return FNN_44(jx,NaFS1pp(Y(Er,ANa)),NaFS2pp(Y(Er,ANa)),Er,ANa);}
double NaFnp_44(double jx, double Er){return FNN_44(jx,NaFS1np(Y(Er,ANa)),NaFS2np(Y(Er,ANa)),Er,ANa);}

double IFnn_44(double jx, double Er){return FNN_44(jx,IFS1nn(Y(Er,AI)),IFS2nn(Y(Er,AI)),Er,AI);}
double IFpp_44(double jx, double Er){return FNN_44(jx,IFS1pp(Y(Er,AI)),IFS2pp(Y(Er,AI)),Er,AI);}
double IFnp_44(double jx, double Er){return FNN_44(jx,IFS1np(Y(Er,AI)),IFS2np(Y(Er,AI)),Er,AI);}

// mathematica:
// 23Na
double NaF00_44(double jx, double Er){return FNN_44(jx,NaF00_S1(Y(Er,ANa)),NaF00_S2(Y(Er,ANa)),Er,ANa);}
double NaF11_44(double jx, double Er){return FNN_44(jx,NaF11_S1(Y(Er,ANa)),NaF11_S2(Y(Er,ANa)),Er,ANa);}
double NaF01_44(double jx, double Er){return FNN_44(jx,NaF01_S1(Y(Er,ANa)),NaF01_S2(Y(Er,ANa)),Er,ANa);}
double NaF10_44(double jx, double Er){return FNN_44(jx,NaF10_S1(Y(Er,ANa)),NaF10_S2(Y(Er,ANa)),Er,ANa);}

// // 23Na
// double NaF00_44(double jx, double Er){return 0.25*(NaFnn_44(jx, Er)+NaFpp_44(jx,Er)+2*NaFnp_44(jx, Er));}
// double NaF11_44(double jx, double Er){return 0.25*(NaFnn_44(jx,Er)+NaFpp_44(jx,Er)-2*NaFnp_44(jx,Er));}
// double NaF01_44(double jx, double Er){return 0.25*(-NaFnn_44(jx,Er)+NaFpp_44(jx,Er));}
// double NaF10_44(double jx, double Er){return 0.25*(-NaFnn_44(jx,Er)+NaFpp_44(jx,Er));}

// 23Na total F44
double NaF44(double jx, double Er){
    double term1 = c0_4*c0_4*NaF00_44(jx,Er);
    double term2 = c1_4*c1_4*NaF11_44(jx,Er);
    double term3 = c0_4*c1_4*NaF01_44(jx,Er);
    double term4 = c1_4*c0_4*NaF10_44(jx,Er);
    return term1+term2+term3+term4;
}

// mathematica:
// 127I
double IF00_44(double jx, double Er){return FNN_44(jx,IF00_S1(Y(Er,AI)),IF00_S2(Y(Er,AI)),Er,AI);}
double IF11_44(double jx, double Er){return FNN_44(jx,IF11_S1(Y(Er,AI)),IF11_S2(Y(Er,AI)),Er,AI);}
double IF01_44(double jx, double Er){return FNN_44(jx,IF01_S1(Y(Er,AI)),IF01_S2(Y(Er,AI)),Er,AI);}
double IF10_44(double jx, double Er){return FNN_44(jx,IF10_S1(Y(Er,AI)),IF10_S2(Y(Er,AI)),Er,AI);}

// // 127I
// double IF00_44(double jx, double Er){return 0.25*(IFnn_44(jx,Er)+IFpp_44(jx,Er)+2*IFnp_44(jx,Er));}
// double IF11_44(double jx, double Er){return 0.25*(IFnn_44(jx,Er)+IFpp_44(jx,Er)-2*IFnp_44(jx,Er));}
// double IF01_44(double jx, double Er){return 0.25*(-IFnn_44(jx,Er)+IFpp_44(jx,Er));}
// double IF10_44(double jx, double Er){return 0.25*(-IFnn_44(jx,Er)+IFpp_44(jx,Er));}

// 127I total F44
double IF44(double jx, double Er){
    double term1 = c0_4*c0_4*IF00_44(jx,Er);
    double term2 = c1_4*c1_4*IF11_44(jx,Er);
    double term3 = c0_4*c1_4*IF01_44(jx,Er);
    double term4 = c1_4*c0_4*IF10_44(jx,Er);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O45--------------------------------------------------------
double FNN_45(double jx, double FS1D, double Er, double A){return C(jx)*pow(Q(Er,A)/mp,2.)*FS1D/(4.);}                // 1308.6288 version
//12/02/2020 changed 4 in denom to 8
//double FNN_45(double jx, double FS1D, double Er, double A){return C(jx)*pow(Q(Er,A),2)*FS1D/(4*mp);}

double NaFnn_45(double jx, double Er){return FNN_45(jx,NaFS1Dnn(Y(Er,ANa)),Er,ANa);}
double NaFpp_45(double jx, double Er){return FNN_45(jx,NaFS1Dpp(Y(Er,ANa)),Er,ANa);}
double NaFpn_45(double jx, double Er){return FNN_45(jx,NaFS1Dpn(Y(Er,ANa)),Er,ANa);}
double NaFnp_45(double jx, double Er){return FNN_45(jx,NaFS1Dnp(Y(Er,ANa)),Er,ANa);}

double IFnn_45(double jx, double Er){return FNN_45(jx,IFS1Dnn(Y(Er,AI)),Er,AI);}
double IFpp_45(double jx, double Er){return FNN_45(jx,IFS1Dpp(Y(Er,AI)),Er,AI);}
double IFpn_45(double jx, double Er){return FNN_45(jx,IFS1Dpn(Y(Er,AI)),Er,AI);}
double IFnp_45(double jx, double Er){return FNN_45(jx,IFS1Dnp(Y(Er,AI)),Er,AI);}

// mathematica:
// 23Na
double NaF00_45(double jx, double Er){return FNN_45(jx,NaF00_S1D(Y(Er,ANa)),Er,ANa);}
double NaF11_45(double jx, double Er){return FNN_45(jx,NaF11_S1D(Y(Er,ANa)),Er,ANa);}
double NaF01_45(double jx, double Er){return FNN_45(jx,NaF01_S1D(Y(Er,ANa)),Er,ANa);}
double NaF10_45(double jx, double Er){return FNN_45(jx,NaF10_S1D(Y(Er,ANa)),Er,ANa);}

// // 23Na
// double NaF00_45(double jx, double Er){return 0.25*(NaFnn_45(jx,Er)+NaFpp_45(jx,Er)+NaFpn_45(jx,Er)+NaFnp_45(jx,Er));}
// double NaF11_45(double jx, double Er){return 0.25*(NaFnn_45(jx,Er)+NaFpp_45(jx,Er)-NaFpn_45(jx,Er)-NaFnp_45(jx,Er));}
// double NaF01_45(double jx, double Er){return 0.25*(-NaFnn_45(jx,Er)+NaFpp_45(jx,Er)-NaFpn_45(jx,Er)+NaFnp_45(jx,Er));}
// double NaF10_45(double jx, double Er){return 0.25*(-NaFnn_45(jx,Er)+NaFpp_45(jx,Er)+NaFpn_45(jx,Er)-NaFnp_45(jx,Er));}

// 23Na total F45
double NaF45(double jx, double Er){
    double term1 = c0_4*c0_5*NaF00_45(jx,Er);
    double term2 = c1_4*c1_5*NaF11_45(jx,Er);
    double term3 = c0_4*c1_5*NaF10_45(jx,Er);
    double term4 = c1_4*c0_5*NaF01_45(jx,Er);
    return term1+term2+term3+term4;
}

// mathematica:
// 127I
double IF00_45(double jx, double Er){return FNN_45(jx,IF00_S1D(Y(Er,AI)),Er,AI);}
double IF11_45(double jx, double Er){return FNN_45(jx,IF11_S1D(Y(Er,AI)),Er,AI);}
double IF01_45(double jx, double Er){return FNN_45(jx,IF01_S1D(Y(Er,AI)),Er,AI);}
double IF10_45(double jx, double Er){return FNN_45(jx,IF10_S1D(Y(Er,AI)),Er,AI);}

// // 127I
// double IF00_45(double jx, double Er){return 0.25*(IFnn_45(jx,Er)+IFpp_45(jx,Er)+IFpn_45(jx,Er)+IFnp_45(jx,Er));}
// double IF11_45(double jx, double Er){return 0.25*(IFnn_45(jx,Er)+IFpp_45(jx,Er)-IFpn_45(jx,Er)-IFnp_45(jx,Er));}
// double IF01_45(double jx, double Er){return 0.25*(-IFnn_45(jx,Er)+IFpp_45(jx,Er)-IFpn_45(jx,Er)+IFnp_45(jx,Er));}
// double IF10_45(double jx, double Er){return 0.25*(-IFnn_45(jx,Er)+IFpp_45(jx,Er)+IFpn_45(jx,Er)-IFnp_45(jx,Er));}

// 127I total F45
double IF45(double jx, double Er){
    double term1 = c0_4*c0_5*IF00_45(jx,Er);
    double term2 = c1_4*c1_5*IF11_45(jx,Er);
    double term3 = c0_4*c1_5*IF10_45(jx,Er);
    double term4 = c1_4*c0_5*IF01_45(jx,Er);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O46--------------------------------------------------------
// double FNN_46(double jx, double FS2, double Er, double A){return C(jx)*pow(Q(Er,A)/mp,2.)*FS2/16;}
double FNN_46(double jx, double FS2, double Er, double A){return C(jx)*pow(Q(Er,A)/mp,2.)*FS2/8.;}

double NaFnn_46(double jx, double Er){return FNN_46(jx,NaFS2nn(Y(Er,ANa)),Er,ANa);}
double NaFpp_46(double jx, double Er){return FNN_46(jx,NaFS2pp(Y(Er,ANa)),Er,ANa);}
double NaFpn_46(double jx, double Er){return FNN_46(jx,NaFS2np(Y(Er,ANa)),Er,ANa);}
double NaFnp_46(double jx, double Er){return FNN_46(jx,NaFS2np(Y(Er,ANa)),Er,ANa);}

double IFnn_46(double jx, double Er){return FNN_46(jx,IFS2nn(Y(Er,AI)),Er,AI);}
double IFpp_46(double jx, double Er){return FNN_46(jx,IFS2pp(Y(Er,AI)),Er,AI);}
double IFpn_46(double jx, double Er){return FNN_46(jx,IFS2np(Y(Er,AI)),Er,AI);}
double IFnp_46(double jx, double Er){return FNN_46(jx,IFS2np(Y(Er,AI)),Er,AI);}

// mathematica:
// 23Na
double NaF00_46(double jx, double Er){return FNN_46(jx,NaF00_S2(Y(Er,ANa)),Er,ANa);}
double NaF11_46(double jx, double Er){return FNN_46(jx,NaF11_S2(Y(Er,ANa)),Er,ANa);}
double NaF01_46(double jx, double Er){return FNN_46(jx,NaF01_S2(Y(Er,ANa)),Er,ANa);}
double NaF10_46(double jx, double Er){return FNN_46(jx,NaF10_S2(Y(Er,ANa)),Er,ANa);}

// // 23Na
// double NaF00_46(double jx, double Er){return 0.25*(NaFnn_46(jx,Er)+NaFpp_46(jx,Er)+NaFpn_46(jx,Er)+NaFnp_46(jx,Er));}
// double NaF11_46(double jx, double Er){return 0.25*(NaFnn_46(jx,Er)+NaFpp_46(jx,Er)-NaFpn_46(jx,Er)-NaFnp_46(jx,Er));}
// double NaF01_46(double jx, double Er){return 0.25*(-NaFnn_46(jx,Er)+NaFpp_46(jx,Er)-NaFpn_46(jx,Er)+NaFnp_46(jx,Er));}
// double NaF10_46(double jx, double Er){return 0.25*(-NaFnn_46(jx,Er)+NaFpp_46(jx,Er)+NaFpn_46(jx,Er)-NaFnp_46(jx,Er));}

// 23Na total F46
double NaF46(double jx, double Er){
    double term1 = c0_4*c0_6*NaF00_46(jx,Er);
    double term2 = c1_4*c1_6*NaF11_46(jx,Er);
    double term3 = c0_4*c1_6*NaF01_46(jx,Er);
    double term4 = c1_4*c0_6*NaF10_46(jx,Er);
    return term1+term2+term3+term4;
}

// mathematica:
// 127I
double IF00_46(double jx, double Er){return FNN_46(jx,IF00_S2(Y(Er,AI)),Er,AI);}
double IF11_46(double jx, double Er){return FNN_46(jx,IF11_S2(Y(Er,AI)),Er,AI);}
double IF01_46(double jx, double Er){return FNN_46(jx,IF01_S2(Y(Er,AI)),Er,AI);}
double IF10_46(double jx, double Er){return FNN_46(jx,IF10_S2(Y(Er,AI)),Er,AI);}

// // 127I
// double IF00_46(double jx, double Er){return 0.25*(IFnn_46(jx,Er)+IFpp_46(jx,Er)+IFpn_46(jx,Er)+IFnp_46(jx,Er));}
// double IF11_46(double jx, double Er){return 0.25*(IFnn_46(jx,Er)+IFpp_46(jx,Er)-IFpn_46(jx,Er)-IFnp_46(jx,Er));}
// double IF01_46(double jx, double Er){return 0.25*(-IFnn_46(jx,Er)+IFpp_46(jx,Er)-IFpn_46(jx,Er)+IFnp_46(jx,Er));}
// double IF10_46(double jx, double Er){return 0.25*(-IFnn_46(jx,Er)+IFpp_46(jx,Er)+IFpn_46(jx,Er)-IFnp_46(jx,Er));}

// 127I total F46
double IF46(double jx, double Er){
    double term1 = c0_4*c0_6*IF00_46(jx,Er);
    double term2 = c1_4*c1_6*IF11_46(jx,Er);
    double term3 = c0_4*c1_6*IF01_46(jx,Er);
    double term4 = c1_4*c0_6*IF10_46(jx,Er);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O5--------------------------------------------------------
//double gFNN_55(double Mw, double jx, double FM, double FD, double Er, double A, double vmin){return 0.25*C(jx)*(FM*pow(Q(Er,A),2.)/pow(mp,2));}       //old
//double hFNN_55(double Mw, double jx, double FM, double FD, double Er, double A){return 0.25*C(jx)*(pow(Q(Er,A)/mp,2.)*FM);}     //old

double gFNN_55(double Mw, double jx, double FM, double FD, double Er, double A, double vmin){return 0.25*C(jx)*(FD*pow(Q(Er,A)/mp,4.)-FM*pow(vmin*Q(Er,A)/(mp),2.));}
double hFNN_55(double Mw, double jx, double FM, double FD, double Er, double A){return 0.25*C(jx)*(pow(Q(Er,A)/mp,2.)*FM);}


double gNaFnn_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaFMnn(Y(Er,ANa)),NaFDnn(Y(Er,ANa)),Er,ANa,vmin);}
double gNaFpp_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaFMpp(Y(Er,ANa)),NaFDpp(Y(Er,ANa)),Er,ANa,vmin);}
double gNaFnp_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaFMnp(Y(Er,ANa)),NaFDnp(Y(Er,ANa)),Er,ANa,vmin);}

double gIFnn_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IFMnn(Y(Er,AI)),IFDnn(Y(Er,AI)),Er,AI,vmin);}
double gIFpp_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IFMpp(Y(Er,AI)),IFDpp(Y(Er,AI)),Er,AI,vmin);}
double gIFnp_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IFMnp(Y(Er,AI)),IFDnp(Y(Er,AI)),Er,AI,vmin);}

double hNaFnn_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaFMnn(Y(Er,ANa)),NaFDnn(Y(Er,ANa)),Er,ANa);}
double hNaFpp_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaFMpp(Y(Er,ANa)),NaFDpp(Y(Er,ANa)),Er,ANa);}
double hNaFnp_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaFMnp(Y(Er,ANa)),NaFDnp(Y(Er,ANa)),Er,ANa);}

double hIFnn_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IFMnn(Y(Er,AI)),IFDnn(Y(Er,AI)),Er,AI);}
double hIFpp_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IFMpp(Y(Er,AI)),IFDpp(Y(Er,AI)),Er,AI);}
double hIFnp_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IFMnp(Y(Er,AI)),IFDnp(Y(Er,AI)),Er,AI);}

// from the mathematica code:
double gNaF00_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaF00_M(Y(Er,ANa)),NaF00_D(Y(Er,ANa)),Er,ANa,vmin);}
double gNaF11_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaF11_M(Y(Er,ANa)),NaF11_D(Y(Er,ANa)),Er,ANa,vmin);}
double gNaF01_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaF01_M(Y(Er,ANa)),NaF01_D(Y(Er,ANa)),Er,ANa,vmin);}
double gNaF10_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,NaF10_M(Y(Er,ANa)),NaF10_D(Y(Er,ANa)),Er,ANa,vmin);}

double hNaF00_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaF00_M(Y(Er,ANa)),NaF00_D(Y(Er,ANa)),Er,ANa);}
double hNaF11_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaF11_M(Y(Er,ANa)),NaF11_D(Y(Er,ANa)),Er,ANa);}
double hNaF01_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaF01_M(Y(Er,ANa)),NaF01_D(Y(Er,ANa)),Er,ANa);}
double hNaF10_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,NaF10_M(Y(Er,ANa)),NaF10_D(Y(Er,ANa)),Er,ANa);}

// // 23Na
// double gNaF00_55(double jx, double Er, double Mw, double vmin){return 0.25*(gNaFnn_55(jx,Er,Mw,vmin)+gNaFpp_55(jx,Er,Mw,vmin)+2*gNaFnp_55(jx,Er,Mw,vmin));}
// double gNaF11_55(double jx, double Er, double Mw, double vmin){return 0.25*(gNaFnn_55(jx,Er,Mw,vmin)+gNaFpp_55(jx,Er,Mw,vmin)-2*gNaFnp_55(jx,Er,Mw,vmin));}
// double gNaF01_55(double jx, double Er, double Mw, double vmin){return 0.25*(-gNaFnn_55(jx,Er,Mw,vmin)+gNaFpp_55(jx,Er,Mw,vmin));}
// double gNaF10_55(double jx, double Er, double Mw, double vmin){return 0.25*(-gNaFnn_55(jx,Er,Mw,vmin)+gNaFpp_55(jx,Er,Mw,vmin));}
//
// double hNaF00_55(double jx, double Er, double Mw){return 0.25*(hNaFnn_55(jx,Er,Mw)+hNaFpp_55(jx,Er,Mw)+2*hNaFnp_55(jx,Er,Mw));}
// double hNaF11_55(double jx, double Er, double Mw){return 0.25*(hNaFnn_55(jx,Er,Mw)+hNaFpp_55(jx,Er,Mw)-2*hNaFnp_55(jx,Er,Mw));}
// double hNaF01_55(double jx, double Er, double Mw){return 0.25*(-hNaFnn_55(jx,Er,Mw)+hNaFpp_55(jx,Er,Mw));}
// double hNaF10_55(double jx, double Er, double Mw){return 0.25*(-hNaFnn_55(jx,Er,Mw)+hNaFpp_55(jx,Er,Mw));}

// 23Na total F55
double gNaF55(double jx, double Er, double Mw, double vmin){
    double term1 = c0_5*c0_5*gNaF00_55(jx,Er,Mw,vmin);
    double term2 = c1_5*c1_5*gNaF11_55(jx,Er,Mw,vmin);
    double term3 = c0_5*c1_5*gNaF01_55(jx,Er,Mw,vmin);
    double term4 = c1_5*c0_5*gNaF10_55(jx,Er,Mw,vmin);
    return term1+term2+term3+term4;
}

double hNaF55(double jx, double Er, double Mw){
    double term1 = c0_5*c0_5*hNaF00_55(jx,Er,Mw);
    double term2 = c1_5*c1_5*hNaF11_55(jx,Er,Mw);
    double term3 = c0_5*c1_5*hNaF01_55(jx,Er,Mw);
    double term4 = c1_5*c0_5*hNaF10_55(jx,Er,Mw);
    return term1+term2+term3+term4;
}

// from the mathematica code:
double gIF00_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IF00_M(Y(Er,AI)),IF00_D(Y(Er,AI)),Er,AI,vmin);}
double gIF11_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IF11_M(Y(Er,AI)),IF11_D(Y(Er,AI)),Er,AI,vmin);}
double gIF01_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IF01_M(Y(Er,AI)),IF01_D(Y(Er,AI)),Er,AI,vmin);}
double gIF10_55(double jx, double Er, double Mw, double vmin){return gFNN_55(Mw,jx,IF10_M(Y(Er,AI)),IF10_D(Y(Er,AI)),Er,AI,vmin);}

double hIF00_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IF00_M(Y(Er,AI)),IF00_D(Y(Er,AI)),Er,AI);}
double hIF11_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IF11_M(Y(Er,AI)),IF11_D(Y(Er,AI)),Er,AI);}
double hIF01_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IF01_M(Y(Er,AI)),IF01_D(Y(Er,AI)),Er,AI);}
double hIF10_55(double jx, double Er, double Mw){return hFNN_55(Mw,jx,IF10_M(Y(Er,AI)),IF10_D(Y(Er,AI)),Er,AI);}

// // 127I
// double gIF00_55(double jx, double Er, double Mw, double vmin){return 0.25*(gIFnn_55(jx,Er,Mw,vmin)+gIFpp_55(jx,Er,Mw,vmin)+2*gIFnp_55(jx,Er,Mw,vmin));}
// double gIF11_55(double jx, double Er, double Mw, double vmin){return 0.25*(gIFnn_55(jx,Er,Mw,vmin)+gIFpp_55(jx,Er,Mw,vmin)-2*gIFnp_55(jx,Er,Mw,vmin));}
// double gIF01_55(double jx, double Er, double Mw, double vmin){return 0.25*(-gIFnn_55(jx,Er,Mw,vmin)+gIFpp_55(jx,Er,Mw,vmin));}
// double gIF10_55(double jx, double Er, double Mw, double vmin){return 0.25*(-gIFnn_55(jx,Er,Mw,vmin)+gIFpp_55(jx,Er,Mw,vmin));}
//
// double hIF00_55(double jx, double Er, double Mw){return 0.25*(hIFnn_55(jx,Er,Mw)+hIFpp_55(jx,Er,Mw)+2*hIFnp_55(jx,Er,Mw));}
// double hIF11_55(double jx, double Er, double Mw){return 0.25*(hIFnn_55(jx,Er,Mw)+hIFpp_55(jx,Er,Mw)-2*hIFnp_55(jx,Er,Mw));}
// double hIF01_55(double jx, double Er, double Mw){return 0.25*(-hIFnn_55(jx,Er,Mw)+hIFpp_55(jx,Er,Mw));}
// double hIF10_55(double jx, double Er, double Mw){return 0.25*(-hIFnn_55(jx,Er,Mw)+hIFpp_55(jx,Er,Mw));}

// 127I total F55
double gIF55(double jx, double Er, double Mw, double vmin){
    double term1 = c0_5*c0_5*gIF00_55(jx,Er,Mw,vmin);
    double term2 = c1_5*c1_5*gIF11_55(jx,Er,Mw,vmin);
    double term3 = c0_5*c1_5*gIF01_55(jx,Er,Mw,vmin);
    double term4 = c1_5*c0_5*gIF10_55(jx,Er,Mw,vmin);
    return term1+term2+term3+term4;
}

double hIF55(double jx, double Er, double Mw){
    double term1 = c0_5*c0_5*hIF00_55(jx,Er,Mw);
    double term2 = c1_5*c1_5*hIF11_55(jx,Er,Mw);
    double term3 = c0_5*c1_5*hIF01_55(jx,Er,Mw);
    double term4 = c1_5*c0_5*hIF10_55(jx,Er,Mw);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O6--------------------------------------------------------
double FNN_66(double jx, double FS2, double Er, double A){return C(jx)*pow(Q(Er,A)/mp,4.)*FS2/16.;}

double NaFnn_66(double jx, double Er){return FNN_66(jx,NaFS2nn(Y(Er,ANa)),Er,ANa);}
double NaFpp_66(double jx, double Er){return FNN_66(jx,NaFS2pp(Y(Er,ANa)),Er,ANa);}
double NaFnp_66(double jx, double Er){return FNN_66(jx,NaFS2np(Y(Er,ANa)),Er,ANa);}

double IFnn_66(double jx, double Er){return FNN_66(jx,IFS2nn(Y(Er,AI)),Er,AI);}
double IFpp_66(double jx, double Er){return FNN_66(jx,IFS2pp(Y(Er,AI)),Er,AI);}
double IFnp_66(double jx, double Er){return FNN_66(jx,IFS2np(Y(Er,AI)),Er,AI);}

// mathematica:
// 23Na
// NB S2 already defined
double NaF00_66(double jx, double Er){return FNN_66(jx,NaF00_S2(Y(Er,ANa)),Er,ANa);}
double NaF11_66(double jx, double Er){return FNN_66(jx,NaF11_S2(Y(Er,ANa)),Er,ANa);}
double NaF01_66(double jx, double Er){return FNN_66(jx,NaF01_S2(Y(Er,ANa)),Er,ANa);}
double NaF10_66(double jx, double Er){return FNN_66(jx,NaF10_S2(Y(Er,ANa)),Er,ANa);}

// // 23Na
// double NaF00_66(double jx, double Er){return 0.25*(NaFnn_66(jx,Er)+NaFpp_66(jx,Er)+2*NaFnp_66(jx,Er));}
// double NaF11_66(double jx, double Er){return 0.25*(NaFnn_66(jx,Er)+NaFpp_66(jx,Er)-2*NaFnp_66(jx,Er));}
// double NaF01_66(double jx, double Er){return 0.25*(-NaFnn_66(jx,Er)+NaFpp_66(jx,Er));}
// double NaF10_66(double jx, double Er){return 0.25*(-NaFnn_66(jx,Er)+NaFpp_66(jx,Er));}

// 23Na total F66
double NaF66(double jx, double Er){
    double term1 = c0_6*c0_6*NaF00_66(jx,Er);
    double term2 = c1_6*c1_6*NaF11_66(jx,Er);
    double term3 = c0_6*c1_6*NaF01_66(jx,Er);
    double term4 = c1_6*c0_6*NaF10_66(jx,Er);
    return term1+term2+term3+term4;
}

// mathematica:
// 127I
// NB S2 already defined
double IF00_66(double jx, double Er){return FNN_66(jx,IF00_S2(Y(Er,AI)),Er,AI);}
double IF11_66(double jx, double Er){return FNN_66(jx,IF11_S2(Y(Er,AI)),Er,AI);}
double IF01_66(double jx, double Er){return FNN_66(jx,IF01_S2(Y(Er,AI)),Er,AI);}
double IF10_66(double jx, double Er){return FNN_66(jx,IF10_S2(Y(Er,AI)),Er,AI);}

// // 127I
// double IF00_66(double jx, double Er){return 0.25*(IFnn_66(jx,Er)+IFpp_66(jx,Er)+2*IFnp_66(jx,Er));}
// double IF11_66(double jx, double Er){return 0.25*(IFnn_66(jx,Er)+IFpp_66(jx,Er)-2*IFnp_66(jx,Er));}
// double IF01_66(double jx, double Er){return 0.25*(-IFnn_66(jx,Er)+IFpp_66(jx,Er));}
// double IF10_66(double jx, double Er){return 0.25*(-IFnn_66(jx,Er)+IFpp_66(jx,Er));}

// 127I total F66
double IF66(double jx, double Er){
    double term1 = c0_6*c0_6*IF00_66(jx,Er);
    double term2 = c1_6*c1_6*IF11_66(jx,Er);
    double term3 = c0_6*c1_6*IF01_66(jx,Er);
    double term4 = c1_6*c0_6*IF10_66(jx,Er);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O7-------------------------------------------------------
double gFNN_77(double jx, double FS1, double Er, double A, double Mw, double vmin){return -vmin*vmin*FS1/(8.);}
double hFNN_77(double jx, double FS1, double Er, double A, double Mw){return FS1/(8.);}

double gNaFnn_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaFS1nn(Y(Er,ANa)),Er,ANa,Mw,vmin);}
double gNaFpp_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaFS1pp(Y(Er,ANa)),Er,ANa,Mw,vmin);}
double gNaFnp_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaFS1np(Y(Er,ANa)),Er,ANa,Mw,vmin);}

double hNaFnn_77(double jx, double Er, double Mw){return hFNN_77(jx,NaFS1nn(Y(Er,ANa)),Er,ANa,Mw);}
double hNaFpp_77(double jx, double Er, double Mw){return hFNN_77(jx,NaFS1pp(Y(Er,ANa)),Er,ANa,Mw);}
double hNaFnp_77(double jx, double Er, double Mw){return hFNN_77(jx,NaFS1np(Y(Er,ANa)),Er,ANa,Mw);}

double gIFnn_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IFS1nn(Y(Er,AI)),Er,AI,Mw,vmin);}
double gIFpp_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IFS1pp(Y(Er,AI)),Er,AI,Mw,vmin);}
double gIFnp_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IFS1np(Y(Er,AI)),Er,AI,Mw,vmin);}

double hIFnn_77(double jx, double Er, double Mw){return hFNN_77(jx,IFS1nn(Y(Er,AI)),Er,AI,Mw);}
double hIFpp_77(double jx, double Er, double Mw){return hFNN_77(jx,IFS1pp(Y(Er,AI)),Er,AI,Mw);}
double hIFnp_77(double jx, double Er, double Mw){return hFNN_77(jx,IFS1np(Y(Er,AI)),Er,AI,Mw);}

// from the mathematica code:
double gNaF00_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaF00_S1(Y(Er,ANa)),Er,ANa,Mw,vmin);}
double gNaF11_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaF11_S1(Y(Er,ANa)),Er,ANa,Mw,vmin);}
double gNaF01_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaF01_S1(Y(Er,ANa)),Er,ANa,Mw,vmin);}
double gNaF10_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,NaF10_S1(Y(Er,ANa)),Er,ANa,Mw,vmin);}
double hNaF00_77(double jx, double Er, double Mw){return hFNN_77(jx,NaF00_S1(Y(Er,ANa)),Er,ANa,Mw);}
double hNaF11_77(double jx, double Er, double Mw){return hFNN_77(jx,NaF11_S1(Y(Er,ANa)),Er,ANa,Mw);}
double hNaF01_77(double jx, double Er, double Mw){return hFNN_77(jx,NaF01_S1(Y(Er,ANa)),Er,ANa,Mw);}
double hNaF10_77(double jx, double Er, double Mw){return hFNN_77(jx,NaF10_S1(Y(Er,ANa)),Er,ANa,Mw);}

// // 23Na
// double gNaF00_77(double jx, double Er, double Mw, double vmin){return 0.25*(gNaFnn_77(jx,Er,Mw,vmin)+gNaFpp_77(jx,Er,Mw,vmin)+2*gNaFnp_77(jx,Er,Mw,vmin));}
// double gNaF11_77(double jx, double Er, double Mw, double vmin){return 0.25*(gNaFnn_77(jx,Er,Mw,vmin)+gNaFpp_77(jx,Er,Mw,vmin)-2*gNaFnp_77(jx,Er,Mw,vmin));}
// double gNaF01_77(double jx, double Er, double Mw, double vmin){return 0.25*(-gNaFnn_77(jx,Er,Mw,vmin)+gNaFpp_77(jx,Er,Mw,vmin));}
// double gNaF10_77(double jx, double Er, double Mw, double vmin){return 0.25*(-gNaFnn_77(jx,Er,Mw,vmin)+gNaFpp_77(jx,Er,Mw,vmin));}
//
//
// double hNaF00_77(double jx, double Er, double Mw){return 0.25*(hNaFnn_77(jx,Er,Mw)+hNaFpp_77(jx,Er,Mw)+2*hNaFnp_77(jx,Er,Mw));}
// double hNaF11_77(double jx, double Er, double Mw){return 0.25*(hNaFnn_77(jx,Er,Mw)+hNaFpp_77(jx,Er,Mw)-2*hNaFnp_77(jx,Er,Mw));}
// double hNaF01_77(double jx, double Er, double Mw){return 0.25*(-hNaFnn_77(jx,Er,Mw)+hNaFpp_77(jx,Er,Mw));}
// double hNaF10_77(double jx, double Er, double Mw){return 0.25*(-hNaFnn_77(jx,Er,Mw)+hNaFpp_77(jx,Er,Mw));}


// 23Na total F77
double gNaF77(double jx, double Er, double Mw, double vmin){
    double term1 = c0_7*c0_7*gNaF00_77(jx,Er,Mw,vmin);
    double term2 = c1_7*c1_7*gNaF11_77(jx,Er,Mw,vmin);
    double term3 = c0_7*c1_7*gNaF01_77(jx,Er,Mw,vmin);
    double term4 = c1_7*c0_7*gNaF10_77(jx,Er,Mw,vmin);
    return term1+term2+term3+term4;
}

double hNaF77(double jx, double Er, double Mw){
    double term1 = c0_7*c0_7*hNaF00_77(jx,Er,Mw);
    double term2 = c1_7*c1_7*hNaF11_77(jx,Er,Mw);
    double term3 = c0_7*c1_7*hNaF01_77(jx,Er,Mw);
    double term4 = c1_7*c0_7*hNaF10_77(jx,Er,Mw);
    return term1+term2+term3+term4;
}

// from the mathematica code:
double gIF00_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IF00_S1(Y(Er,AI)),Er,AI,Mw,vmin);}
double gIF11_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IF11_S1(Y(Er,AI)),Er,AI,Mw,vmin);}
double gIF01_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IF01_S1(Y(Er,AI)),Er,AI,Mw,vmin);}
double gIF10_77(double jx, double Er, double Mw, double vmin){return gFNN_77(jx,IF10_S1(Y(Er,AI)),Er,AI,Mw,vmin);}
double hIF00_77(double jx, double Er, double Mw){return hFNN_77(jx,IF00_S1(Y(Er,AI)),Er,AI,Mw);}
double hIF11_77(double jx, double Er, double Mw){return hFNN_77(jx,IF11_S1(Y(Er,AI)),Er,AI,Mw);}
double hIF01_77(double jx, double Er, double Mw){return hFNN_77(jx,IF01_S1(Y(Er,AI)),Er,AI,Mw);}
double hIF10_77(double jx, double Er, double Mw){return hFNN_77(jx,IF10_S1(Y(Er,AI)),Er,AI,Mw);}

// // 127I
// double gIF00_77(double jx, double Er, double Mw, double vmin){return 0.25*(gIFnn_77(jx,Er,Mw,vmin)+gIFpp_77(jx,Er,Mw,vmin)+2*gIFnp_77(jx,Er,Mw,vmin))               ;}
// double gIF11_77(double jx, double Er, double Mw, double vmin){return 0.25*(gIFnn_77(jx,Er,Mw,vmin)+gIFpp_77(jx,Er,Mw,vmin)-2*gIFnp_77(jx,Er,Mw,vmin));}
// double gIF01_77(double jx, double Er, double Mw, double vmin){return 0.25*(-gIFnn_77(jx,Er,Mw,vmin)+gIFpp_77(jx,Er,Mw,vmin));}
// double gIF10_77(double jx, double Er, double Mw, double vmin){return 0.25*(-gIFnn_77(jx,Er,Mw,vmin)+gIFpp_77(jx,Er,Mw,vmin));}
//
// double hIF00_77(double jx, double Er, double Mw){return 0.25*(hIFnn_77(jx,Er,Mw)+hIFpp_77(jx,Er,Mw)+2*hIFnp_77(jx,Er,Mw));}
// double hIF11_77(double jx, double Er, double Mw){return 0.25*(hIFnn_77(jx,Er,Mw)+hIFpp_77(jx,Er,Mw)-2*hIFnp_77(jx,Er,Mw));}
// double hIF01_77(double jx, double Er, double Mw){return 0.25*(-hIFnn_77(jx,Er,Mw)+hIFpp_77(jx,Er,Mw));}
// double hIF10_77(double jx, double Er, double Mw){return 0.25*(-hIFnn_77(jx,Er,Mw)+hIFpp_77(jx,Er,Mw));}

// 127I total F77
double gIF77(double jx, double Er, double Mw, double vmin){
    double term1 = c0_7*c0_7*gIF00_77(jx,Er,Mw,vmin);
    double term2 = c1_7*c1_7*gIF11_77(jx,Er,Mw,vmin);
    double term3 = c0_7*c1_7*gIF01_77(jx,Er,Mw,vmin);
    double term4 = c1_7*c0_7*gIF10_77(jx,Er,Mw,vmin);
    return term1+term2+term3+term4;
}

double hIF77(double jx, double Er, double Mw){
    double term1 = c0_7*c0_7*hIF00_77(jx,Er,Mw);
    double term2 = c1_7*c1_7*hIF11_77(jx,Er,Mw);
    double term3 = c0_7*c1_7*hIF01_77(jx,Er,Mw);
    double term4 = c1_7*c0_7*hIF10_77(jx,Er,Mw);
    return term1+term2+term3+term4;
}

//--------------------------------------------------------Operator O10--------------------------------------------------------
double FNN_10(double jx, double FS2, double Er, double A){return pow(Q_conv*Q(Er,A),2.)*FS2/4.;}

double NaFnn_10(double jx, double Er){return FNN_10(jx,NaFS2nn(Y(Er,ANa)),Er,ANa);}
double NaFpp_10(double jx, double Er){return FNN_10(jx,NaFS2pp(Y(Er,ANa)),Er,ANa);}
double NaFnp_10(double jx, double Er){return FNN_10(jx,NaFS2np(Y(Er,ANa)),Er,ANa);}

double IFnn_10(double jx, double Er){return FNN_10(jx,IFS2nn(Y(Er,AI)),Er,AI);}
double IFpp_10(double jx, double Er){return FNN_10(jx,IFS2pp(Y(Er,AI)),Er,AI);}
double IFnp_10(double jx, double Er){return FNN_10(jx,IFS2np(Y(Er,AI)),Er,AI);}

// 23Na
double NaF00_10(double jx, double Er){return 0.25*(NaFnn_10(jx, Er)+NaFpp_10(jx, Er)+2*NaFnp_10(jx, Er));}
double NaF11_10(double jx, double Er){return 0.25*(NaFnn_10(jx, Er)+NaFpp_10(jx, Er)-2*NaFnp_10(jx, Er));}
double NaF01_10(double jx, double Er){return 0.25*(-NaFnn_10(jx, Er)+NaFpp_10(jx, Er));}
double NaF10_10(double jx, double Er){return 0.25*(-NaFnn_10(jx, Er)+NaFpp_10(jx, Er));}

// 23Na total F10
double NaF10(double jx, double Er){
    double term1 = c0_10*c0_10*NaF00_10(jx,Er);
    double term2 = c1_10*c1_10*NaF11_10(jx,Er);
    double term3 = c0_10*c1_10*NaF01_10(jx,Er);
    double term4 = c1_10*c0_10*NaF10_10(jx,Er);
    return term1+term2+term3+term4;
}

// 127I
double IF00_10(double jx, double Er){return 0.25*(IFnn_10(jx, Er)+IFpp_10(jx, Er)+2*IFnp_10(jx, Er));}
double IF11_10(double jx, double Er){return 0.25*(IFnn_10(jx, Er)+IFpp_10(jx, Er)-2*IFnp_10(jx, Er));}
double IF01_10(double jx, double Er){return 0.25*(-IFnn_10(jx, Er)+IFpp_10(jx, Er));}
double IF10_10(double jx, double Er){return 0.25*(-IFnn_10(jx, Er)+IFpp_10(jx, Er));}

// 127I total F10
double IF10(double jx, double Er){
    double term1 = c0_10*c0_10*IF00_10(jx,Er);
    double term2 = c1_10*c1_10*IF11_10(jx,Er);
    double term3 = c0_10*c1_10*IF01_10(jx,Er);
    double term4 = c1_10*c0_10*IF10_10(jx,Er);
    return term1+term2+term3+term4;
}

//-------------------------------------------------------- HELM FORM FAC --------------------------------------------------------
//For the purposes of getting this integration going properly, lets go back to the OG form factors for a spin indep so we know when we have recreated
double c1(double A){return 1.23*pow(A,(1./3.)) -0.6;}
double rn(double A){return pow(( c1(A)*c1(A) + (7./3.)*0.52*0.52*TMath::Pi()*TMath::Pi() - 5.*(0.9)*(0.9)) , 0.5);}
double HelmF(double Er, double A){ return 3.* (exp(-0.5*(Q_conv*Q_conv*Q(Er,A)*Q(Er,A)*0.9*0.9))) * ( sin(Q_conv*Q(Er,A)*rn(A)) - (Q_conv*Q(Er,A)*rn(A)*cos(Q_conv*Q(Er,A)*rn(A))) )/(pow(Q_conv*Q(Er,A)*rn(A),3.));}

//-------------------------------------------------------- Total Form factor sum --------------------------------------------------------

double gNa_total_form(double jx, double Er, double Mw, double vmin){
 return (NaF11(Er)+NaF44(jx,Er)+NaF45(jx,Er)+NaF46(jx,Er)+gNaF55(jx,Er,Mw,vmin)+NaF66(jx,Er)+gNaF77(jx,Er,Mw,vmin)+NaF10(jx,Er));//+c_Helm*pow(HelmF(Er,ANa),2);
}

double hNa_total_form(double jx, double Er, double Mw){
 return (hNaF55(jx,Er,Mw)+hNaF77(jx,Er,Mw));
}

double gI_total_form(double jx, double Er, double Mw, double vmin){
 return (IF11(Er)+IF44(jx,Er)+IF45(jx,Er)+IF46(jx,Er)+gIF55(jx,Er,Mw,vmin)+IF66(jx,Er)+gIF77(jx,Er,Mw,vmin)+IF10(jx,Er));//+c_Helm*pow(HelmF(Er,AI),2);
}

double hI_total_form(double jx, double Er, double Mw){
 return (hIF55(jx,Er,Mw)+hIF77(jx,Er,Mw));
}



#endif
