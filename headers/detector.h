// This header file define a number of detector dependent terms, such as quenching factor, resolution function etc.
// As it stand, this is equipped for NaI detectors SABRE and DAMA/LIBRA, assuming they have the same efficiency and resolution.


#ifndef qf_H
#define qf_H

const double QvalueI=0.09;					//Iodine quenching factor
const double DAMAQuenchgNA=0.3;				//DAMA na quenching value
bool poiss;
double threshval;

//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                                      Q u e n c h i n g   F a c t o r
//--------------------------------------------------------------------------------------------------------------------------------------------------

// Select which measured quenching factor will be used:
// 0=DAMA quenching values, 1= Princeton quenching values, 2=University of Texas quenching values, 
//int whichNaQuench =1;


//Quenching factor measured by Princeton as a function of the energy (nuclear recoil):
const double PrincetonQuench[11]={0.133,0.129,0.162,0.159,0.16,0.168,0.171,0.188,0.191,0.204,0.207};
const double PrincetonErrorQuench[11]={0.018,0.014,0.012,0.019,0.010,0.009,0.010,0.008,0.011,0.008,0.010};
const double PrincetonEnergy_nr[11]={5.7,8.8,9.1,14.3,15.,19.4,24.9,29.,33.3,43.,51.8};
const double PrincetonErrorEnergy_nr[11]={0.7,1.2,1.2,2.4,1.4,1.6,2.4,1.9,2.8,2.2,2.6};

//Quenching factor measured by the University of Texas as a function of the energy (nuclear recoil):
const double TexasQuench[11]={0.08,0.056,0.068,0.08,0.105,0.125,0.1425,0.1525,0.14,0.18,0.18};
const double TexasErrorQuench[11]={0.02,0.02,0.015,0.03,0.04,0.04,0.03,0.02,0.04,0.03,0.04};
const double TexasEnergy_nr[11]={7.31,8.39,9.46,11.6,17.7,20.8,30.2,31.4,34.4,39.7,47.7};
const double TexasErrorEnergy_nr[11]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

//ANU
const double ANUQuench[16]={0.160,0.198,0.201,0.210,0.235,0.216,0.194,0.223,0.233,0.220,0.252,0.237,0.250,0.260,0.309,0.305};
const double ANUErrorQuench[16]={0.013,0.009,0.008,0.007,0.008,0.011,0.007,0.009,0.007,0.006,0.010,0.009,0.017,0.018,0.006,0.006};
const double ANUEnergy_nr[16]={36,58,65,65,71,79,86,96,107,107,147,147,179,179,360,401};
const double ANUErrorEnergy_nr[16]={5,1,8,8,10,3,12,3,7,7,5,5,4,4,24,15};

TF1 *ENuclearRecoilfunction;
TF1 *QuenchNuclearRecoilfunction;


void setquenching(int Q){
//Set the quenching values selected in the configuration file:
  
  //if(Q!=0) ENuclearRecoilfunction = new TF1("ENuclearRecoilfunction","((-[0]+pow([0]*[0]+4.*[1]*x,0.5))/(2.*[1]))",0.,100.); // E_nr=E_er/Q(E_er) with an energy dependent value of quenching (get E_nr from E_er = Q(E_nr)*E_nr, where Q(E_nr)= fit experimental data with [0]+[1]*E_nr)
  if(Q!=0) ENuclearRecoilfunction = new TF1("ENuclearRecoilfunction","(pow(x/[0],1./(1.+[1])))",0.,100.); // E_nr=E_er/Q(E_er) with an energy dependent value of quenching (get E_nr from E_er = Q(E_nr)*E_nr, where Q(E_nr)= fit experimental data with a power law)
  else if(Q==0){
    ENuclearRecoilfunction = new TF1("ENuclearRecoilfunction","x/[0]",0.,100.); // E_nr=E_er/Q(E_er) with a costant value of quenching (DAMA: Q(E_er)=0.3 energy independent)
    ENuclearRecoilfunction->SetParameter(0,DAMAQuenchgNA);
  }
  TGraphErrors *QuenchingValuesGraph;
  if(Q==1){QuenchingValuesGraph = new TGraphErrors(11,PrincetonEnergy_nr,PrincetonQuench,PrincetonErrorEnergy_nr,PrincetonErrorQuench);}
  if(Q==2) {QuenchingValuesGraph = new TGraphErrors(11,TexasEnergy_nr,TexasQuench,TexasErrorEnergy_nr,TexasErrorQuench);}
  if(Q==3) {QuenchingValuesGraph = new TGraphErrors(16,ANUEnergy_nr,ANUQuench,ANUErrorEnergy_nr,ANUErrorQuench);}
  if(Q!=0){
    // QuenchNuclearRecoilfunction = new TF1("QuenchNuclearRecoilfunction","[0]+[1]*x",0.,100.);
    QuenchNuclearRecoilfunction = new TF1("QuenchNuclearRecoilfunction","[0]*pow(x,[1])",0.,100.);
    QuenchingValuesGraph->Fit("QuenchNuclearRecoilfunction","QR");
    ENuclearRecoilfunction->SetParameter(0,QuenchNuclearRecoilfunction->GetParameter(0));
    ENuclearRecoilfunction->SetParameter(1,QuenchNuclearRecoilfunction->GetParameter(1));
  }
  //delete QuenchingValuesGraph;
  //delete QuenchNuclearRecoilfunction;
}



double Quench(double A, double E2){
	if(A==22.99){
        return (1E3)*ENuclearRecoilfunction->Eval(E2);
		}
	else return E2*1.E3/QvalueI;
}

double DerivQuench(double A, double E2){
	if(A==22.99){
	     return ENuclearRecoilfunction->Derivative(E2);
		}
	else return 1./QvalueI;
}


//--------------------------------------------------------------------------------------------------------------------------------------------------
//                                                  R e s o l u t i o n   f u n c t i o n
//--------------------------------------------------------------------------------------------------------------------------------------------------

// Efficiency of a detector
double Threshold(double E){
    if (E>=8.) return 1.;
    // else return 0.075*E+0.4;          // value for DAMA-phase1
    else return 0.0429*E+0.657;         // value for DAMA-phase2
}  
    

// Resolution function
double DE(double Eee, int x){ 
  if(x==0){ //DAMA resolution function 
    return ((0.488*pow(Eee,0.5))+(0.0091*Eee)); 
  } 
  if(x==1){ //COSINE resolution function 
    return ((0.3171*pow(Eee,0.5))+(0.008189*Eee)); 
  } 
  if(x==2){ //SABRE resolution function 
    return (0.632*pow(Eee,0.5)); 
  } 
  else{
    cout<<"invalid resolution function choice";
    return 0;
  }
}

double ResFunction(double E, double Eee, int x) { // Resolution function 
  double A = 1./(sqrt(2.*TMath::Pi())*DE(Eee,x)); 
  double eps;
  if(x==0){eps = 8.0;} //Light yield of DAMA
  if(x==1){eps = 15;} //Light yield of COSINE
  if(x==2){eps = 11.0;} //Light yield of NaI-33
  double poiss_av = Eee*eps;
  double poiss_x = E*eps;
  if(poiss==false){
    return A*exp(-0.5*pow((E- Eee)/DE(Eee,x), 2.));
  }
  else{
    if(Eee<threshval){  
    return eps*TMath::Poisson(poiss_x,poiss_av);
    }
    else{return A*exp(-0.5*pow((E- Eee)/DE(Eee,x), 2.));}
  }
}



#endif
