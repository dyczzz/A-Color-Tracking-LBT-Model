///////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Linear Boltzmann Transport Model
//
//  Updated on 6/1/2017: clean up code and extract input parameters out into a input file
//
//  Run by ./LBT parameter_file input_parton_list HQ_output light_positive_output light_negative_output
//  If initHardFlag==0, no input_parton_list
//  If heavyOut==0, no HQ_output
//  If lightOut==0, no light_positive/negative_output
//
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<fstream>
#include<iostream>
#include<iomanip>

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<cstdlib>

using namespace std;

#include <math.h>
#include <vector>

// declare functions for using OSU/CCNU hydro

#include<cstring>

#include "common.h"
#include "setParameter.cpp"
#include "LBT.cpp"
#include "fnc.cpp"

int main(int argc, char* argv[]){  					

    double initPT,rangePT;
    string parameterFile;
    ifstream fpList;
    ofstream outHQ,positive,negative;

    ofstream ffqperp;                                       
    ffqperp.open("qperp.dat");

    if(argc<2) {
        cout<<"Please add file name for input parameter ./LBT parameter_file ..." << endl;
        exit(EXIT_FAILURE);
    }
    parameterFile = argv[1];
    setParameter(parameterFile); // re-set parameters from input parameter file	

    if(checkParameter(argc)==0) { // check whether the input parameters are all correct
        cout << "Parameter check passed" << endl;
    } else {
        cout << "Parameter check failed" << endl;
        exit(EXIT_FAILURE);
    }
	
    if(vacORmed==1) read_tables(); // initialize various tables

//...define derived quantities
    temp00=temp0;		
    dt=dtau;
    timend=tauend;			
    time0=tau0;	  	
//...alphas	
    alphas=alphas0(Kalphas,temp0);
//...Debye Mass square
    qhat0=DebyeMass2(Kqhat0,alphas,temp0);	
	 
//...multiple scattering process


//...initialize the random number generator
    srand((unsigned)time(NULL));
    NUM1=-1*rand();
//    NUM1=-3;
    struct tm *local_start;
    time_t time_start;
    time_start=time(NULL);
    local_start=localtime(&time_start);
    
    char buf1[80];
    strftime(buf1,80,"Current Time: %Y-%m-%d %H:%M:%S",local_start);
    cout << "the program starts at:" <<endl;
    cout << buf1 << endl;

//...input and output files
    if(initHardFlag==1) { // initialize hard partons within LBT
        if(heavyOut==1 && lightOut==1) { // both heavy and light output
            outHQ.open(argv[2]);
            positive.open(argv[3]);
            negative.open(argv[4]);
        } else if(heavyOut==1) { // only heavy output
            outHQ.open(argv[2]);
        } else if(lightOut==1) { // only light output
            positive.open(argv[2]);
            negative.open(argv[3]);
        } else {
            cout << "Warning: no output file will be generated!" << endl;
        }
    } else { // need initial parton list 
        fpList.open(argv[2]);
        if(heavyOut==1 && lightOut==1) { // both heavy and light output
            outHQ.open(argv[3]);
            positive.open(argv[4]);
            negative.open(argv[5]);
        } else if(heavyOut==1) { // only heavy output
            outHQ.open(argv[3]);
        } else if(lightOut==1) { // only light output
            positive.open(argv[3]);
            negative.open(argv[4]);
        } else {
            cout << "Warning: no output file will be generated!" << endl;
        }
    }

//...Begin event + time loop for parton evolution	  

    int numEvent = 0;	   
    long num1;
    for(int n=1; ; ++n) {
        num1=NUM1;
        if(initHardFlag==1) {  // initialize within LBT
            if(numEvent>=ncall) break;
        } else if(initHardFlag==2) { // initialize by reading particle list
            if(fpList.eof()) break;
            if(numEvent>=ncall) break;
        } // no other possibility, already checked by checkParameter function

//        eGluon=0.0;       
//        nGluon=0.0;

        jetClean();
        double eweight=0.0;

//...initilization of jet parton

        if(initHardFlag==1) {  // initialize within LBT
            jetInitialize(numInitXY); // initialize jet partons
        } else { // initialize by reading particle list
            int dummyInt;
            fpList >> dummyInt;
            if(fpList.eof()) break;
            fpList >> nj;

            // read particle information from file, may move to a function outside in future
            double EiTot=0.0;
            for(int i=1; i<=nj; i=i+1) {

               if(colorFlag==1) { // SC:color

                    fpList >> dummyInt >> KATT1[i] >> P[1][i] >> P[2][i] >> P[3][i] >> P[0][i] >> P[4][i] >> Vfrozen[1][i] >> Vfrozen[2][i] >> Vfrozen[3][i] >> Vfrozen[0][i] >> color[i] >> anticolor[i];
                    if(abs(color[i]) > maxColor) maxColor = abs(color[i]);
                    if(abs(anticolor[i] > maxColor)) maxColor = abs(anticolor[i]);

                } else {
                    fpList >> dummyInt >> KATT1[i] >> P[1][i] >> P[2][i] >> P[3][i] >> P[0][i] >> P[4][i] >> Vfrozen[1][i] >> Vfrozen[2][i] >> Vfrozen[3][i] >> Vfrozen[0][i];
                }

                Vfrozen[1][i]=Vfrozen[1][i]+P[1][i]/P[0][i]*Vfrozen[0][i];
                Vfrozen[2][i]=Vfrozen[2][i]+P[2][i]/P[0][i]*Vfrozen[0][i];
                Vfrozen[3][i]=0.0;

                V[1][i]=Vfrozen[1][i];
                V[2][i]=Vfrozen[2][i];
                V[3][i]=Vfrozen[3][i];
                       
                V[0][i]=-log(1.0-ran0(&NUM1));
  
                // adjust momentum to fit energy and mass

                if(abs(KATT1[i])==1||abs(KATT1[i])==2||abs(KATT1[i])==3||abs(KATT1[i])==21)
                {
                P[4][i]=0.0;
                P[0][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]+P[3][i]*P[3][i]+P[4][i]*P[4][i]);
                P[5][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]);
                WT[i]=1.0;
                }


                if(abs(KATT1[i])==4)
                {
                P[4][i]=1.27;
                P[0][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]+P[3][i]*P[3][i]+P[4][i]*P[4][i]);
                P[5][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]);
                WT[i]=1.0;
                }


                if(abs(KATT1[i])==5)
                {
                P[4][i]=4.19;
                P[0][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]+P[3][i]*P[3][i]+P[4][i]*P[4][i]);
                P[5][i]=sqrt(P[1][i]*P[1][i]+P[2][i]*P[2][i]);
                WT[i]=1.0;
                }

                fake0[i]=1; // SC: no matching negative partons for initial jet partons

                EiTot=EiTot+P[0][i];

                for(int j=0;j<=3;j++) Prad[j][i]=P[j][i];
  
                Tfrozen[i]=0.0;
                vcfrozen[1][i]=0.0;
                vcfrozen[2][i]=0.0;
                vcfrozen[3][i]=0.0;
            }

//            cout<<"EiTot: "<<EiTot<<endl;

            // reset position information if necessary
            if(flagJetX==1) setJetX(numInitXY);

        } 

        np=nj;
        
//...end initilization of jet parton

//...time evolution in LBT if in medium

        if(vacORmed==1) {
		
            for(double ti=time0+dt;ti<=timend+epsilon;ti=ti+dt) {

                LBT(n,ti);

                ffqperp << ti << "    " << nqperp << endl;
                for(int kk=0; kk<nqperp; kk++)
                {
                        ffqperp << qperp[kk] << endl;
                        qperp[kk]=0.0;
                }
                nqperp=0;

            }
        }

//...end of time evolution in LBT

        numEvent=numEvent+1;		

//...write outputs

        double EfTot=0.0; // check energy conservation

        int ip=0;
        int iHQ=0;
        for(int i=1;i<=np;i++) {
//            if(P[0][i] < cutOut && abs(KATT1[i])!=4 && abs(KATT1[i])!=5) ip+=1;
            if(abs(KATT1[i])==4 || abs(KATT1[i])==5) iHQ+=1;
            EfTot=EfTot+P[0][i];
        }	
        int nnn=np-ip;
        int ip0=0;
        for(int i=2;i<=np;i++) {
            if(P0[0][i] < cutOut) ip0+=1;
            EfTot=EfTot-P0[0][i];
        }	
        int nnn0=np-ip0;

//        cout<<"EfTot: "<<EfTot<<endl;
        if (etaflag == 0) {
        if(outFormat==1){

            if(heavyOut==1) outHQ << n << "    " << iHQ << endl;
            if(lightOut==1) positive << n << "    " << nnn << "    " << nj << "    " << num1 << endl;
            for(int i=1;i<=np;i++) { // if only want to write out leading partons, change np->nj

//                if(P[0][i] < cutOut && abs(KATT1[i])!=4) continue; // throw away light partons below cut
                if(abs(KATT1[i])==4 || abs(KATT1[i])==5){ // write out heavy quark
                    if(heavyOut==1 && lightOut==1) {
                        outHQ << i << "    "  << KATT1[i] << "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << P[0][i] << "    " << P[4][i] << "    " << Vfrozen[1][i] << "    " << Vfrozen[2][i] << "    " << Vfrozen[3][i] << "    " << Vfrozen[0][i] << "    " << color[i] << "    " << anticolor[i]  << endl;
                        positive << i << "    "  << KATT1[i] << "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << P[0][i] << "    " << P[4][i] << "    " << Vfrozen[1][i] << "    " << Vfrozen[2][i] << "    " << Vfrozen[3][i] << "    " << Vfrozen[0][i] << "    " << color[i] << "    " << anticolor[i]  << endl;
                    }
                    if(heavyOut==1 && lightOut!=1) outHQ << i << "    "  << KATT1[i] << "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << P[0][i] << "    " << Vfrozen[1][i] << "    " << Vfrozen[2][i] << "    " << Vfrozen[3][i] << "    " << Vfrozen[0][i] << endl;
                } else { // write out light parton
                    if(lightOut==1) positive << i << "    " << KATT1[i] <<  "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << P[0][i] << "    " << P[4][i] << "    " << Vfrozen[1][i] << "    " << Vfrozen[2][i] << "    " << Vfrozen[3][i] << "    " << Vfrozen[0][i] << "    " << color[i] << "    " << anticolor[i]  << endl;
                }
            }	



            if(lightOut==1) negative << n << "    " << nnn0-1 << endl;
            for(int i=2;i<=np;i++) {
                if(fake0[i]==1) continue; // SC: only throw away paritcles that don't exist -- negative that match emitted gluons
//                if(P0[0][i] < cutOut) continue;
                if(lightOut==1) negative << i << "   " << KATT10[i] << "    " << P0[1][i] << "    " << P0[2][i] << "    " << P0[3][i] << "    " << P0[0][i] << "    " << P0[4][i] << "    " << vcfrozen0[1][i] << "    " << vcfrozen0[2][i] << "    " << vcfrozen0[3][i] << "    " << Vfrozen0[0][i] << "    " << color0[i] << "    " << anticolor0[i]   << endl;
            }	

        } else if(outFormat==2) { // for JETSCAPE patch code

            if(heavyOut==1) outHQ << n << "    " << iHQ << "    " << initPT << "    " << rangePT << endl;
            if(lightOut==1) positive << n << "    " << nnn << "    " << initPT << "    " << rangePT << endl;
            for(int i=1;i<=np;i++) { // if only want to write out leading partons, change np->nj
                if(P[0][i] < cutOut && abs(KATT1[i])!=4) continue; // throw away light partons below cut
                if(abs(KATT1[i])==4){ // write out heavy quark
                    if(heavyOut==1) outHQ << i << "    " << KATT1[i] <<  "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << P[0][i] << "    " << P[4][i] << "    " << Vfrozen[1][i] << "    " << Vfrozen[2][i] << "    " << Vfrozen[3][i] << "    " << Vfrozen[0][i] << endl;
                } else { // write out light parton
                    if(lightOut==1) positive << i << "    " << KATT1[i] <<  "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << P[0][i] << "    " << P[4][i] << "    " << Vfrozen[1][i] << "    " << Vfrozen[2][i] << "    " << Vfrozen[3][i] << "    " << Vfrozen[0][i] << endl;
                }
            }	

            if(lightOut==1) negative << n << "    " << nnn0-1 << "    " << initPT << "    " << rangePT << endl;
            for(int i=2;i<=np;i++) {
                if(P0[0][i] < cutOut) continue;
                if(lightOut==1) negative << i << "   " << KATT10[i] << "    " << P0[1][i] << "    " << P0[2][i] << "    " << P0[3][i] << "    " << P0[0][i] << "    " << P0[4][i] << "    " << Vfrozen0[1][i] << "    " << Vfrozen0[2][i] << "    " << Vfrozen0[3][i] << "    " << Vfrozen0[0][i] << endl;
            }	

        }
        } else {

         positive << n << "    " << nnn << "    " << nj << "    " << num1 << endl;     

        }

//...end of output

        int print=n%nprint;

	if(print==0) {
            cout << "  " << endl;
            cout<< "n" << "  " << "np" << endl;      
            cout << n << "  " << np << endl;
	}


    } //for(int n=1;n<=ncall;n++)

//...end of event loop (e.g. ncall)

    if(lightOut==1){ 
        positive.close();
        negative.close();
    }
    if(heavyOut==1){
        outHQ.close();
    }
    ffqperp.close(); 

//...check system time, computational cost, etc.

    struct tm *local_end;
    time_t time_end;
    time_end=time(NULL);
    local_end=localtime(&time_end);
    
    char buf2[80];
    strftime(buf2,80,"Current Time: %Y-%m-%d %H:%M:%S",local_end);
    cout << "the program ends at:" << endl;
    cout << buf2 << endl;
    
    unsigned cost,nh,nm,ns;
    cost=difftime(time_end,time_start);
    
    nh=cost/3600;
    nm=(cost%3600)/60;
    ns=(cost%3600)%60;
    
    cout << "the program costs:" << endl;
    cout << cost << "s:" << " " << nh << "h" << " " << nm << "m" << " " << ns << "s" << endl;	
    cout << "counth100" << " " << counth100 <<endl;
    
}



