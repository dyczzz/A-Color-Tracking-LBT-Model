///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Define functions for LBT
//
///////////////////////////////////////////////////////////////////////////////////////////////////

void read_tables(){ // intialize various tables for LBT

     if(bulkFlag==1) { // read in OSU hydro profiles
         int dataID_in=1;
         char dataFN_in[]="hydroProfile/JetData.dat";
         int ctlID_in=2;
         char ctlFN_in[]="hydroProfile/JetCtl.dat";
         int bufferSize=1000;
         int len1=strlen(dataFN_in);
         int len2=strlen(ctlFN_in);
         sethydrofilesez_(&dataID_in, dataFN_in, &ctlID_in, ctlFN_in, &bufferSize, len1, len2);
     } else if(bulkFlag==2) { // read in CCNU hydro profiles
         char dataFN_in[]="hydroProfile/bulk.dat";
         int len1=strlen(dataFN_in);
         read_ccnu_(dataFN_in,len1);
     }
 
     if(fixPosition!=1) read_xyMC(numInitXY);

//...read scattering rate
     int it,ie;
     int n=500;
     ifstream f1("tables/ratedata");
     if(!f1.is_open())
     {
      cout<<"Erro openning date file1!\n";
     }
     else
     {
      for(int i=1;i<=n;i++)
	  {
	  f1>>it>>ie;
	  f1>>qhatG[it][ie]>>Rg[it][ie]>>Rg1[it][ie]>>Rg2[it][ie]>>Rg3[it][ie]>>qhatLQ[it][ie]>>Rq[it][ie]>>Rq3[it][ie]>>Rq4[it][ie]>>Rq5[it][ie]>>Rq6[it][ie]>>Rq7[it][ie]>>Rq8[it][ie];	      
	  }
     }
     f1.close();

// duplicate for heavy quark
     ifstream f11("tables/ratedata-HQ");
//     ifstream f11("b-tables/ratedata-HQ");
     if(!f11.is_open())
     {
      cout<<"Erro openning HQ data file!\n";
     }
     else
     {
      for(int i=1;i<=n;i++)
	  {
	  f11>>it>>ie;
	  f11>>RHQ[it][ie]>>RHQ11[it][ie]>>RHQ12[it][ie]>>qhatHQ[it][ie];	      
	  }
     }
     f11.close();

// read radiation table for heavy quark
     ifstream f12("tables/dNg_over_dt_cD6.dat");
//     ifstream f12("b-tables/dNg_over_dt_bD6.dat");
     ifstream f13("tables/dNg_over_dt_qD6.dat");
     ifstream f14("tables/dNg_over_dt_gD6.dat");
     if(!f12.is_open()||!f13.is_open()||!f14.is_open())
     {
      cout<<"Erro openning HQ radiation table file!\n";
     }
     else
     {
     for(int k=1; k<=t_gn; k++) {
         char dummyChar[100];
         long double dummyD;
         f12 >> dummyChar >> dummyChar >> dummyChar >> dummyChar;
         f13 >> dummyChar >> dummyChar >> dummyChar >> dummyChar;
         f14 >> dummyChar >> dummyChar >> dummyChar >> dummyChar;
         for(int i=1; i<=temp_gn; i++) {
            dNg_over_dt_c[k+1][i][0]=0.0;
            dNg_over_dt_q[k+1][i][0]=0.0;
            dNg_over_dt_g[k+1][i][0]=0.0;
            max_dNgfnc_c[k+1][i][0]=0.0;
            max_dNgfnc_q[k+1][i][0]=0.0;
            max_dNgfnc_g[k+1][i][0]=0.0;
            for(int j=1; j<=HQener_gn; j++) {
                f12 >> dNg_over_dt_c[k+1][i][j] >> max_dNgfnc_c[k+1][i][j];
                f13 >> dNg_over_dt_q[k+1][i][j] >> max_dNgfnc_q[k+1][i][j];
                f14 >> dNg_over_dt_g[k+1][i][j] >> max_dNgfnc_g[k+1][i][j];
            }
         }
     }
     }
//     cout << dNg_over_dt_c[t_gn+1][temp_gn][HQener_gn] << "    " << max_dNgfnc_c[t_gn+1][temp_gn][HQener_gn] << endl;
//     cout << dNg_over_dt_q[t_gn+1][temp_gn][HQener_gn] << "    " << max_dNgfnc_q[t_gn+1][temp_gn][HQener_gn] << endl;
//     cout << dNg_over_dt_g[t_gn+1][temp_gn][HQener_gn] << "    " << max_dNgfnc_g[t_gn+1][temp_gn][HQener_gn] << endl;
     f12.close();
     f13.close();
     f14.close();

     for(int i=1; i<=temp_gn; i++) {
         for(int j=1; j<=HQener_gn; j++) {
             dNg_over_dt_c[1][i][j]=0.0;
             dNg_over_dt_q[1][i][j]=0.0;
             dNg_over_dt_g[1][i][j]=0.0;
             max_dNgfnc_c[1][i][j]=0.0;
             max_dNgfnc_q[1][i][j]=0.0;
             max_dNgfnc_g[1][i][j]=0.0;
         }
     }

// preparation for HQ 2->2
     ifstream fileB("tables/distB.dat");
//     ifstream fileB("b-tables/distB.dat");
     if(!fileB.is_open()) {
        cout << "Erro openning data file distB.dat!" << endl;
     } else {
        for(int i=0;i<N_T;i++) {
           for(int j=0;j<N_p1;j++) {
              double dummy_T,dummy_p1;
              fileB>>dummy_T>>dummy_p1;
              if(fabs(min_T+(0.5+i)*bin_T-dummy_T)>1.0e-5 || fabs(min_p1+(0.5+j)*bin_p1-dummy_p1)>1.0e-5) {
                  cout << "Erro in reading data file distB.dat!" << endl;
                  exit (EXIT_FAILURE);
              }
              fileB>>distFncBM[i][j];
              for(int k=0;k<N_e2;k++) fileB>>distFncB[i][j][k];
              for(int k=0;k<N_e2;k++) fileB>>distMaxB[i][j][k];
           }
        }
     }
     fileB.close();
   
     ifstream fileF("tables/distF.dat");
//     ifstream fileF("b-tables/distF.dat");
     if(!fileF.is_open()) {
        cout << "Erro openning data file distF.dat!" << endl;
     } else {
        for(int i=0;i<N_T;i++) {
           for(int j=0;j<N_p1;j++) {
              double dummy_T,dummy_p1;
              fileF>>dummy_T>>dummy_p1;
              if(fabs(min_T+(0.5+i)*bin_T-dummy_T)>1.0e-5 || fabs(min_p1+(0.5+j)*bin_p1-dummy_p1)>1.0e-5) {
                  cout << "Erro in reading data file distF.dat!" << endl;
                  exit (EXIT_FAILURE);
              }
              fileF>>distFncFM[i][j];
              for(int k=0;k<N_e2;k++) fileF>>distFncF[i][j][k];
              for(int k=0;k<N_e2;k++) fileF>>distMaxF[i][j][k];
           }
        }
     }
     fileF.close();
}

//..............................................................subroutine
//..............................................................

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//#################
#include <time.h>
//###############

float ran0(long *idum)

{
  int j;
  long k;
  //##############
  static long seed_w=(unsigned) time(NULL);
 
  static long idum2=seed_w;
  //##############
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) { 
	if (-(*idum) < 1) *idum=1; 
	else *idum = -(*idum);
	for (j=NTAB+7;j>=0;j--) { 
	  k=(*idum)/IQ1;
	  *idum=IA1*(*idum-k*IQ1)-k*IR1;
	  if (*idum < 0) *idum += IM1;
	  if (j < NTAB) iv[j] = *idum;
	}
	iy=iv[0];
  }
  k=(*idum)/IQ1; 
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1; 
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 
  iy=iv[j]-idum2; 
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}
	
//.........................................................................
double alphas0(int &Kalphas,double temp0){
      
  double X;
  if(Kalphas==1) 
	{      
//      X=fixAlphas;      
        X=0.3;
	}
  return X;
}
//.........................................................................
double DebyeMass2(int &Kqhat0,double alphas,double temp0){

  double Y;
  if(Kqhat0==1)
	{      
      Y=4.0*pi*alphas*pow(temp0,2);      
	}
  if(Kqhat0==2)
	{      
      Y=(3.0/2.0)*4.0*pi*alphas*pow(temp0,2);      
	}
  if(Kqhat0==3)
	{      
      Y=1.0;      
	}
  return Y;

}
//.........................................................................

//.........................................................................	
void titau(double ti,double vf[4],double vp[4],double p0[4],double &Vx,double &Vy,double &Veta,double &Xtau){   	

  //..............................................................test part
  //		  cout<<"ti"<<" "<<ti<<" "<<"vf"<<" "<<vf[1]<<" "<<vf[2]<<" "<<vf[3]<<endl;
  //		  cout<<"vp[4]"<<" "<<vp[1]<<" "<<vp[2]<<" "<<vp[3]<<" "<<vp[0]<<endl;	
  //		  cout<<"p0[4]"<<" "<<p0[1]<<" "<<p0[2]<<" "<<p0[3]<<" "<<p0[0]<<endl;		  
  //..............................................................test part

  //....notice the form of vf
  double gamma=1.0/sqrt(1-(vf[1]*vf[1]+vf[2]*vf[2]));
  double mt=sqrt(p0[1]*p0[1]+p0[2]*p0[2]);
  double Yp=1.0/2.0*log((p0[0]+p0[3])/(p0[0]-p0[3]));	
  double etas=vp[3];
  double etaf=atanh(vf[3])+etas;
  double pper=sqrt(p0[1]*p0[1]+p0[2]*p0[2]);
  double vper=sqrt(vf[1]*vf[1]+vf[2]*vf[2]);
  double pvper=p0[1]*vf[1]+p0[2]*vf[2];

  Vx=p0[1]/pper/cosh(Yp-etas);	  
  Vy=p0[2]/pper/cosh(Yp-etas);
  Veta=(1.0/ti)*tanh(Yp-etas);
	  
  Xtau=(gamma*mt*cosh(Yp-etaf)-pvper*vper)/(mt*cosh(Yp-etas));

  //..............................................................test part
  //		  cout<<"gamma"<<" "<<gamma<<" "<<"mt"<<" "<<mt<<" "<<"Yp"<<" "<<Yp<<endl;
  //		  cout<<"etas"<<" "<<etas<<" "<<"etaf"<<" "<<etaf<<" "<<"pper"<<" "<<pper<<endl;	
  //		  cout<<"vper"<<" "<<vper<<" "<<"pvper"<<" "<<pvper<<endl;
  //		  cout<<"Vx"<<" "<<Vx<<" "<<"Vy"<<" "<<Vy<<" "<<"Veta"<<" "<<Veta<<endl;
  //		  cout<<"Xtau"<<" "<<Xtau<<endl;		  
  //..............................................................test part	  
}
//.........................................................................	
	
//.........................................................................
	
void lam(int KATT0,double &RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2){   
	  
  double dtemp=0.02;
  iT1=(int)((T-0.1)/dtemp);
  iT2=iT1+1;
  iE1=(int)(log(E)+2);
  if(iE1<1) iE1=1;
  iE2=iE1+1;

  T1=0.12+(iT1-1)*0.02;
  T2=T1+dtemp;
  E1=exp(iE1-2.0);
  E2=exp(iE2-2.0);
	  
  if(KATT0==21) {	
	  double RTE1=(Rg[iT2][iE1]-Rg[iT1][iE1])*(T-T1)/(T2-T1)+Rg[iT1][iE1];
	  double RTE2=(Rg[iT2][iE2]-Rg[iT1][iE2])*(T-T1)/(T2-T1)+Rg[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
  } else if (KATT0==4||KATT0==-4) { // add heavy quark channel
	  double RTE1=(RHQ[iT2][iE1]-RHQ[iT1][iE1])*(T-T1)/(T2-T1)+RHQ[iT1][iE1];
	  double RTE2=(RHQ[iT2][iE2]-RHQ[iT1][iE2])*(T-T1)/(T2-T1)+RHQ[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	   
  } else {    
	  double RTE1=(Rq[iT2][iE1]-Rq[iT1][iE1])*(T-T1)/(T2-T1)+Rq[iT1][iE1];
	  double RTE2=(Rq[iT2][iE2]-Rq[iT1][iE2])*(T-T1)/(T2-T1)+Rq[iT1][iE2];
	  RTE=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
//          cout<<"RTE2,RTE1,E,E1,E2,RTE: "<<RTE2<<"  "<<RTE1<<"  "<<E<<"  "<<E1<<"  "<<E2<<"  "<<RTE<<endl;
  }

}

//.........................................................................

void flavor(int &CT,int &KATT0,int &KATT2,int &KATT3,double RTE,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2){  

  double RTEg;
  double RTEg1;
  double RTEg2;
  double RTEg3;
  double RTEq;
  double RTEq3;
  double RTEq4;
  double RTEq5;
  double RTEq6;
  double RTEq7;
  double RTEq8;
  double RTEHQ;
  double RTEHQ11;
  double RTEHQ12;
  double qhatTP;

  int vb[7]={0};
  int b=0;
  int KATT00=KATT0;


  vb[1]=1; 
  vb[2]=2; 
  vb[3]=3; 
  vb[4]=-1; 
  vb[5]=-2; 
  vb[6]=-3;

	  
  //.....
  linear(KATT0,E,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2,RTEg,RTEg1,RTEg2,RTEg3,RTEq,RTEq3,RTEq4,RTEq5,RTEq6,RTEq7,RTEq8,RTEHQ,RTEHQ11,RTEHQ12,qhatTP);

  if(KATT00==21) {   //.....for gluon
	  double R0=RTE;
	  double R1=RTEg1;
	  double R2=RTEg2;
	  double R3=RTEg3;
	   	
	  double a=ran0(&NUM1);
	   
	  if(a<=R1/R0) {
	      CT=1;
	      KATT3=21;
	      KATT2=21;
	      KATT0=21;
	  }	  
	  
	  if(a>R1/R0 && a<=(R1+R2)/R0) {
	      CT=2;
	      b=floor(ran0(&NUM1)*6+1);
	      if(b==7) {
		  b=6;
	      }
	      KATT3=21;
	      KATT2=vb[b];
	      KATT0=-KATT2;
	  }	  
	  
	  if(a>(R1+R2)/R0 && a<=1.0) {
	      CT=3;
	      b=floor(ran0(&NUM1)*6+1);
	      if(b==7) {
	          b=6;
	      }
	      KATT3=vb[b];
	      KATT2=KATT3;
	      KATT0=21;
	  }
  } else if(KATT00==4||KATT00==-4) { // for heavy quark
          double R0=RTE;
	  double R1=RTEHQ11;
	  double R2=RTEHQ12;
	   	
	  double a=ran0(&NUM1);

//          qhat_over_T3=qhatTP;  // what is read in is qhat/T^3 of quark
//          D2piT=8.0*pi/qhat_over_T3;

          if(a<=R1/R0) { //Qq->Qq
              CT=11;
	      b=floor(ran0(&NUM1)*6+1);
	      if(b==7) {
	          b=6;
	      }
	      KATT3=vb[b];
              KATT2=KATT3;
          } else { //Qg->Qg
	      CT=12;
	      KATT3=21;
	      KATT2=KATT3;
	  }	  
	  
  } else { //.....for quark and antiquark (light)
	  double R00=RTE;
	  double R3=RTEq3;
	  double R4=RTEq4;
	  double R5=RTEq5;
	  double R6=RTEq6;
	  double R7=RTEq7;
	  double R8=RTEq8;

	  double a=ran0(&NUM1);
	  if(a<=R3/R00)
	    { 
	      CT=13;
	      KATT3=21;
          KATT2=21; 
	      KATT0=KATT0;
	    }
	  	  
	  if(a>R3/R00 && a<=(R3+R4)/R00)
		{ 
	      CT=4;
		f1:	      b=floor(ran0(&NUM1)*6+1);
		  if(b==7)
		    {
			  b=6;
	        }
		  KATT3=vb[b]; 
		  if(KATT3==KATT0)
			{
			  goto f1;
	        }
		  KATT2=KATT3;
		  //	      KATT0=KATT0
	    }
	  	  
	  if(a>(R3+R4)/R00 && a<=(R3+R4+R5)/R00)
	    { 
	      CT=5;
	      KATT3=KATT0;
	      KATT2=KATT0;
	    }
	  //.....the only difference between quark and antiquark	   	   
	    
	  if(a>(R3+R4+R5)/R00 && a<=(R3+R4+R5+R6)/R00)
		{
          CT=6;
	      KATT3=-KATT0;  
		f2:	      b=floor(ran0(&NUM1)*3+1);
		  if(b==4)
		    {
			  b=3;
	        }
		  KATT2=-KATT0/abs(KATT0)*vb[b]; 
		  if(abs(KATT2)==abs(KATT3))
			{
			  goto f2;
            }
		  KATT0=-KATT2;
	    }
	   	   
	  if(a>(R3+R4+R5+R6)/R00 && a<=(R3+R4+R5+R6+R7)/R00)
        {	   
	      CT=7;
	      KATT3=-KATT0; 
	      KATT2=KATT3;
		  //	      KATT0=KATT0
	    }	   
	   
	  if(a>(R3+R4+R5+R6+R7)/R00 && a<=1.0)
	    {
		  CT=8;
		  KATT3=-KATT0; 
		  KATT2=21; 
		  KATT0=21;
	    }	   	  
   }	
}
 	
  

//.........................................................................

void linear(int KATT,double E,double T,double &T1,double &T2,double &E1,double &E2,int &iT1,int &iT2,int &iE1,int &iE2,double &RTEg,double &RTEg1,double &RTEg2,double &RTEg3,double &RTEq,double &RTEq3,double &RTEq4,double &RTEq5,double &RTEq6,double &RTEq7,double &RTEq8,double &RTEHQ,double &RTEHQ11,double &RTEHQ12,double &qhatTP){
  if(KATT==21){
	  double RTE1=(Rg1[iT2][iE1]-Rg1[iT1][iE1])*(T-T1)/(T2-T1)+Rg1[iT1][iE1];
	  double RTE2=(Rg1[iT2][iE2]-Rg1[iT1][iE2])*(T-T1)/(T2-T1)+Rg1[iT1][iE2];
	  RTEg1=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
	
	  RTE1=(Rg2[iT2][iE1]-Rg2[iT1][iE1])*(T-T1)/(T2-T1)+Rg2[iT1][iE1];
	  RTE2=(Rg2[iT2][iE2]-Rg2[iT1][iE2])*(T-T1)/(T2-T1)+Rg2[iT1][iE2];
	  RTEg2=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	
	
	  RTE1=(Rg3[iT2][iE1]-Rg3[iT1][iE1])*(T-T1)/(T2-T1)+Rg3[iT1][iE1];
	  RTE2=(Rg3[iT2][iE2]-Rg3[iT1][iE2])*(T-T1)/(T2-T1)+Rg3[iT1][iE2];
	  RTEg3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;		

//	  RTE1=(qhatG[iT2][iE1]-qhatG[iT1][iE1])*(T-T1)/(T2-T1)+qhatG[iT1][iE1];
//	  RTE2=(qhatG[iT2][iE2]-qhatG[iT1][iE2])*(T-T1)/(T2-T1)+qhatG[iT1][iE2];
//	  qhatTP=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	

  } else if(KATT==4||KATT==-4){  // add heavy quark channel
	  double RTE1=(RHQ11[iT2][iE1]-RHQ11[iT1][iE1])*(T-T1)/(T2-T1)+RHQ11[iT1][iE1];
	  double RTE2=(RHQ11[iT2][iE2]-RHQ11[iT1][iE2])*(T-T1)/(T2-T1)+RHQ11[iT1][iE2];
	  RTEHQ11=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
	
	  RTE1=(RHQ12[iT2][iE1]-RHQ12[iT1][iE1])*(T-T1)/(T2-T1)+RHQ12[iT1][iE1];
	  RTE2=(RHQ12[iT2][iE2]-RHQ12[iT1][iE2])*(T-T1)/(T2-T1)+RHQ12[iT1][iE2];
	  RTEHQ12=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	

//	  RTE1=(qhatHQ[iT2][iE1]-qhatHQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE1];
//	  RTE2=(qhatHQ[iT2][iE2]-qhatHQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE2];
//	  qhatTP=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	

  } else {
	  double RTE1=(Rq3[iT2][iE1]-Rq3[iT1][iE1])*(T-T1)/(T2-T1)+Rq3[iT1][iE1];
	  double RTE2=(Rq3[iT2][iE2]-Rq3[iT1][iE2])*(T-T1)/(T2-T1)+Rq3[iT1][iE2];
	  RTEq3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq4[iT2][iE1]-Rq4[iT1][iE1])*(T-T1)/(T2-T1)+Rq4[iT1][iE1];
	  RTE2=(Rq4[iT2][iE2]-Rq4[iT1][iE2])*(T-T1)/(T2-T1)+Rq4[iT1][iE2];
	  RTEq4=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq5[iT2][iE1]-Rq5[iT1][iE1])*(T-T1)/(T2-T1)+Rq5[iT1][iE1];
	  RTE2=(Rq5[iT2][iE2]-Rq5[iT1][iE2])*(T-T1)/(T2-T1)+Rq5[iT1][iE2];
	  RTEq5=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq6[iT2][iE1]-Rq6[iT1][iE1])*(T-T1)/(T2-T1)+Rq6[iT1][iE1];
	  RTE2=(Rq6[iT2][iE2]-Rq6[iT1][iE2])*(T-T1)/(T2-T1)+Rq6[iT1][iE2];
	  RTEq6=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq7[iT2][iE1]-Rq7[iT1][iE1])*(T-T1)/(T2-T1)+Rq7[iT1][iE1];
	  RTE2=(Rq7[iT2][iE2]-Rq7[iT1][iE2])*(T-T1)/(T2-T1)+Rq7[iT1][iE2];
	  RTEq7=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq8[iT2][iE1]-Rq8[iT1][iE1])*(T-T1)/(T2-T1)+Rq8[iT1][iE1];
	  RTE2=(Rq8[iT2][iE2]-Rq8[iT1][iE2])*(T-T1)/(T2-T1)+Rq8[iT1][iE2];
	  RTEq8=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

//	  RTE1=(qhatLQ[iT2][iE1]-qhatLQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE1];
//	  RTE2=(qhatLQ[iT2][iE2]-qhatLQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE2];
//	  qhatTP=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	

  }	
}

//.........................................................................

void twflavor(int &CT,int &KATT0,int &KATT2,double E,double T){  

  double RTEg;
  double RTEg1;
  double RTEg2;
  double RTEg3;
  double RTEq;
  double RTEq3;
  double RTEq4;
  double RTEq5;
  double RTEq6;
  double RTEq7;
  double RTEq8;

  int vb[7]={0};
  int b=0;
  int KATT00=KATT0;
  int KATT20=KATT2;
	  
  vb[1]=1;
  vb[2]=2; 
  vb[3]=3; 
  vb[4]=-1; 
  vb[5]=-2; 
  vb[6]=-3; 

  twlinear(KATT0,E,T,RTEg1,RTEg2,RTEq6,RTEq7,RTEq8);

  //.....for gluon
  if(KATT00==21)
	{
	  //	R0  =RTE
	  double R1  =RTEg1;
	  double R2  =RTEg2;
	  //	R3  =RTEg3

	  if(KATT20==21)
		{
		  double a=ran0(&NUM1);
		  if(a<=R1/(R1+R2))
			{
			  CT=1;
			  //	        KATT3=KATT2
			  //		    KATT2=21 
			  //		    KATT0=21
			}
	     
		  if(a>R1/(R1+R2))
			{
			  CT=2;
			  //	        KATT3=KATT2 
			  b=floor(ran0(&NUM1)*6+1);
			  if(b==7)
				{
				  b=6;
				}
			  KATT2=vb[b]; 
			  KATT0=-KATT2;
			}
		}	  
         
	  if(KATT20!=21)
		{ 
	      CT=3;
		  //	      KATT3=KATT2 
		  //		  KATT2=KATT2 
		  //		  KATT0=21
		}	     
	}   	            
	  
  //.....for quark and antiquark
  if(KATT00!=21)
	{      

	  //	R00 =RTE
	  //	R3  =RTEq3
	  //	R4  =RTEq4
	  //	R5  =RTEq5
	  double R6  =RTEq6;
	  double R7  =RTEq7;
	  double R8  =RTEq8;
	  double R00 =R6+R7+R8;

	  if(KATT20==21)
		{
	      CT=13;
		  //	      KATT3=KATT2 
		  //		  KATT2=21 
		  //		  KATT0=KATT0
		}
	  	  
	  if(KATT20!=21)
	    {

	      if(abs(KATT20)!=abs(KATT00))
			{
			  CT=4;
			  //	      KATT3=KATT2 
			  //		  KATT2=KATT3 
			  //		  KATT0=KATT0
			}
	  	  
	      if(KATT20==KATT00)
			{
			  CT=5;
			  //	      KATT3=KATT2
			  //	      KATT0=KATT0
			}
   	      
	      if(KATT20==-KATT00)
			{
			  double a=ran0(&NUM1);
			  if(a<=(R6)/R00)
				{
				  CT=6;
				  //	         KATT3=KATT2
				tf2:	     b=floor(ran0(&NUM1)*3+1);
				  if(b==4)
					{
					  b=3;
					}  
				  KATT2=-KATT0/abs(KATT0)*vb[b]; 
				  if(abs(KATT2)==abs(KATT0))
					{
					  goto tf2;			
					} 
				  KATT0=-KATT2;	     
				} 
	   	   
			  if(a>(R6)/R00 && a<=(R6+R7)/R00)
				{
				  CT=7;
				  //	         KATT3=KATT2 
				  //		     KATT2=KATT3 
				  //		     KATT0=KATT0
				}	   
	   
			  if(a>(R6+R7)/R00 && a<=1.0)
				{
				  CT=8;
				  //	         KATT3=KATT2 
				  KATT2=21; 
				  KATT0=21;
				}	      
			}	   	  	   
	    } 
    }      	

}  

//.........................................................................

void twlinear(int KATT,double E,double T,double &RTEg1,double &RTEg2,double &RTEq6,double &RTEq7,double &RTEq8){ 

  //.....    
  double dtemp=0.02;
  int iT1=floor((T-0.1)/dtemp);
  int iT2=iT1+1;
  int iE1=floor(log(E)+2);
  int iE2=iE1+1;
  //
  double T1=0.12+(iT1-1)*0.02;
  double T2=T1+dtemp;
  double E1=exp(iE1-2.0);
  double E2=exp(iE2-2.0);
  //
  if(KATT==21)
	{
	  double RTE1=(Rg1[iT2][iE1]-Rg1[iT1][iE1])*(T-T1)/(T2-T1)+Rg1[iT1][iE1];
	  double RTE2=(Rg1[iT2][iE2]-Rg1[iT1][iE2])*(T-T1)/(T2-T1)+Rg1[iT1][iE2];
	  RTEg1=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;
	
	  RTE1=(Rg2[iT2][iE1]-Rg2[iT1][iE1])*(T-T1)/(T2-T1)+Rg2[iT1][iE1];
	  RTE2=(Rg2[iT2][iE2]-Rg2[iT1][iE2])*(T-T1)/(T2-T1)+Rg2[iT1][iE2];
	  RTEg2=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;	
	
	  //	   RTE1=(Rg3[iT2][iE1]-Rg3[iT1][iE1])*(T-T1)/(T2-T1)+Rg3[iT1][iE1];
	  //	   RTE2=(Rg3[iT2][iE2]-Rg3[iT1][iE2])*(T-T1)/(T2-T1)+Rg3[iT1][iE2];
	  //	   RTEg3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;		
     
	}

  if(KATT!=21)
	{
	  //	   RTE1=(Rq3[iT2][iE1]-Rq3[iT1][iE1])*(T-T1)/(T2-T1)+Rq3[iT1][iE1];
	  //	   RTE2=(Rq3[iT2][iE2]-Rq3[iT1][iE2])*(T-T1)/(T2-T1)+Rq3[iT1][iE2];
	  //	   RTEq3=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  //	   RTE1=(Rq4[iT2][iE1]-Rq4[iT1][iE1])*(T-T1)/(T2-T1)+Rq4[iT1][iE1];
	  //	   RTE2=(Rq4[iT2][iE2]-Rq4[iT1][iE2])*(T-T1)/(T2-T1)+Rq4[iT1][iE2];
	  //	   RTEq4=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1

	  //	   RTE1=(Rq5[iT2][iE1]-Rq5[iT1][iE1])*(T-T1)/(T2-T1)+Rq5[iT1][iE1];
	  //	   RTE2=(Rq5[iT2][iE2]-Rq5[iT1][iE2])*(T-T1)/(T2-T1)+Rq5[iT1][iE2];
	  //	   RTEq5=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  double RTE1=(Rq6[iT2][iE1]-Rq6[iT1][iE1])*(T-T1)/(T2-T1)+Rq6[iT1][iE1];
	  double RTE2=(Rq6[iT2][iE2]-Rq6[iT1][iE2])*(T-T1)/(T2-T1)+Rq6[iT1][iE2];
	  RTEq6=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq7[iT2][iE1]-Rq7[iT1][iE1])*(T-T1)/(T2-T1)+Rq7[iT1][iE1];
	  RTE2=(Rq7[iT2][iE2]-Rq7[iT1][iE2])*(T-T1)/(T2-T1)+Rq7[iT1][iE2];
	  RTEq7=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	  RTE1=(Rq8[iT2][iE1]-Rq8[iT1][iE1])*(T-T1)/(T2-T1)+Rq8[iT1][iE1];
	  RTE2=(Rq8[iT2][iE2]-Rq8[iT1][iE2])*(T-T1)/(T2-T1)+Rq8[iT1][iE2];
	  RTEq8=(RTE2-RTE1)*(E-E1)/(E2-E1)+RTE1;

	}	
}	

//.........................................................................

void trans(double v[4],double p[4]){	
  double vv=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);	
  double ga=1.0/sqrt(1.0-vv*vv);
  double ppar=p[1]*v[1]+p[2]*v[2]+p[3]*v[3];
  double gavv=(ppar*ga/(1.0+ga)-p[0])*ga;
  p[0]=ga*(p[0]-ppar);
  p[1]=p[1]+v[1]*gavv;
  p[2]=p[2]+v[2]*gavv;
  p[3]=p[3]+v[3]*gavv;
}

void transback(double v[4],double p[4]){
  double vv=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  double ga=1.0/sqrt(1.0-vv*vv); 
  double ppar=p[1]*v[1]+p[2]*v[2]+p[3]*v[3];
  double gavv=(-ppar*ga/(1.0+ga)-p[0])*ga;
  p[0]=ga*(p[0]+ppar);
  p[1]=p[1]-v[1]*gavv;
  p[2]=p[2]-v[2]*gavv;
  p[3]=p[3]-v[3]*gavv;
}	

//.........................................................................
void colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt){
  //
  //    p0 initial jet momentum, output to final momentum
  //    p2 final thermal momentum,p3 initial termal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	

  //    transform to local comoving frame of the fluid
//  cout << endl;
//  cout << "flow  "<< v0[1] << " " << v0[2] << " " << v0[3] << " "<<" Elab " << p0[0] << endl;
  
  trans(v0,p0);
//  cout << p0[0] << " " << sqrt(qhat0ud) << endl;
  
//  cout << sqrt(pow(p0[1],2)+pow(p0[2],2)+pow(p0[3],2)) << " " << p0[1] << " " << p0[2] << " " << p0[3] << endl;

  //************************************************************
  trans(v0,p4);
  //************************************************************


  //    sample the medium parton thermal momentum in the comoving frame


  double xw;
  double razim;
  double rcos;
  double rsin;
	  
  double ss;  
  double tmin;
  double tmid;
  double tmax;
	  
  double rant;
  double tt;
	  
  double uu;	  
  double ff;
  double rank;
	  
  double mmax;
  double msq;
	  
  double f1;
  double f2;
  
  double p0ex[4]={0.0};

  //    Initial 4-momentum of jet
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	  

  int ic=0;

  do
	{
      do{	  
	  th:   xw=15.0*ran0(&NUM1);
		razim=2.0*pi*ran0(&NUM1);
		rcos=1.0-2.0*ran0(&NUM1);
		rsin=sqrt(1.0-rcos*rcos);
		//
		p2[0]=xw*temp;
		p2[3]=p2[0]*rcos;
		p2[1]=p2[0]*rsin*cos(razim);
		p2[2]=p2[0]*rsin*sin(razim);
	  
		f1=pow(xw,3 )/(exp(xw)-1)/1.4215;
		f2=pow(xw,3)/(exp(xw)+1)/1.2845;
		//
		//    cms energy
		//
		ss=2.0*(p0[0]*p2[0]-p0[1]*p2[1]-p0[2]*p2[2]-p0[3]*p2[3]);
	  
		//	if(ss.lt.2.d0*qhat0ud) goto 14

		tmin=qhat0ud;
		tmid=ss/2.0;
		tmax=ss-qhat0ud;
	  
		//    use (s^2+u^2)/(t+qhat0ud)^2 as scattering cross section in the
		//
		rant=ran0(&NUM1);
		tt=rant*ss;	  

//		ic+=1;
//		cout << p0[0] << "  " << p2[0] <<  endl;
//		cout << tt << "  " << ss <<  "" << qhat0ud <<endl;
//		cout << ic << endl;

	  }while((tt<qhat0ud) || (tt>(ss-qhat0ud))); 

	  uu=ss-tt;	  
	  
	  if(CT==1)
		{	
		  ff=f1;	  
		  mmax=4.0/pow(ss,2)*(3.0-tmin*(ss-tmin)/pow(ss,2)+(ss-tmin)*ss/pow(tmin,2)+tmin*ss/pow((ss-tmin),2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(3.0-tt*uu/pow(ss,2)+uu*ss/pow(tt,2)+tt*ss/pow(uu,2))/mmax;
		}

	
	  if(CT==2)
		{
		  ff=f1;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
		}

	
	  if(CT==3)
		{
		  ff=f2;	  
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin));
			}
		  else
			{	  
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		}

	
	  if(CT==13)
		{
		  ff=f1;
	  
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin));
			}
		  else 
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		}

	
	  if(CT==4)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(ss,2)+pow(uu,2))/pow(tt,2))/mmax;
		}

	
	  if(CT==5)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(ss,2)+pow(tmin,2))/pow((ss-tmin),2)-2.0/3.0*pow(ss,2)/tmin/(ss-tmin));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*((pow(ss,2)+pow(uu,2))/pow(tt,2)+(pow(ss,2)+pow(tt,2))/pow(uu,2)-2.0/3.0*pow(ss,2)/tt/uu))/mmax;	  
		}

	
	  if(CT==6)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+0.5);
		}

	
	  if(CT==7)
		{
		  ff=f2;	  
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2)+2.0/3.0*pow((ss-tmin),2)/ss/tmin);
		  msq=(pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(((pow(ss,2)+pow(uu,2))/pow(tt,2))+(pow(tt,2)+pow(uu,2))/pow(ss,2)+2.0/3.0*pow(uu,2)/ss/tt)))/mmax;	  
		}

	
	  if(CT==8)
		{
		  ff=f2;	 
		  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
		  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);	  
		}
	
	  rank=ran0(&NUM1);
	}while(rank>(msq*ff));
	
  //
  p3[1]=p2[1]; 
  p3[2]=p2[2];
  p3[3]=p2[3];
  p3[0]=p2[0];

  //    velocity of the center-of-mass
  //
  vc[1]=(p0[1]+p2[1])/(p0[0]+p2[0]);
  vc[2]=(p0[2]+p2[2])/(p0[0]+p2[0]);
  vc[3]=(p0[3]+p2[3])/(p0[0]+p2[0]);
  //
  //    transform into the cms frame
  //
  trans(vc,p0);
  trans(vc,p2);
  //
  //    cm momentum
  //
  double pcm=p2[0];
  //
  //    sample transverse momentum transfer with respect to jet momentum
  //    in cm frame
  //
  double ranp=2.0*pi*ran0(&NUM1);
  //
  //    transverse momentum transfer
  //
  qt=sqrt(pow(pcm,2)-pow((tt/2.0/pcm-pcm),2));
  double qx=qt*cos(ranp);
  double qy=qt*sin(ranp);

  //
  //    longitudinal momentum transfer
  //
  double qpar=tt/2.0/pcm;
  //
  //    qt is perpendicular to pcm, need to rotate back to the cm frame
  //
  double upt=sqrt(p2[1]*p2[1]+p2[2]*p2[2])/p2[0];
  double upx=p2[1]/p2[0];
  double upy=p2[2]/p2[0];
  double upz=p2[3]/p2[0];
  //
  //    momentum after collision in cm frame
  //
  p2[1]=p2[1]-qpar*upx;
  p2[2]=p2[2]-qpar*upy;
  if(upt!=0.0) 
	{
      p2[1]=p2[1]+(upz*upx*qy+upy*qx)/upt;
      p2[2]=p2[2]+(upz*upy*qy-upx*qx)/upt;
	}
 s2:  p2[3]=p2[3]-qpar*upz-upt*qy;

  p0[1]=-p2[1];
  p0[2]=-p2[2];
  p0[3]=-p2[3];
  //
  //    transform from cm back to the comoving frame
  //
  transback(vc,p2);
  transback(vc,p0);

  //************************************************************
  //
  //     calculate qt in the rest frame of medium
  //
//  if(p0[4]>p2[4])
//	{
	  rotate(p4[1],p4[2],p4[3],p0,1);
	  qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
	  rotate(p4[1],p4[2],p4[3],p0,-1);
//	}
//  else
//	{
//	  rotate(p4[1],p4[2],p4[3],p2,1);
//	  qt=sqrt(pow(p2[1],2)+pow(p2[2],2));
//	  rotate(p4[1],p4[2],p4[3],p2,-1);
//	}
  //************************************************************	  
	  
	  
	  
  //
  //    transform from comoving frame to the lab frame
  //
  transback(v0,p2);
  transback(v0,p0);
  transback(v0,p3);
	  
  //************************************************************
  transback(v0,p4);
  //************************************************************	  
	  
}

//.........................................................................
void collHQ22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt){
  //
  //    HQ 2->2 scatterings
  //    p0 initial HQ momentum, output to final momentum
  //    p2 final thermal momentum, p3 initial thermal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************


  // transform to local comoving frame of the fluid
  trans(v0,p0);

  //************************************************************

  //    sample the medium parton thermal momentum in the comoving frame

  double xw;
  double razim;
  double rcos;
  double rsin;
	  
  double ss;  
	  
  double rant;
  double tt;
	  
  double uu;	  
  double ff;
  double rank;
	  
  double msq;
	  
  double e2,theta2,theta4,phi24;   // the four independent variables      
  double e1,e4,p1,cosTheta24,downFactor,sigFactor; // other useful variables
  double HQmass,fBmax,fFmax,fB,fF,maxValue;
  int index_p1,index_T,index_e2;
  int ct1_loop,ct2_loop,flag1,flag2;

  flag1=0;
  flag2=0;

// continue this function for HQ scattering

  HQmass=p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3];
  if(HQmass>1e-12) {
      HQmass=sqrt(HQmass);
  } else {
      HQmass = 0.0;
  }

  //    Initial 4-momentum of HQ
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	  

  p1=sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]);
  index_p1=(int)((p1-min_p1)/bin_p1);
  index_T=(int)((temp-min_T)/bin_T);
  if(index_p1>=N_p1) {
      index_p1=N_p1-1;
      cout << "warning: p1 is over p_max: " << p1 << endl;
  }
  if(index_T>=N_T) {
      index_T=N_T-1;
      cout << "warning: T is over T_max: " << temp << endl;
  } 
  if(index_T<0) {
      index_T=0;
      cout << "warning: T is below T_min: " << temp << endl;
  } 

  fBmax=distFncBM[index_T][index_p1];
  fFmax=distFncFM[index_T][index_p1];  // maximum of f(xw) at given p1 and T

  maxValue=10.0;  // need actual value later
	
  ct1_loop=0;
  do {   // sample p2 (light parton) using distribution integrated over 3 angles
        ct1_loop++;
        if(ct1_loop>1e6) {
//            cout << "cannot sample light parton for HQ scattering ..." << endl;
            flag1=1;
            break;
        }
        xw=max_e2*ran0(&NUM1);
        index_e2=(int)((xw-min_e2)/bin_e2);
        if(index_e2>=N_e2) index_e2=N_e2-1;
        if(CT==11) { // qc->qc
            ff=distFncF[index_T][index_p1][index_e2]/fFmax;
            maxValue=distMaxF[index_T][index_p1][index_e2];
        } else if(CT==12) { // gc->gc
            ff=distFncB[index_T][index_p1][index_e2]/fBmax;
            maxValue=distMaxB[index_T][index_p1][index_e2];
        } else {
            cout << "Wrong HQ channel ID" << endl;
            exit(EXIT_FAILURE);
        }
  } while(ran0(&NUM1)>ff);
      
  e2=xw*temp;
  e1=p0[0];

  // now e2 is fixed, need to sample the remaining 3 variables
  ct2_loop=0;
  do {
          ct2_loop++;
          if(ct2_loop>1e6) {
              cout << "cannot sample final states for HQ scattering ..." << endl;
              flag2=1;
              break;
          }

          theta2=pi*ran0(&NUM1);
          theta4=pi*ran0(&NUM1);
          phi24=2.0*pi*ran0(&NUM1);

          cosTheta24=sin(theta2)*sin(theta4)*cos(phi24)+cos(theta2)*cos(theta4);
          downFactor=e1-p1*cos(theta4)+e2-e2*cosTheta24;
          e4=(e1*e2-p1*e2*cos(theta2))/downFactor;
          sigFactor=sin(theta2)*sin(theta4)*e2*e4/downFactor; 

          // calculate s,t,u, different definition from light quark -- tt, uu are negative
          ss=2.0*e1*e2+HQmass*HQmass-2.0*p1*e2*cos(theta2);
          tt=-2.0*e2*e4*(1.0-cosTheta24);
          uu=2.0*HQmass*HQmass-ss-tt;

          // re-sample if the kinematic cuts are not satisfied
          if(ss<=2.0*qhat0ud || tt>=-qhat0ud || uu>=-qhat0ud) {
	       rank=ran0(&NUM1);
               sigFactor=0.0;
               msq=0.0;
               continue;
          }

	  if(CT==11) {  // qc->qc
               ff=(1.0/(exp(e2/temp)+1.0))*(1.0-1.0/(exp(e4/temp)+1.0));
               sigFactor=sigFactor*ff;
 	       msq=Mqc2qc(ss,tt,HQmass)/maxValue;
	  }

          if(CT==12) {  // gc->gc
               ff=(1.0/(exp(e2/temp)-1.0))*(1.0+1.0/(exp(e4/temp)-1.0));
               sigFactor=sigFactor*ff;
 	       msq=Mgc2gc(ss,tt,HQmass)/maxValue;
	  }

	  rank=ran0(&NUM1);

  } while(rank>(msq*sigFactor));
	
  if(flag1==0 && flag2==0) {

      // pass p2 value to p3 for initial thermal parton
      p3[1]=e2*sin(theta2);
      p3[2]=0.0;
      p3[3]=e2*cos(theta2);
      p3[0]=e2;
    
      // calculate momenta of outgoing particles
      // here p2 is for p4 (light parton) in my note
    
      p2[1]=e4*sin(theta4)*cos(phi24);
      p2[2]=e4*sin(theta4)*sin(phi24);
      p2[3]=e4*cos(theta4);
      p2[0]=e4;
    
      // Because we treated p0 (p1 in my note for heavy quark) as the z-direction, proper rotations are necessary here
      rotate(p4[1],p4[2],p4[3],p2,-1);
      rotate(p4[1],p4[2],p4[3],p3,-1);
	
      p0[1]=p4[1]+p3[1]-p2[1];
      p0[2]=p4[2]+p3[2]-p2[2];
      p0[3]=p4[3]+p3[3]-p2[3];
      p0[0]=sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]+HQmass*HQmass);
    
      // Debug
      if(fabs(p0[0]+p2[0]-p3[0]-p4[0])>0.00001) {
         cout << "Violation of energy conservation in HQ 2->2 scattering:  " << fabs(p0[0]+p2[0]-p3[0]-p4[0]) << endl;
      }
      
      // calculate qt in the rest frame of medium
      rotate(p4[1],p4[2],p4[3],p0,1);
      qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
      rotate(p4[1],p4[2],p4[3],p0,-1);

      // transform from comoving frame to the lab frame
      transback(v0,p2);
      transback(v0,p0);
      transback(v0,p3);
      transback(v0,p4);

  } else { // no scattering
      transback(v0,p0);
      transback(v0,p4);
      qt=0;
      p2[0]=0;
      p2[1]=0;
      p2[2]=0;
      p2[3]=0;
      p3[0]=0;
      p3[1]=0;
      p3[2]=0;
      p3[3]=0;
  }

	  	 	  
}

void twcoll(int CT,double qhat0ud,double v0[4],double p0[4],double p2[4]){    
  //	
  //     p0 initial jet momentum, output to final momentum
  //     p2 final thermal momentum,p3 initial thermal energy
  //
  //     amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //     transform to local comoving frame of the fluid
       
  trans(v0,p0);
  trans(v0,p2);

  //    velocity of the center-of-mass
			

  vc[1]=(p0[1]+p2[1])/(p0[0]+p2[0]);
  vc[2]=(p0[2]+p2[2])/(p0[0]+p2[0]);
  vc[3]=(p0[3]+p2[3])/(p0[0]+p2[0]);
  //
  //    transform into the cms frame
  //

  trans(vc,p0);
  trans(vc,p2);	   

  //
  //     cm momentum
  //
  double pcm=p2[0];	   	   
  //
  //     sample transverse momentum transfer with respect to jet momentum
  //     in cm frame
  //
  //
  //     Gaussian distribution
  //
  //     qt=sqrt(-dt*qhat0ud*log(1-rant+rant*exp(-scm/(4.d0*dt*qhat0ud))))
  //
  //     static potential distribution
  //	   
  double ss=4.0*pow(pcm,2);
  //
  double tmin=qhat0ud;
  double tmid=ss/2.0;
  double tmax=ss-qhat0ud;
	   
  double rant;
  double tt;
  double uu;
  double mmax;
  double msq;
  double rank;	

  /////////////////////////////////////////////
  double ranp;
  double qt;
  double qx;
  double qy;
  double qpar;
  double upt;
  double upx;
  double upy;
  double upz;
  /////////////////////////////////////////////	   
	   
  //
  //	   CT is a variable notated different collision types.
  //
        
  do
	{	  
	tw:     rant=ran0(&NUM1);
	  tt=rant*ss;
		
	  if((tt<qhat0ud) || (tt>(ss-qhat0ud)))
		break;
	  //      if((tt<qhat0ud) || (tt>(ss-qhat0ud)))
	  //      {		
	  //		goto t59;
	  //      }		
	  uu=ss-tt;

	  //	   gg to gg
	  if(CT==1)
		{
		  //		
		  mmax=3.0-tmin*(ss-tmin)/pow(ss,2)+(ss-tmin)*ss/pow(tmin,2)+tmin*ss/pow((ss-tmin),2);
		  msq=(3.0-tt*uu/pow(ss,2)+uu*ss/pow(tt,2)+tt*ss/pow(uu,2))/mmax;
		  //
		}
	

	  //	   gg to qqbar
	  if(CT==2)
		{
		  //
		  mmax=4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2);
		  msq=(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
		  //
		}

	  //	   gq to gq, gqbar to gqbar
	  if(CT==3)
		{
		  //
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin);
			}
		  else
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		  //
		}
	   
	  //	   qg to qg, qbarg to qbarg
	  if(CT==13)
		{	   
		  //
		  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
			{
			  mmax=(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin);		
			}
		  else
			{
			  mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
			}
		  //
		  msq=((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
		  //	    
		}

	  //	   qiqj to qiqj, qiqjbar to qiqjbar, qibarqj to qibarqj, qibarqjbar to qibarqjbar
	  //	   for i not equal j
	  if(CT==4)
		{	  
		  //		
		  mmax=4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2);
		  msq=(4.0/9.0*(pow(ss,2)+pow(uu,2))/pow(tt,2))/mmax;
		  //
		}

	  //	   qiqi to qiqi, qibarqibar to qibarqibar
	  if(CT==5)
		{
		  //	 
		  mmax=4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(ss,2)+pow(tmin,2))/pow((ss-tmin),2)-2.0/3.0*pow(ss,2)/tmin/(ss-tmin);
		  msq=(4.0/9.0*((pow(ss,2)+pow(uu,2))/pow(tt,2)+(pow(ss,2)+pow(tt,2))/pow(uu,2)-2.0/3.0*pow(ss,2)/tt/uu))/mmax;
		  //	  
		}

	  //     qiqibar to qjqjbar for i not equal j
	  if(CT==6)
		{
		  //	 
		  mmax=4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2);
		  msq=(4.0/9.0*(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+0.5);
		  //
		}

	  //	   qiqibar to qiqibar
	  if(CT==7)
		{
		  //	 
		  mmax=4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2)+2.0/3.0*pow((ss-tmin),2)/ss/tmin;
		  msq=(4.0/9.0*(((pow(ss,2)+pow(uu,2))/pow(tt,2))+(pow(tt,2)+pow(uu,2))/pow(ss,2)+2.0/3.0*pow(uu,2)/ss/tt))/mmax;
		  //
		}

	  //	   qqbar to gg
	  if(CT==8)
		{
		  //	 
		  mmax=4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2);
		  msq=(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
		  //
		}
	   
	  rank=ran0(&NUM1);
      

	}while(rank>msq);	   	  
	  	  
  ///////////////////////////////////////////////////////////////////////
  //    transverse momentum transfer

  if((tt>qhat0ud) && (tt<(ss-qhat0ud)))
	{

      ranp=2.0*pi*ran0(&NUM1);
	  //
	  //
	  //
      qt=sqrt(pow(pcm,2)-pow((tt/2.0/pcm-pcm),2));
      qx=qt*cos(ranp);
      qy=qt*sin(ranp);

	  //
	  //    longitudinal momentum transfer
	  //
      qpar=tt/2.0/pcm;
	  //
	  //    qt is perpendicular to pcm, need to rotate back to the cm frame
	  //
      upt=sqrt(p2[1]*p2[1]+p2[2]*p2[2])/p2[0];
      upx=p2[1]/p2[0];
      upy=p2[2]/p2[0];
      upz=p2[3]/p2[0];
	  //
	  //    momentum after collision in cm frame
	  //

      p2[1]=p2[1]-qpar*upx;
      p2[2]=p2[2]-qpar*upy;
      if(upt!=0.0) 
		{
		  p2[1]=p2[1]+(upz*upx*qy+upy*qx)/upt;
		  p2[2]=p2[2]+(upz*upy*qy-upx*qx)/upt;
		}
	s3:   p2[3]=p2[3]-qpar*upz-upt*qy;

      p0[1]=-p2[1];
      p0[2]=-p2[2];
      p0[3]=-p2[3];
	  
	}
  //
  //    transform from cm back to the comoving frame
  //
 t59:  transback(vc,p2);
       transback(vc,p0);
  //
  //    transform from comoving frame to the lab frame
  //

  transback(v0,p2);
  transback(v0,p0);
  transback(v0,p3);
	  
}


void collHQ23(int parID, double temp_med, double qhat0ud, double v0[4], double p0[4], double p2[4], double p3[4], double p4[4], double qt, int &ic, double Tdiff, double HQenergy, double max_Ng, double xLow, double xInt) {

  //    p0 initial jet momentum, output to final momentum
  //    p3 initial thermal momentum
  //    p2 initial thermal momentum, output to final thermal momentum
  //    p4 radiated gluon momentum
  //    qt transverse momentum transfer in the rest frame of medium
  //    q0,ql energy and longitudinal momentum transfer
  //    i=0: 2->3 finished; i=1: can not find a gluon when nloop<=30
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
	
  double randomX,randomY;
  int count_sample,nloopOut,flagOut;
  double theta_gluon,kperp_gluon;
  double kpGluon[4];
  double zDirection[4];
  double HQmass=sqrt(p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]);

  if(abs(parID)!=4) {
      HQmass = 0.0;
      p0[0] = sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]+HQmass*HQmass);
  }   

  ic=0;
  nloopOut=0;
  flagOut=0;

  //    initial thermal parton momentum in lab frame
  p2[1]=p3[1];
  p2[2]=p3[2];
  p2[3]=p3[3];
  p2[0]=p3[0];

//  cout<<"check1 p0[0]: "<<p0[0]<<endl;

  //    transform to local comoving frame of the fluid
  trans(v0,p0);
  trans(v0,p2);

  if(p0[0]<2.0*sqrt(qhat0ud)){
      ic=1;
      return;
  }

  // define z-direction
  zDirection[1]=p0[1];
  zDirection[2]=p0[2];
  zDirection[3]=p0[3];
  zDirection[0]=p0[0];

  rotate(zDirection[1],zDirection[2],zDirection[3],p2,1); // rotate p2 into p1 system

//  cout<<"check2 p0[0]: "<<p0[0]<<endl;

  do {

      do {
          randomX=xLow+xInt*ran0(&NUM1);
          randomY=ran0(&NUM1);
      } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
     
      count_sample=0;
      while(max_Ng*ran0(&NUM1)>dNg_over_dxdydt(parID,randomX,randomY,HQenergy,HQmass,temp_med,Tdiff)) {
          count_sample=count_sample+1;
          if(count_sample>1e+5) {
              cout << "give up loop at point 1 ..." << endl;
              kpGluon[1]=0.0;
              kpGluon[2]=0.0;
              kpGluon[3]=0.0;
              kpGluon[0]=0.0;
              ic=1;
    //          break;
              return;
          }
    
          do {
             randomX=xLow+xInt*ran0(&NUM1);
             randomY=ran0(&NUM1);
          } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
      }
    
      if(parID==21&&randomX>0.5) randomX=1.0-randomX;
      theta_gluon=2.0*pi*ran0(&NUM1);
      kperp_gluon=randomX*randomY*HQenergy;
      kpGluon[1]=kperp_gluon*cos(theta_gluon);
      kpGluon[2]=kperp_gluon*sin(theta_gluon);
      kpGluon[3]=randomX*HQenergy*sqrt(1.0-randomY*randomY);
      kpGluon[0]=sqrt(kpGluon[1]*kpGluon[1]+kpGluon[2]*kpGluon[2]+kpGluon[3]*kpGluon[3]);
      
//      rotate(p4[1],p4[2],p4[3],kpGluon,-1);
    
      if(kpGluon[0]>(HQenergy-HQmass)) {
          kpGluon[1]=0.0;
          kpGluon[2]=0.0;
          kpGluon[3]=0.0;
          kpGluon[0]=0.0;
          nloopOut++;
          continue;
      }
    
      // update mometum
    
      p4[1]=kpGluon[1];
      p4[2]=kpGluon[2];
      p4[3]=kpGluon[3];
      p4[0]=kpGluon[0];
    
//      rotate(zDirection[1],zDirection[2],zDirection[3],p2,1); // rotate p2 into p1 system

    // solve energy-momentum conservation
      double sE1=p0[0];
//      double sp1x=p0[1];
//      double sp1y=p0[2];
//      double sp1z=p0[3];
      double sp1x=0.0;
      double sp1y=0.0;
      double sp1z=sqrt(p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]);
      double sE2=p2[0];
      double sp2x=p2[1];
      double sp2y=p2[2];
      double sp2z=p2[3];
      double sk0=p4[0];
      double skx=p4[1];
      double sky=p4[2];
      double skz=p4[3];
      double sqx,sqy,sqzA,sq0A,sqzB,sq0B,sqz,sq0,sqtheta;
      double sAA,sBB,sCC,aaa,bbb,ccc,abc,abc2;
      int nloop1=0;
      int nloop2=0;
      int flagDone=0;
      int yesA,yesB;
     
      do {
          sqtheta=2.0*pi*ran0(&NUM1);
          sqx=qt*cos(sqtheta);
          sqy=qt*sin(sqtheta);
          sAA=(sE1+sE2-sk0)/(sp1z+sp2z-skz);
          sBB=(pow(sE2,2)-pow((sE1-sk0),2)+pow((sp1x-sqx-skx),2)+pow((sp1y-sqy-sky),2)+pow((sp1z-skz),2)+pow(HQmass,2)-pow((sp2x+sqx),2)-pow((sp2y+sqy),2)-pow(sp2z,2))/2.0/(sp1z+sp2z-skz);
          aaa=sAA*sAA-1.0;
          bbb=2.0*(sAA*sp2z+sAA*sBB-sE2);
          ccc=pow((sp2x+sqx),2)+pow((sp2y+sqy),2)+sp2z*sp2z+2.0*sp2z*sBB+sBB*sBB-sE2*sE2;
          abc2=bbb*bbb-4.0*aaa*ccc;
     
          if(abc2<0.0) {
              nloop1++;
              continue;
          } else {
              nloop1=0;
          }
     
          abc=sqrt(abc2);
          sq0A=(-bbb+abc)/2.0/aaa;
          sq0B=(-bbb-abc)/2.0/aaa;
          sqzA=sAA*sq0A+sBB;
          sqzB=sAA*sq0B+sBB;
    
    // require space-like and E_final > M;
          if(sq0A*sq0A-sqx*sqx-sqy*sqy-sqzA*sqzA<0 && sE1-sq0A-sk0>HQmass && sE2+sq0A>0) {
              yesA=1; 
          } else {
              yesA=0;
          }
          if(sq0B*sq0B-sqx*sqx-sqy*sqy-sqzB*sqzB<0 && sE1-sq0B-sk0>HQmass && sE2+sq0B>0) {
              yesB=1; 
          } else {
              yesB=0;
          }
                   
    // select appropriate solution
          if(yesA==0 && yesB==0) {
    //          cout << "Solutions fail ..." << endl;
              nloop2++;
              continue;
          } else if(yesA==1 && yesB==1) {
    //          cout << "Both solutions work!" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              if(abs(sq0A)<abs(sq0B)) {
                  sq0=sq0A;
                  sqz=sqzA;
              } else {
                  sq0=sq0B;
                  sqz=sqzB;
              }
          } else if(yesA==1) {
    //          cout << "pass A ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sq0=sq0A;
              sqz=sqzA;
          } else {
    //          cout << "pass B ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sq0=sq0B;
              sqz=sqzB;
          }
    
          flagDone=1;   
      } while (flagDone==0 && nloop1<loopN && nloop2<loopN);

      if(flagDone==0) { 
    //      ic=1;  // no appropriate solution
    //      cout << "solution fails ..." << endl;
          nloopOut++;
          continue;
      } else {
          p0[0]=sE1-sq0-sk0;
          p0[1]=sp1x-sqx-skx;
          p0[2]=sp1y-sqy-sky;
          p0[3]=sp1z-sqz-skz;
     
          p2[0]=sE2+sq0;
          p2[1]=sp2x+sqx;
          p2[2]=sp2y+sqy;
          p2[3]=sp2z+sqz;
        
    //           cout<<"sE1,sE2,p0[0],p2[0],p4[0]: "<<sE1<<"  "<<sE2<<"  "<<p0[0]<<"  "<<p2[0]<<"  "<<p4[0]<<endl;

          rotate(zDirection[1],zDirection[2],zDirection[3],p0,-1); // rotate p0 into global system
          rotate(zDirection[1],zDirection[2],zDirection[3],p2,-1); // rotate p2 into global system
          rotate(zDirection[1],zDirection[2],zDirection[3],p4,-1); // rotate p4 into global system
    
          transback(v0,p0);
          transback(v0,p2);
          transback(v0,p4);
        
    // debug: check on-shell condition
          if(abs(p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]-HQmass*HQmass)>0.000001 || abs(p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]-p2[3]*p2[3])>0.000001) {
              cout << "Wrong solution -- not on shell" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sq0 << "  " << sqx << "  " << sqy << "  " << sqz << endl;
              cout << abs(p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]-HQmass*HQmass) << "  " << abs(p2[0]*p2[0]-p2[1]*p2[1]-p2[2]*p2[2]-p2[3]*p2[3]) << "  " << HQmass << endl;
          }

//  cout<<"in colljet23 -- init thermal, final jet, final thermal, gluon: "<<p3[0]<<"  "<<p0[0]<<"  "<<p2[0]<<"  "<<p4[0]<<endl;
          flagOut=1;
      }
  } while(flagOut==0 && nloopOut<loopN);

  if(flagOut==0) ic=1;   

}

void radiationHQ(int parID, double qhat0ud, double v0[4], double P2[4], double P3[4], double P4[4], double Pj0[4], int &ic, double Tdiff, double HQenergy, double max_Ng, double temp_med, double xLow, double xInt){  

  //    work in the rest frame of medium
  //    return the 4-momentum of final states in 1->3 radiation
  //    input: P2(4)-momentum of radiated gluon from 2->3
  //           P3(4)-momentum of daughter parton from 2->3
  //           Pj0(4)-inital momentum of jet before 2->3
  //           v0-local velocity of medium
  //    output:P2(4)-momentum of 1st radiated gluon
  //           P3(4)-momentum of daughter parton
  //           P4(4)-momentum of 2nd radiated gluon
  //           i=1: no radiation; 0:successful radiation


  double randomX,randomY;
  int count_sample,nloopOut,flagOut;
  double theta_gluon,kperp_gluon;
  double kpGluon[4];
  double HQmass=sqrt(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]);

  if(abs(parID)!=4) {
      HQmass = 0.0;
      P3[0] = sqrt(P3[1]*P3[1]+P3[2]*P3[2]+P3[3]*P3[3]);
  }   


  double P50[4]={0.0};
  double P2i[4]={0.0};
  double P3i[4]={0.0};	  
	  
  ic=0;

  for(int i=0;i<=3;i++){
      P2i[i]=P2[i];
      P3i[i]=P3[i];
  }
	  
  // transform to local comoving frame of the fluid
  trans(v0,P2);
  trans(v0,P3);
  trans(v0,Pj0);

  xInt=(P3[0]-HQmass)/HQenergy-xLow;
  if(xInt<=0.0) {
      ic=1;
      return;
  }

  double px0=P2[1]+P3[1];
  double py0=P2[2]+P3[2];
  double pz0=P2[3]+P3[3];

  // rotate to the frame in which jet moves along z-axis
  rotate(px0,py0,pz0,P3,1);
  rotate(px0,py0,pz0,P2,1);
  rotate(px0,py0,pz0,Pj0,1);
  
  for(int i=0;i<=3;i++){
      P50[i]=P2[i]+P3[i];
  }
	
  nloopOut=0;
  flagOut=0;
  // sample the second gluon
  do {
      do {
          randomX=xLow+xInt*ran0(&NUM1);
          randomY=ran0(&NUM1);
      } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
     
      count_sample=0;
      while(max_Ng*ran0(&NUM1)>dNg_over_dxdydt(parID,randomX,randomY,HQenergy,HQmass,temp_med,Tdiff)) {
          count_sample=count_sample+1;
          if(count_sample>1e+5) {
              cout << "give up loop at point 1 ..." << endl;
              kpGluon[1]=0.0;
              kpGluon[2]=0.0;
              kpGluon[3]=0.0;
              kpGluon[0]=0.0;
              ic=1;
    //          break;
              return;
          }
    
          do {
             randomX=xLow+xInt*ran0(&NUM1);
             randomY=ran0(&NUM1);
          } while(tau_f(randomX,randomY,HQenergy,HQmass)<1.0/pi/temp_med);   
      }
    
      if(parID==21&&randomX>0.5) randomX=1.0-randomX;
      theta_gluon=2.0*pi*ran0(&NUM1);
      kperp_gluon=randomX*randomY*HQenergy;
      kpGluon[1]=kperp_gluon*cos(theta_gluon);
      kpGluon[2]=kperp_gluon*sin(theta_gluon);
      kpGluon[3]=randomX*HQenergy*sqrt(1.0-randomY*randomY);
      kpGluon[0]=sqrt(kpGluon[1]*kpGluon[1]+kpGluon[2]*kpGluon[2]+kpGluon[3]*kpGluon[3]);
      
      rotate(Pj0[1],Pj0[2],Pj0[3],kpGluon,-1);
    
      if(kpGluon[0]>(P3[0]-HQmass)) { // which should be impossible due to reset of xInt above
          kpGluon[1]=0.0;
          kpGluon[2]=0.0;
          kpGluon[3]=0.0;
          kpGluon[0]=0.0;
          nloopOut++;
          continue;
      }
    
      // update mometum
    
      P4[1]=kpGluon[1];
      P4[2]=kpGluon[2];
      P4[3]=kpGluon[3];
      P4[0]=kpGluon[0];
    
    // solve energy-momentum conservation
    // my notation in the note
    // p0 is re-constructed off-shell parton, correponding to P5
    // p1 is the final heavy quark from 2->3, corresponding to P3 (input). P3 (output) is the final heavy quark after 1->3.
    // k1 is the first gluon, corresponding to P2
    // k2 is the second gluon, corresponding to P4
    // assume k10 and p10 unchanged and modify their other components while p0 and k2 are fixed
      double sp0z=sqrt(P50[1]*P50[1]+P50[2]*P50[2]+P50[3]*P50[3]);
      double sk10=P2[0];
      double sk1zOld=P2[3];
      double sk1z,sk1p,sk1x,sk1y,sktheta; // unknown
      double sp10=P3[0];
      double sk20=P4[0];
      double sk2x=P4[1];
      double sk2y=P4[2];
      double sk2z=P4[3];
      double sk2p=sqrt(sk2x*sk2x+sk2y*sk2y);
    
      double stheta12,cos12Min2;
      double sAA,aaa,bbb,ccc,abc,abc2;
      double sk1z1,sk1z2,sk1p1,sk1p2;
      int nloop1=0;
      int nloop2=0;
      int flagDone=0;
      int yesA,yesB;
      
      sAA=sk10*sk10+sk2p*sk2p+sp0z*sp0z+sk2z*sk2z-2.0*sp0z*sk2z+HQmass*HQmass-(sp10-sk20)*(sp10-sk20);
      cos12Min2=(sAA*sAA/4.0/sk10/sk10-(sp0z-sk2z)*(sp0z-sk2z))/sk2p/sk2p;
    //  cout << "cos^2 min: " << cos12Min2 << endl;
    
      if(cos12Min2>1.0) {
          nloopOut++;
          continue;
    //      cout << "cos^2 min is outside the kinematic range: " << cos12Min2 << endl;
    //      return;
      }       
    
      do {
          stheta12=2.0*pi*ran0(&NUM1); // theta between k1 and k2
          aaa=4.0*((sp0z-sk2z)*(sp0z-sk2z)+sk2p*sk2p*cos(stheta12)*cos(stheta12));
          bbb=-4.0*sAA*(sp0z-sk2z);
          ccc=sAA*sAA-4.0*sk10*sk10*sk2p*sk2p*cos(stheta12)*cos(stheta12);
    
          abc2=bbb*bbb-4.0*aaa*ccc;
    
          if(abc2<0.0) {
              nloop1++;
              continue;
          } else {
              nloop1=0;
          }
    
          abc=sqrt(abc2);
          sk1z1=(-bbb+abc)/2.0/aaa;
          sk1z2=(-bbb-abc)/2.0/aaa;
          sk1p1=sqrt(sk10*sk10-sk1z1*sk1z1);
          sk1p2=sqrt(sk10*sk10-sk1z2*sk1z2);
    
    // Since we have squared both sides of the original equation during solving, the solutions may not satisfy the original equation. Double check is needed.
    
    
    // require time-like of p1 and k10>k1z;
          if(2.0*sk1p1*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z1+sAA<0.000001 && sk10>abs(sk1z1) && sp10*sp10-(sp0z-sk1z1)*(sp0z-sk1z1)-(sk10*sk10-sk1z1*sk1z1)>HQmass*HQmass) {
              yesA=1; 
          } else {
              yesA=0;
          }
          if(2.0*sk1p2*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z2+sAA<0.000001 && sk10>abs(sk1z2) && sp10*sp10-(sp0z-sk1z2)*(sp0z-sk1z2)-(sk10*sk10-sk1z2*sk1z2)>HQmass*HQmass) {
              yesB=1; 
          } else {
              yesB=0;
          }
              
    // select appropriate solution
          if(yesA==0 && yesB==0) {
    //          cout << "Solutions fail ..." << endl;
              nloop2++;
              continue;
          } else if(yesA==1 && yesB==1) {
    //          cout << "Both solutions work!" << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              if(abs(sk1z1-sk1zOld)<abs(sk1z2-sk1zOld)) {
                  sk1z=sk1z1;
              } else {
                  sk1z=sk1z2;
              }
          } else if(yesA==1) {
    //          cout << "pass A ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sk1z=sk1z1;
          } else {
    //          cout << "pass B ..." << "  " << sE1 << "  " << sp1x << "  " << sp1y << "  " << sp1z << "  " << sE2 << "  " << sp2x << "  " << sp2y << "  " << sp2z << "  " << sk0 << "  " << skx << "  " << sky << "  " << skz << "  " << sqx << "  " << sqy << "  " << sq0A << "  " << sqzA << "  " << sq0B << "  " << sqzB << endl;
              sk1z=sk1z2;
          }
    
          flagDone=1;   
    
          sk1p=sqrt(sk10*sk10-sk1z*sk1z);
    //           cout << "check solution: " << pow(sp10-sk20,2)-pow(sk1p,2)-pow(sk2p,2)-2.0*sk1p*sk2p*cos(stheta12)-pow(sp0z-sk1z-sk2z,2)-pow(HQmass,2) <<"  "<< aaa*sk1z*sk1z+bbb*sk1z+ccc <<"  "<<2.0*sk1p*sk2p*cos(stheta12)-2.0*(sp0z-sk2z)*sk1z+sAA<<"  "<<2.0*sk1p*sk2p*cos(stheta12)+2.0*(sp0z-sk2z)*sk1z-sAA << endl;
    
      } while (flagDone==0 && nloop1<loopN && nloop2<loopN);
    
      if(flagDone==0) { 
   //       ic=1;  // no appropriate solution
   //       cout << "solution fails ...  " << nloop1 <<"  "<<nloop2<< endl;
          nloopOut++;
          continue;
      } else {
          sktheta=atan2(sk2y,sk2x); // theta for k2
          sktheta=sktheta+stheta12; // theta for k1
          sk1p=sqrt(sk10*sk10-sk1z*sk1z);
          sk1x=sk1p*cos(sktheta);
          sk1y=sk1p*sin(sktheta);
    
          P2[0]=sk10;
          P2[1]=sk1x;
          P2[2]=sk1y;
          P2[3]=sk1z;
          
          P3[0]=sp10-sk20;
          P3[1]=-sk1x-sk2x;
          P3[2]=-sk1y-sk2y;
          P3[3]=sp0z-sk1z-sk2z;
    
        // rotate back to the local coordinate
          rotate(px0,py0,pz0,P2,-1);
          rotate(px0,py0,pz0,P3,-1);
          rotate(px0,py0,pz0,P4,-1);
          rotate(px0,py0,pz0,Pj0,-1);
    
        // boost back to the global frame 
          transback(v0,P2);
          transback(v0,P3);
          transback(v0,P4);
          transback(v0,Pj0);
        
        // debug: check on-shell condition of P2, P3 and P4
          if(abs(P2[0]*P2[0]-P2[1]*P2[1]-P2[2]*P2[2]-P2[3]*P2[3])>0.000001 || abs(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]-HQmass*HQmass)>0.000001 || abs(P4[0]*P4[0]-P4[1]*P4[1]-P4[2]*P4[2]-P4[3]*P4[3])>0.000001) {
              cout << "Wrong solution -- not on shell" << "  " << sk10 << "  " << sk1x << "  " << sk1y << "  " << sk1z << "  " << sk20 << "  " << sk2x << "  " << sk2y << "  " << sk2z << "  " << stheta12 << "  " << sp10 << "  " << sp10-sk20 << "  " << -sk1x-sk2x << "  " << -sk1y-sk2y << "  " << sp0z << "  " << sp0z-sk1z-sk2z << "  " <<HQmass<< "  "<<pow(sp10-sk20,2)-pow(sk1x+sk2x,2)-pow(sk1y+sk2y,2)-pow(sp0z-sk1z-sk2z,2)-pow(HQmass,2)<<endl;
              cout << abs(P2[0]*P2[0]-P2[1]*P2[1]-P2[2]*P2[2]-P2[3]*P2[3]) <<"  "<< abs(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]-HQmass*HQmass) <<"  "<< abs(P4[0]*P4[0]-P4[1]*P4[1]-P4[2]*P4[2]-P4[3]*P4[3]) <<endl;
          }
    
    //       cout<<"in radiation HQ -- final gluon1, gluon2 & HQ: "<<P2[0]<<"  "<<P4[0]<<"  "<<P3[0]<<endl;
    
          flagOut=1;

          // SC: check energy-momentum conservation
          for(int i=0; i<=3; i++) {
              if(ic==0 && abs(P2i[i]+P3i[i]-P2[i]-P3[i]-P4[i])>0.000001) {
                  cout << "Warning: Violation of E.M. conservation!  " << i << " " << abs(P2i[i]+P3i[i]-P2[i]-P3[i]-P4[i]) << endl;
              }
          }
        
          // SC: check on-shell condition
          double shell2=abs(P2[0]*P2[0]-P2[1]*P2[1]-P2[2]*P2[2]-P2[3]*P2[3]);
          double shell3=abs(P3[0]*P3[0]-P3[1]*P3[1]-P3[2]*P3[2]-P3[3]*P3[3]-HQmass*HQmass);
          double shell4=abs(P4[0]*P4[0]-P4[1]*P4[1]-P4[2]*P4[2]-P4[3]*P4[3]);
          if(ic==0 && (shell2>0.000001 || shell3>0.000001 || shell4>0.000001)) {
              cout << "Warning: Violation of on-shell: " << shell2 << "  " << shell3 << "  " << shell4 << endl;
          }

     }

  } while(flagOut==0 && nloopOut<loopN);

  if(flagOut==0) ic=1;
 
}


void rotate(double px,double py,double pz,double pr[4],int icc){
  //     input:  (px,py,pz), (wx,wy,wz), argument (i)
  //     output: new (wx,wy,wz)
  //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
  //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

	   
  double wx,wy,wz,E,pt,w,cosa,sina,cosb,sinb;
  double wx1,wy1,wz1;	   	   

  wx=pr[1];
  wy=pr[2];
  wz=pr[3];

  E=sqrt(px*px+py*py+pz*pz);
  pt=sqrt(px*px+py*py);

  w=sqrt(wx*wx+wy*wy+wz*wz);

  if(pt==0)
	{
	  cosa=1;
	  sina=0;
	} 
  else
	{
	  cosa=px/pt;
	  sina=py/pt;
	}

  cosb=pz/E;
  sinb=pt/E;

  if(icc==1)
	{ 
	  wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
	  wy1=-wx*sina+wy*cosa;
	  wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
	}  

  else
	{
	  wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
	  wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
	  wz1=-wx*sinb+wz*cosb;
	}

  wx=wx1;
  wy=wy1;
  wz=wz1;

  pr[1]=wx;
  pr[2]=wy;
  pr[3]=wz;      

//  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);

}

	  
int KPoisson(double alambda){
  //....Generate numbers according to Poisson distribution
  //    P(lambda)=lambda**k exp(-lambda)/k!
  //    input: average number of radiated gluons lambda
  //    output: number of radiated gluons in this event

  double target,p;		

  double KKPoisson=0;
  target=exp(-alambda);
  p=ran0(&NUM1);
		
  while(p>target)
	{
	  p=p*ran0(&NUM1);
	  KKPoisson=KKPoisson+1;
	}		
  return KKPoisson;
}
	  
double nHQgluon(int parID,double dtLRF,double &time_gluon,double &temp_med,double &HQenergy,double &max_Ng){
  // gluon radiation probability for heavy quark       

  int time_num,temp_num,HQenergy_num;
  double delta_Ng;
  double rate_EGrid_low,rate_EGrid_high,max_EGrid_low,max_EGrid_high;
  double rate_T1E1,rate_T1E2,rate_T2E1,rate_T2E2,max_T1E1,max_T1E2,max_T2E1,max_T2E2;

  if(time_gluon>t_max) {
     cout << "accumulated time exceeds t_max" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
     time_gluon=t_max;
  }

  if(temp_med>temp_max) {
     cout << "temperature exceeds temp_max -- extrapolation is used" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
//     temp_med=temp_max;
  }

  if(HQenergy>HQener_max) {
     cout << "HQenergy exceeds HQener_max -- extrapolation is used" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
//     HQenergy=HQener_max;
  }

  if(temp_med<temp_min) {
     cout << "temperature drops below temp_min" << endl;
     cout << time_gluon << "    " << temp_med << "    " << HQenergy << endl;
     temp_med=temp_min;
  }

  time_num=(int)(time_gluon/delta_tg+0.5)+1;
//  temp_num=(int)((temp_med-temp_min)/delta_temp+0.5);
//  HQenergy_num=(int)(HQenergy/delta_HQener+0.5); // use linear interpolation instead of finding nearest point for E and T dimensions
  temp_num=(int)((temp_med-temp_min)/delta_temp);
  HQenergy_num=(int)(HQenergy/delta_HQener); // normal interpolation
  if(HQenergy_num >= HQener_gn) HQenergy_num=HQener_gn-1; // automatically become extrapolation
  if(temp_num >= temp_gn) temp_num=temp_gn-1;

  if(parID==21) {
      rate_T1E1 = dNg_over_dt_g[time_num][temp_num][HQenergy_num];
      rate_T1E2 = dNg_over_dt_g[time_num][temp_num][HQenergy_num+1];
      rate_T2E1 = dNg_over_dt_g[time_num][temp_num+1][HQenergy_num];
      rate_T2E2 = dNg_over_dt_g[time_num][temp_num+1][HQenergy_num+1];
      max_T1E1 = max_dNgfnc_g[time_num][temp_num][HQenergy_num];
      max_T1E2 = max_dNgfnc_g[time_num][temp_num][HQenergy_num+1];
      max_T2E1 = max_dNgfnc_g[time_num][temp_num+1][HQenergy_num];
      max_T2E2 = max_dNgfnc_g[time_num][temp_num+1][HQenergy_num+1];
  } else if(abs(parID)==4) {
      rate_T1E1 = dNg_over_dt_c[time_num][temp_num][HQenergy_num];
      rate_T1E2 = dNg_over_dt_c[time_num][temp_num][HQenergy_num+1];
      rate_T2E1 = dNg_over_dt_c[time_num][temp_num+1][HQenergy_num];
      rate_T2E2 = dNg_over_dt_c[time_num][temp_num+1][HQenergy_num+1];
      max_T1E1 = max_dNgfnc_c[time_num][temp_num][HQenergy_num];
      max_T1E2 = max_dNgfnc_c[time_num][temp_num][HQenergy_num+1];
      max_T2E1 = max_dNgfnc_c[time_num][temp_num+1][HQenergy_num];
      max_T2E2 = max_dNgfnc_c[time_num][temp_num+1][HQenergy_num+1];
  } else {
      rate_T1E1 = dNg_over_dt_q[time_num][temp_num][HQenergy_num];
      rate_T1E2 = dNg_over_dt_q[time_num][temp_num][HQenergy_num+1];
      rate_T2E1 = dNg_over_dt_q[time_num][temp_num+1][HQenergy_num];
      rate_T2E2 = dNg_over_dt_q[time_num][temp_num+1][HQenergy_num+1];
      max_T1E1 = max_dNgfnc_q[time_num][temp_num][HQenergy_num];
      max_T1E2 = max_dNgfnc_q[time_num][temp_num][HQenergy_num+1];
      max_T2E1 = max_dNgfnc_q[time_num][temp_num+1][HQenergy_num];
      max_T2E2 = max_dNgfnc_q[time_num][temp_num+1][HQenergy_num+1];
  } 

  rate_EGrid_low = rate_T1E1+(temp_med-temp_min-temp_num*delta_temp)/delta_temp*(rate_T2E1-rate_T1E1);
  rate_EGrid_high = rate_T1E2+(temp_med-temp_min-temp_num*delta_temp)/delta_temp*(rate_T2E2-rate_T1E2);
  max_EGrid_low = max_T1E1+(temp_med-temp_min-temp_num*delta_temp)/delta_temp*(max_T2E1-max_T1E1);
  max_EGrid_high = max_T1E2+(temp_med-temp_min-temp_num*delta_temp)/delta_temp*(max_T2E2-max_T1E2);

  delta_Ng = rate_EGrid_low+(HQenergy-HQenergy_num*delta_HQener)/delta_HQener*(rate_EGrid_high-rate_EGrid_low);
  max_Ng = max_EGrid_low+(HQenergy-HQenergy_num*delta_HQener)/delta_HQener*(max_EGrid_high-max_EGrid_low);
 
  delta_Ng=delta_Ng*6.0/D2piT*dtLRF;
  max_Ng=max_Ng*6.0/D2piT;

//  if(delta_Ng>1) {
//     cout << "Warning: Ng greater than 1   " << time_gluon << "  " << delta_Ng << endl;
//  }

  return delta_Ng;

//  write(6,*) temp_num,HQenergy_num,delta_Ng

}	  

//
  double Mqc2qc(double s, double t, double M) {

      double m2m=M*M;
      double u=2.0*m2m-s-t;
      double MM;

      MM=64.0/9.0*(pow((m2m-u),2)+pow((s-m2m),2)+2.0*m2m*t)/t/t;

      return(MM);

  }

  double Mgc2gc(double s, double t, double M) {

      double m2m=M*M;
      double u=2.0*m2m-s-t;
      double MM;

      MM=32.0*(s-m2m)*(m2m-u)/t/t;
      MM=MM+64.0/9.0*((s-m2m)*(m2m-u)+2.0*m2m*(s+m2m))/pow((s-m2m),2);
      MM=MM+64.0/9.0*((s-m2m)*(m2m-u)+2.0*m2m*(u+m2m))/pow((u-m2m),2);
      MM=MM+16.0/9.0*m2m*(4.0*m2m-t)/((s-m2m)*(m2m-u));
      MM=MM+16.0*((s-m2m)*(m2m-u)+m2m*(s-u))/(t*(s-m2m));
      MM=MM+16.0*((s-m2m)*(m2m-u)-m2m*(s-u))/(t*(u-m2m));

      return(MM);

  }

////////////////////////////////////////////////////////////////////////////////////////
// functions for gluon radiation from heavy quark

  double alphasHQ(double kTFnc, double tempFnc) {

//      double kTEff,error_para,resultFnc;
//
//      error_para=1.0;
//
//      if(kTFnc<pi*tempFnc*error_para) {
//         kTEff=pi*tempFnc*error_para;
//      } else {
//         kTEff=kTFnc;
//      }
//
//      resultFnc=4.0*pi/(11.0-2.0*nflavor(kTEff)/3.0)/2.0/log(kTEff/lambdas(kTEff));
//
//      return(resultFnc);
      return(alphas);
  }

////////////////////////////////////////////////////////////////////////////////////////

  double nflavor(double kTFnc) {

      double resultFnc;
      double cMass=1.27;
      double bMass=4.19;

      if (kTFnc<cMass) {
         resultFnc=3.0;
      } else if (kTFnc<bMass) {
         resultFnc=4.0;
      } else {
         resultFnc=5.0;
      }

      return(resultFnc);
  }

////////////////////////////////////////////////////////////////////////////////////////

  double lambdas(double kTFnc) {

      double resultFnc;
      double cMass=1.27;
      double bMass=4.19;

      if (kTFnc<cMass) {
         resultFnc=0.2;
      } else if (kTFnc<bMass) {
         resultFnc=0.172508;
      } else {
         resultFnc=0.130719;
      }

      return(resultFnc);
  }

////////////////////////////////////////////////////////////////////////////////////////

  double splittingP(int parID, double z0g) {

      double resultFnc;

//      resultFnc = (2.0-2.0*z0g+z0g*z0g)*CF/z0g;

      if(parID==21) resultFnc = 2.0*pow(1.0-z0g+pow(z0g,2),3)/z0g/(1.0-z0g);
      else resultFnc = (1.0-z0g)*(2.0-2.0*z0g+z0g*z0g)/z0g;

      return(resultFnc);
  }

////////////////////////////////////////////////////////////////////////////////////////

  double tau_f(double x0g, double y0g, double HQenergy, double HQmass) {

      double resultFnc;

      resultFnc = 2.0*HQenergy*x0g*(1.0-x0g)/(pow((x0g*y0g*HQenergy),2)+x0g*x0g*HQmass*HQmass);

      return(resultFnc);
  }

////////////////////////////////////////////////////////////////////////////////////////


  double dNg_over_dxdydt(int parID, double x0g, double y0g, double HQenergy, double HQmass, double temp_med, double Tdiff) {

      double resultFnc,tauFnc,qhatFnc;


      qhatFnc = qhat_over_T3*pow(temp_med,3);   // no longer have CF factor, taken out in splittingP too.
      tauFnc = tau_f(x0g,y0g,HQenergy,HQmass);

      resultFnc = 4.0/pi*CA*alphasHQ(x0g*y0g*HQenergy,temp_med)*splittingP(parID,x0g)*qhatFnc*pow(y0g,5)*pow(sin(Tdiff/2.0/tauFnc/sctr),2)*pow((HQenergy*HQenergy/(y0g*y0g*HQenergy*HQenergy+HQmass*HQmass)),4)/x0g/x0g/HQenergy/HQenergy/sctr;

      return(resultFnc);
  }

/////////////////////////////////////////////////////////////////////////////////////////

  void read_xyMC(int &numXY) {

      numXY=0;

      ifstream fxyMC("hydroProfile/mc_glauber.dat");
      if(!fxyMC.is_open()) {
          cout<<"Erro openning date fxyMC!\n";
          exit (EXIT_FAILURE);
      }

      while(true) {
          if(fxyMC.eof()) {
              numXY--;
              break;
          }
          if(numXY>=maxMC) break;
          fxyMC>>initMCX[numXY]>>initMCY[numXY];
          numXY++;
      }

      cout << "Number of (x,y) points from MC-Glauber: " << numXY << endl;

      fxyMC.close();

  }

  void jetClean(){
                
      maxColor=0; // SC:color

      for(int i=0;i<dimParList;i++) { // clear arrays before initialization for each event
          P[1][i]=0.0;
          P[2][i]=0.0;
          P[3][i]=0.0;
          P[0][i]=0.0;
          //
          P0[1][i]=0.0;
          P0[2][i]=0.0;
          P0[3][i]=0.0;
          P0[0][i]=0.0;
          //
          PP[1][i]=0.0;
          PP[2][i]=0.0;
          PP[3][i]=0.0;
          PP[0][i]=0.0;
          //
          PP0[1][i]=0.0;
          PP0[2][i]=0.0;
          PP0[3][i]=0.0;
          PP0[0][i]=0.0;
          //
          V[1][i]=0.0;
          V[2][i]=0.0;
          V[3][i]=0.0;
          //
          V0[1][i]=0.0;
          V0[2][i]=0.0;
          V0[3][i]=0.0;
          //
          VV[1][i]=0.0;
          VV[2][i]=0.0;
          VV[3][i]=0.0;
          //
          VV0[1][i]=0.0;
          VV0[2][i]=0.0;
          VV0[3][i]=0.0;
          //
          CAT[i]=0;
          CAT0[i]=0;
          //
          radng[i]=0.0;
          tirad[i]=0.0;
          tiscatter[i]=0.0;
          //
          Vfrozen[0][i]=0.0; 
          Vfrozen[1][i]=0.0; 
          Vfrozen[2][i]=0.0; 
          Vfrozen[3][i]=0.0; 
          Vfrozen0[0][i]=0.0;   
          Vfrozen0[1][i]=0.0;   
          Vfrozen0[2][i]=0.0;   
          Vfrozen0[3][i]=0.0;   
          Tfrozen[i]=0.0;
          Tfrozen0[i]=0.0;
          vcfrozen[0][i]=0.0;
          vcfrozen[1][i]=0.0;
          vcfrozen[2][i]=0.0;
          vcfrozen[3][i]=0.0;
          vcfrozen0[0][i]=0.0;
          vcfrozen0[1][i]=0.0;
          vcfrozen0[2][i]=0.0;
          vcfrozen0[3][i]=0.0;

          Tint_lrf[i]=0.0; //for heavy quark

          color[i]=0; // SC:color
          anticolor[i]=0;
          color0[i]=0;
          anticolor0[i]=0;
          fake[i]=0;
          fake0[i]=0;


      }		
  }

  void jetInitialize(int numXY) {

      double ipx,ipy,ipz,ip0,iwt;
      double pT_len,ipT,phi,rapidity;

//...single jet parton input, only useful when fixMomentum=1
      double px0=0.0;
      double py0=0.0;  // initial momemtum
      double pz0=sqrt(ener*ener-amss*amss);
      double en0=ener;
     
      for(int i=1; i<=nj; i=i+2) {

          if(fixPosition == 1) { // initialize position
              V[1][i]=0.0; 
              V[2][i]=0.0;
              V[3][i]=0.0;
          } else {
              int index_xy = (int)(ran0(&NUM1)*numXY);
              if(index_xy>=numXY) index_xy = numXY-1;
              V[1][i]=initMCX[index_xy]; 
              V[2][i]=initMCY[index_xy];
              V[3][i]=0.0;
          }
          V[0][i]=-log(1.0-ran0(&NUM1));

          if(fixMomentum == 1) { // initialize momentum
	      P[1][i]=px0;
	      P[2][i]=py0;
	      P[3][i]=pz0;
	      P[0][i]=en0;
	      P[4][i]=amss;
              P[5][i]=sqrt(px0*px0+py0*py0);
              WT[i]=1.0;
          } else {
              pT_len=ipTmax-ipTmin;
              ipT=ipTmin+ran0(&NUM1)*pT_len;
              phi=ran0(&NUM1)*2.0*pi;
              ipx=ipT*cos(phi);
              ipy=ipT*sin(phi);
              rapidity=2.0*eta_cut*ran0(&NUM1)-eta_cut;
              ipz=sqrt(ipT*ipT+amss*amss)*sinh(rapidity);
              ip0=sqrt(ipT*ipT+ipz*ipz+amss*amss);
              P[1][i]=ipx;
              P[2][i]=ipy;
              P[3][i]=ipz;
              P[0][i]=ip0;
              P[4][i]=amss;
              P[5][i]=ipT;
              WT[i]=pT_len;
          }

          KATT1[i]=Kjet;

          // let jet partons stream freely for tau_0 in x-y plane;
          V[1][i]=V[1][i]+P[1][i]/P[0][i]*tau0;
          V[2][i]=V[2][i]+P[2][i]/P[0][i]*tau0;
                                                  
          for(int j=0;j<=3;j++) Prad[j][i]=P[j][i];

          Vfrozen[0][i]=tau0;
          Vfrozen[1][i]=V[1][i];
          Vfrozen[2][i]=V[2][i];
          Vfrozen[3][i]=V[3][i];
          Tfrozen[i]=0.0;
          vcfrozen[1][i]=0.0;
          vcfrozen[2][i]=0.0;
          vcfrozen[3][i]=0.0;

          // now its anti-particle (if consider pair production)
	  if(KATT1[i]==21) KATT1[i+1]=21; else KATT1[i+1]=-KATT1[i];
          P[1][i+1]=-P[1][i];
          P[2][i+1]=-P[2][i];
          P[3][i+1]=-P[3][i];
          P[0][i+1]=P[0][i];
          P[4][i+1]=P[4][i];
          P[5][i+1]=P[5][i];
	  WT[i+1]=WT[i];
          V[1][i+1]=V[1][i];
          V[2][i+1]=V[2][i];
          V[3][i+1]=V[3][i];
          V[0][i+1]=V[0][i];
          for(int j=0;j<=3;j++) Prad[j][i+1]=P[j][i];
          Vfrozen[0][i+1]=Vfrozen[0][i];
          Vfrozen[1][i+1]=Vfrozen[1][i];
          Vfrozen[2][i+1]=Vfrozen[2][i];
          Vfrozen[3][i+1]=Vfrozen[3][i];
          Tfrozen[i+1]=Tfrozen[i];
          vcfrozen[1][i+1]=vcfrozen[1][i];
          vcfrozen[2][i+1]=vcfrozen[2][i];
          vcfrozen[3][i+1]=vcfrozen[3][i];
      } 
  }

  void setJetX(int numXY) {

      double setX,setY,setZ;
 
      if(fixPosition == 1) { 
          setX=0.0;
          setY=0.0;
          setZ=0.0;
      } else {
          int index_xy = (int)(ran0(&NUM1)*numXY);
          if(index_xy>=numXY) index_xy = numXY-1;
          setX=initMCX[index_xy]; 
          setY=initMCY[index_xy];
          setZ=0.0;
      }

      for(int i=1; i<=nj; i++) {

          V[1][i]=setX; 
          V[2][i]=setY;
          V[3][i]=setZ;
          V[0][i]=-log(1.0-ran0(&NUM1));

          V[1][i]=V[1][i]+P[1][i]/P[0][i]*tau0;
          V[2][i]=V[2][i]+P[2][i]/P[0][i]*tau0;
                                                  
          Vfrozen[0][i]=tau0;
          Vfrozen[1][i]=V[1][i];
          Vfrozen[2][i]=V[2][i];
          Vfrozen[3][i]=V[3][i];
          Tfrozen[i]=0.0;
          vcfrozen[1][i]=0.0;
          vcfrozen[2][i]=0.0;
          vcfrozen[3][i]=0.0;

      } 
  }



// ..........................................................................................................

    double calRapidity(double fourP[4]){

        double p0=fourP[0];
        double pz=fourP[3];
        double epsilon=1e-6;

        if(p0<fabs(pz)) {
            cout << "Energy smaller than pz when calculating rapidity!" << endl;
            exit(EXIT_FAILURE);
        }

        if(fabs(p0-pz)<epsilon) {
           if(pz>0) return 1000.0;
           else return -1000.0;
        } else {
           return 0.5*log((p0+pz)/(p0-pz));
        }
    }
