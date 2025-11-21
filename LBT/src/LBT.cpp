///////////////////////////////////////////////////////////////////////////////////////////////////
//
// This is the key function of LBT
// Update elastic and inelastic evolution of the particle list for one time step
//
///////////////////////////////////////////////////////////////////////////////////////////////////

void LBT(int &n, double &ti){
		
  int CT;                        //collision type 1 2 3 13 4 5 6 7 8
  int KATTC0;                    //flavor code 1/d 2/u 3/s -1/dbar -2/ubar -3/sbar 21/gluon 
  int KATT2;
  int KATT3;
  int KATT30;
		
  double RTE;                    //scattering rate (energy temperature)
  double E;                      //parton energy
  double PLen;                   //parton momentum
  double T;                      //local temperature		

  double T1;                     //index for difference method
  double T2;
  double E1;
  double E2;
  int iT1;
  int iT2;
  int iE1;
  int iE2;	
		
  int nrad;
  int idlead;
  int idlead1;
  int idlead2;
  double Xtau;                   //....................main & titau	
  double Vx;
  double Vy;		
  double Veta;

  double tcar;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar;		
  double ycar;
  double zcar;					

  double tcar0;                   //....................t x y z in Cartesian coordinate this is for the radiation process
  double xcar0;		
  double ycar0;
  double zcar0;

  
  double rans;
		
  int KATTx;
  int KATTy;

  double Ejp;
  double Elab;
		
  double qt;
		
  double px0;
  double py0;		
  double pz0;

  double Ncoll22;                    //...................average number of elastic scattering for particle i in the current step  
  int Nscatter;                   //...................number of elastic scattering for particle i in the current step

  int nnpp=np;                    //...................number of particles in the current np loop
  int np0=np;                     //...................number of particles in the current step (np0)
                                  //...................number of particles in the beginning of current step (np)

  int free=0;
  int free0=0;
  double fraction=0.0;
  double fraction0=0.0;
  double vc0b[4]={0.0};             //flow velocity     
  double pMag,vMag,flowFactor;

  double kt2;
  
  double vp0[4]={0.0};
  double p0temp[4]={0.0}; 

  double p0temp1=0.0; 
  double p0temp2=0.0; 
		
  double pcx[4]={0.0};
  double pcy[4]={0.0};		

  double pcx1[4]={0.0};
  double pcy1[4]={0.0};		

// for heavy quark
  double dt_lrf,maxFncHQ,Tdiff,lim_low,lim_high,lim_int;
  double probCol,probRad,probTot;

//...........................................................................................................NCUT
  ncut=0;
  ncut0=0;
  //...........................................................................................................NCUTEND		  
		  
  //...collision of the jet parton (i<=nj), recoiled partons, radiated gluons with the background thermal partons
  //...thermal partons
  //...np number of partons in current step
  //...nj number of jet partons 

  for(int i=1;i<=np;i++) {
      if(KATT1[i]!=21 && abs(KATT1[i])>5) continue; 
      // if the production time of a parton is ahead of the current time, do nothing
      if(Vfrozen[0][i]>=ti) continue;  

      icl22=0;
      qt=0;
      idlead=i;
		  
//......propagation of each positive particle and its negative counterpart in the current time step (for all particles)		
		
      if(P[0][i]<epsilon) {  //anti-overflow NOT NECESSARY
          V[1][i]=0.0;
          V[2][i]=0.0;
          V[3][i]=0.0;
          CAT[i]=1;		  
      }
      
      if(P0[0][i]<epsilon) {  //anti-overflow NOT NECESSARY
          V0[1][i]=0.0;
          V0[2][i]=0.0;
          V0[3][i]=0.0;
          CAT0[i]=1;
      }
 
      if(CAT0[i] != 1) {	  
//...........................................t-z coordinates	
          if(tauswitch==0) {

              V0[1][i]=V0[1][i]+dt*P0[1][i]/P0[0][i];
              V0[2][i]=V0[2][i]+dt*P0[2][i]/P0[0][i];
	      V0[3][i]=V0[3][i]+dt*P0[3][i]/P0[0][i];

	      tcar0=ti;
	      zcar0=V0[3][i];
	      xcar0=V0[1][i];
	      ycar0=V0[2][i];

              double ed00,sd00,VX00,VY00,VZ00;
              int hydro_ctl0;

	      free0=0;

              if(bulkFlag==1) { // read OSU hydro
                  readhydroinfoshanshan_(&tcar0,&xcar0,&ycar0,&zcar0,&ed00,&sd00,&temp00,&VX00,&VY00,&VZ00,&hydro_ctl0);
              } else if(bulkFlag==2) { // read CCNU hydro
                  hydroinfoccnu_(&tcar0, &xcar0, &ycar0, &zcar0, &temp00, &VX00, &VY00, &VZ00, &hydro_ctl0);
              } else if(bulkFlag==0) { // static medium
                  VX00=0.0;
                  VY00=0.0;
                  VZ00=0.0;
                  hydro_ctl0=0;
              }

              if(hydro_ctl0==0 && temp00>=hydro_Tc) {

	          qhat00=DebyeMass2(Kqhat0,alphas,temp00);
                  fraction0=1.0;

                  Vfrozen0[0][i]=ti;
                  Vfrozen0[1][i]=V0[1][i];
                  Vfrozen0[2][i]=V0[2][i];
                  Vfrozen0[3][i]=V0[3][i];
                  Tfrozen0[i]=temp00;
                  vcfrozen0[1][i]=VX00;
                  vcfrozen0[2][i]=VY00;
                  vcfrozen0[3][i]=VZ00;

              } else {
                  free0=1;
                  fraction0=0.0;
              }

	  }//if(tauswitch==0)
		  
      }//if(CAT0[i] != 1)

      if(CAT[i]!=1) {
//...........................................t-z coordinates	
          if(tauswitch==0) {

              V[1][i]=V[1][i]+dt*P[1][i]/P[0][i];
              V[2][i]=V[2][i]+dt*P[2][i]/P[0][i];
              V[3][i]=V[3][i]+dt*P[3][i]/P[0][i];

	      tcar=ti;
	      zcar=V[3][i];
	      xcar=V[1][i];
	      ycar=V[2][i];

          //...........................................Hydro part		

              double ed,sd,VX,VY,VZ;
              int hydro_ctl;

	      free=0;

              // SC: do not process particles outside the preset rapidity range 
              double fourP[4]={P[0][i],P[1][i],P[2][i],P[3][i]};
              double rapidity=calRapidity(fourP);
              // cout << "check rapidity: " << P[0][i] << "    " << P[1][i] << "    " << P[2][i] << "    " << P[3][i] << "    " << rapidity << endl;
              if(fabs(rapidity) > cutRapidity) free=1;

              if(free==0) {

                  if(bulkFlag==1) { // read OSU hydro
                      readhydroinfoshanshan_(&tcar,&xcar,&ycar,&zcar,&ed,&sd,&temp0,&VX,&VY,&VZ,&hydro_ctl);
                  } else if(bulkFlag==2) { // read CCNU hydro
                      hydroinfoccnu_(&tcar, &xcar, &ycar, &zcar, &temp0, &VX, &VY, &VZ, &hydro_ctl);
                  } else if(bulkFlag==0) { // static medium
                      VX=0.0;
                      VY=0.0;
                      VZ=0.0;
                      hydro_ctl=0;
                  }

                  if(hydro_ctl==0 && temp0>=hydro_Tc) {

                      vc0[1]=VX;
                      vc0[2]=VY;
                      vc0[3]=VZ;

//                      vc0[1]=0.0;
//                      vc0[2]=0.0;
//                      vc0[3]=0.0;

                      vc0b[1]=VX;
                      vc0b[2]=VY;
                      vc0b[3]=VZ;

        
                      //...alphas!	
                      alphas=alphas0(Kalphas,temp0);
                      //...Debye Mass square
                      qhat0=DebyeMass2(Kqhat0,alphas,temp0);	
	
                      fraction=1.0;
                      Vfrozen[0][i]=ti;
                      Vfrozen[1][i]=V[1][i];
                      Vfrozen[2][i]=V[2][i];
                      Vfrozen[3][i]=V[3][i];
                      Tfrozen[i]=temp0;
                      vcfrozen[1][i]=VX;
                      vcfrozen[2][i]=VY;
                      vcfrozen[3][i]=VZ;

                  } else {
                      free=1;
                      fraction=0.0;
                  }
  
       	          pc0[1]=P[1][i];
                  pc0[2]=P[2][i];
       	          pc0[3]=P[3][i];
       	          pc0[0]=P[0][i];

                  vMag=sqrt(vc0b[1]*vc0b[1]+vc0b[2]*vc0b[2]+vc0b[3]*vc0b[3]);
                  flowFactor=(1.0-(pc0[1]*vc0b[1]+pc0[2]*vc0b[2]+pc0[3]*vc0b[3])/pc0[0])/sqrt(1.0-vMag*vMag);
 
	          trans(vc0,pc0);
			  
	          if(pc0[0]<sqrt(qhat0)) {
	              free=1;
	          }
	  
	          if(i>nj && pc0[0]<Ecmcut) {
	              free=1;
                  }   
		  
	          transback(vc0,pc0);

              }//if(free==0)
	  } //if(tauswitch==0)
				

//...........................................tau-eta coordinates		
	  if(tauswitch==1) {
		
		//????????????? why I need to get the medium information before propagation
		
          //...........................................Hydro part		

		  /////////////////////////////////////////////////////////////
		  //           GetInfo( CD & tau, CD & x, CD &y, CD &etas, double & Ed, double & T, double & Vx, double & Vy, double & Veta, double & phase_ratio )        
		  /////////////////////////////////////////////////////////////

          //...........................................Hydro part end
		  
		  vp[0]=0.0;
		  vp[1]=V[1][i];
		  vp[2]=V[2][i];
		  vp[3]=V[3][i];

		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];
		  
		  titau(ti,vc0,vp,pc0,Vx,Vy,Veta,Xtau);

		  V[1][i]=V[1][i]+dt*Vx;
		  V[2][i]=V[2][i]+dt*Vy;
		  V[3][i]=V[3][i]+dt*Veta;

		  ///////////////////////////////////////time and position in Cartesian coordinate			 
		  tcar=ti*cosh(V[3][i]);
		  zcar=ti*sinh(V[3][i]);
		  xcar=V[1][i];
		  ycar=V[2][i];

                  alphas=alphas0(Kalphas,temp0);
		  qhat0=DebyeMass2(Kqhat0,alphas,temp0);


		  if(i>nj)
		  {
		  //......................................negative particle propagation
		  
		  /////////////////////////////////////////////////////////////
		  //           GetInfo( CD & tau, CD & x, CD &y, CD &etas, double & Ed, double & T, double & Vx, double & Vy, double & Veta, double & phase_ratio )        
		  /////////////////////////////////////////////////////////////
		  
		  vp0[0]=0.0;
		  vp0[1]=V0[1][i];
		  vp0[2]=V0[2][i];
		  vp0[3]=V0[3][i];

		  pc00[1]=P0[1][i];
		  pc00[2]=P0[2][i];
		  pc00[3]=P0[3][i];
		  pc00[0]=P0[0][i];
		  
		  titau(ti,vc0,vp,pc00,Vx,Vy,Veta,Xtau);

		  V0[1][i]=V0[1][i]+dt*Vx;
		  V0[2][i]=V0[2][i]+dt*Vy;
		  V0[3][i]=V0[3][i]+dt*Veta;

		  ///////////////////////////////////////time system in Cartesian coordinate			 
		  tcar=ti*cosh(V0[3][i]);
		  zcar=ti*sinh(V0[3][i]);
		  xcar=V0[1][i];
		  ycar=V0[2][i];
			
//                temp00= temp0;
			
//		  qhat00=DebyeMass2(Kqhat0,alphas,temp00);
		  
		  } // if(i>nj)
          } // if(tauswitch==1)			

//......propagation of each positive particle and its negative counterpart in the current time step	end	

//......CAT[i]=1 free streaming particle

//......free streaming when energy  < Ecut
//......free streaming when energy0 < Ecut0......in the next step

          if(free==0) {

              double qhatTP,RTE1,RTE2;

// modification: interplotate rate in l.r.f. 
       	      pc0[1]=P[1][i];
       	      pc0[2]=P[2][i];
       	      pc0[3]=P[3][i];
       	      pc0[0]=P[0][i];
              trans(vc0,pc0);
              E=pc0[0]; //  p4-the initial 4-momentum of the jet parton in the local rest frame
              PLen=sqrt(pc0[1]*pc0[1]+pc0[2]*pc0[2]+pc0[3]*pc0[3]);
              transback(vc0,pc0);

	      T=temp0;
	      KATTC0=KATT1[i];

//              if(E<T) preKT=4.0*pi/9.0/log(scaleAK*T*T/0.04)/alphas/fixKT;
//              else preKT=4.0*pi/9.0/log(scaleAK*E*T/0.04)/alphas/fixKT;

//              runKT=fixAlphas/0.3;
              runKT=4.0*pi/9.0/log(2.0*E*T/0.04)/0.3;
	      lam(KATTC0,RTE,PLen,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2); //modified: use P instead

//              preKT=alphas/0.3;

// calculate p,T-dependence K factor
//              KPfactor=1.0;
//              KTfactor=1.0;
//              if(fabs(KATT1[i])==4||fabs(KATT1[i])==5) {
                  KPfactor=1.0+KPamp*exp(-PLen*PLen/2.0/KPsig/KPsig);
                  KTfactor=1.0+KTamp*exp(-pow((temp0-hydro_Tc),2)/2.0/KTsig/KTsig);
//              }
//              Kfactor=KPfactor*KTfactor*KTfactor*preKT*preKT;      
              Kfactor=KPfactor*KTfactor*KTfactor*runKT*preKT;      

// get qhat from table
              if(KATTC0==21) {
                  RTE1=(qhatG[iT2][iE1]-qhatG[iT1][iE1])*(T-T1)/(T2-T1)+qhatG[iT1][iE1];
                  RTE2=(qhatG[iT2][iE2]-qhatG[iT1][iE2])*(T-T1)/(T2-T1)+qhatG[iT1][iE2];
              } else if(KATTC0==4||KATTC0==-4) {
                  RTE1=(qhatHQ[iT2][iE1]-qhatHQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE1];
                  RTE2=(qhatHQ[iT2][iE2]-qhatHQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatHQ[iT1][iE2];
              } else if(KATTC0==5||KATTC0==-5) {
                  RTE1=(qhatHQb[iT2][iE1]-qhatHQb[iT1][iE1])*(T-T1)/(T2-T1)+qhatHQb[iT1][iE1];
                  RTE2=(qhatHQb[iT2][iE2]-qhatHQb[iT1][iE2])*(T-T1)/(T2-T1)+qhatHQb[iT1][iE2];
              } else {
                  RTE1=(qhatLQ[iT2][iE1]-qhatLQ[iT1][iE1])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE1];
                  RTE2=(qhatLQ[iT2][iE2]-qhatLQ[iT1][iE2])*(T-T1)/(T2-T1)+qhatLQ[iT1][iE2];
              }

              qhatTP=(RTE2-RTE1)*(PLen-E1)/(E2-E1)+RTE1;

              qhatTP=qhatTP*Kfactor;

              qhat_over_T3=qhatTP;  // what is read in is qhat/T^3 of quark/gluon

              // D2piT is diffusion coefficient "just for quark"
              if(KATTC0==21) D2piT=8.0*pi/(qhat_over_T3/CA*CF); else D2piT=8.0*pi/qhat_over_T3;

              qhat=qhatTP*pow(T,3); // for light quark and gluon, need additional table to make it correct

// update interaction time and average gluon number for heavy quark radiation 
// (as long as it's in medium, even though no collision)
//              dt_lrf=dt*sqrt(1.0-vc0[1]*vc0[1]-vc0[2]*vc0[2]-vc0[3]*vc0[3]);
              dt_lrf=dt*flowFactor;
              Tint_lrf[i]=Tint_lrf[i]+dt_lrf;
              Tdiff=Tint_lrf[i];

// get radiation probablity by reading tables -- same for heavy and light partons
              if(KATTC0==21) radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*KTfactor/2.0*runKT;
              else radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*KTfactor*runKT;
//              if(KATTC0==21) radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)/2.0*KTfactor*preKT;
//              else radng[i]+=nHQgluon(KATT1[i],dt_lrf,Tdiff,temp0,E,maxFncHQ)*KTfactor*preKT;
              lim_low=sqrt(6.0*pi*alphas)*temp0/E;
//              lim_low=sqrt(6.0*pi*0.3)*temp0/E;
              if(abs(KATT1[i])==4 || abs(KATT1[i])==5) lim_high=1.0;
              else lim_high=1.0-lim_low;
              lim_int=lim_high-lim_low;
              if(lim_int>0.0) probRad=1.0-exp(-radng[i]);
              else probRad=0.0;

//...........................................t-z coordinates
              if(tauswitch==0) {

		  pc0[1]=P[1][i];
		  pc0[2]=P[2][i];
		  pc0[3]=P[3][i];
		  pc0[0]=P[0][i];
		  
//                V[0][i]=V[0][i]-dt*RTE/0.1970*sqrt(pow(P[1][i],2)+pow(P[2][i],2)+pow(P[3][i],2))/P[0][i];	
//		  Ncoll22=dt*RTE/0.197*((pc0[0]*1-pc0[1]*vc0[1]-pc0[2]*vc0[2]-pc0[3]*vc0[3])/E);		  

//??????@@@
//////////////////
                  V[0][i]=V[0][i]-fraction*dt*RTE/0.1970*sqrt(pow(P[1][i],2)+pow(P[2][i],2)+pow(P[3][i],2))/P[0][i];
                  probCol=fraction*dt_lrf*RTE/0.1970;
              }
//...........................................tau-eta coordinates
              if(tauswitch==1) {

                  V[0][i]=V[0][i]-dt*RTE*Xtau/0.1970*sqrt(pow(P[1][i],2)+pow(P[2][i],2)+pow(P[3][i],2))/P[0][i];
                  probCol=dt_lrf*RTE*Xtau/0.1970;  // ?? fraction
//		  Ncoll22=dt*RTE/0.197*Xtau*((pc0[0]*1-pc0[1]*vc0[1]-pc0[2]*vc0[2]-pc0[3]*vc0[3])/E);  		  
              }


//          cout<<"t, T, PLen, rate, lambda:  "<<ti<<"  "<<T<<"  "<<PLen<<"  "<<RTE<<"  "<<1/RTE*0.197<<endl;
////////////////////////////////////////...new elastic part	  
 
//          if(ran0(&NUM1)<1-exp(-Ncoll22)) // !Yes, collision!
////////////////////////////////////////		  

              if(KINT0==0) probRad=0.0;  // switch off radiation
//              probCol=probCol*KPfactor*KTfactor*preKT;
              probCol=probCol*KPfactor*KTfactor*runKT;
              probCol=(1.0-exp(-probCol))*(1.0-probRad);  // probability of pure elastic scattering
              if(KINT0==2) probCol=0.0;
              probTot=probCol+probRad;

	              		  
              if(ran0(&NUM1)<probTot) { // !Yes, collision! Either elastic or inelastic.

//			  Nscatter=KPoisson(Ncoll22);

                  Nscatter=1; 
			  
                  for(int nsca=1;nsca<=Nscatter;nsca++) {			
			
//                          cout << "scatter!" << endl;
//			  np0=np0+1;
			  pc0[1]=P[1][i];
			  pc0[2]=P[2][i];
			  pc0[3]=P[3][i];
			  pc0[0]=P[0][i];

// use PLen instead			  
			  flavor(CT,KATTC0,KATT2,KATT3,RTE,PLen,T,T1,T2,E1,E2,iT1,iT2,iE1,iE2);	
			 
			  //...........collision between a jet parton or a recoiled parton and a thermal parton from the medium 			 

                          if(CT==11||CT==12) { // for heavy quark scattering
			      collHQ22(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt,KATT1[i]);
                          } else { // for light parton scattering
			      colljet22(CT,temp0,qhat0,vc0,pc0,pc2,pc3,pc4,qt);
                          }

                          if(qt>epsilon) { //SC: update 2024-12-18, do scattering only when qt is properly sampled

                              np0=np0+1; // valid scattering

                              icl22 = 1;
                              n_coll22+=1;
			 
		    	      KATT1[i]=KATTC0;
			      KATT1[np0]=KATT2;
			      KATT10[np0]=KATT3;

			      if(pc0[0]<pc2[0] && abs(KATTC0)!=4 && abs(KATTC0)!=5 && KATTC0==KATT2) { //disable switch for heavy quark, only allow switch for identical particles
			          for(int k=0;k<=3;k++) {
				      p0temp[k]=pc2[k];
				      pc2[k]=pc0[k];
				      pc0[k]=p0temp[k];
			          }
			          KATT1[i]=KATT2;
     			          KATT1[np0]=KATTC0;
			      }

                              // SC:color assign color to jet parton 0, recoil 2 and negative 3 
                              if(colorFlag==1) {
                                  int parentColor=color[i];  // color of the incoming jet parton
                                  int parentAntiColor=anticolor[i];
                                  switch(CT) {
                                      case 1: // gg->gg
                                          maxColor++;
                                          color[i]=parentColor;
                                          anticolor[i]=maxColor;
                                          color[np0]=maxColor;
                                          maxColor++;
                                          anticolor[np0]=maxColor;
                                          color0[np0]=parentAntiColor;
                                          anticolor0[np0]=maxColor;
                                          break;
                                      case 2: // gg->qqbar
                                          if(KATT1[i]>0) { // final jet parton is quark
                                              color[i]=parentColor;
                                              anticolor[i]=0;
                                              maxColor++;
                                              color[np0]=0;
                                              anticolor[np0]=maxColor;
                                              color0[np0]=parentAntiColor;
                                              anticolor0[np0]=maxColor;
                                          } else { // final jet parton is antiquark
                                              color[i]=0;
                                              anticolor[i]=parentAntiColor;
                                              maxColor++;
                                              color[np0]=maxColor;
                                              anticolor[np0]=0;
                                              color0[np0]=maxColor;
                                              anticolor0[np0]=parentColor;
                                          }
                                          break;
                                      case 3: // gq->qq or gqbar->gqbar
                                          if(KATT1[np0]>0) { // gq->gq
                                              maxColor++;
                                              color[i]=parentColor;
                                              anticolor[i]=maxColor;
                                              color[np0]=maxColor;
                                              anticolor[np0]=0;
                                              color0[np0]=parentAntiColor;
                                              anticolor0[np0]=0;
                                          } else { // gqbar->gqbar
                                              maxColor++;
                                              color[i]=maxColor;
                                              anticolor[i]=parentAntiColor;
                                              color[np0]=0;
                                              anticolor[np0]=maxColor;
                                              color0[np0]=0;
                                              anticolor0[np0]=parentColor;
                                          }
                                          break;
                                      case 13: // qg->qg or qbarg->qbarg
                                      case 12: // Qg->Qg or Qbarg->Qbarg
                                          if(KATT1[i]>0) { // qg->qg
                                              maxColor++;
                                              color[i]=maxColor;
                                              anticolor[i]=0;
                                              anticolor[np0]=maxColor;
                                              maxColor++;
                                              color[np0]=maxColor;
                                              color0[np0]=maxColor;
                                              anticolor0[np0]=parentColor;
                                          } else { // qbarg->qbarg
                                              maxColor++;
                                              color[i]=0;
                                              anticolor[i]=maxColor;
                                              color[np0]=maxColor;
                                              maxColor++;
                                              anticolor[np0]=maxColor;
                                              color0[np0]=parentAntiColor;
                                              anticolor0[np0]=maxColor;
                                          }
                                          break;
                                      case 4: // scatterings between different-flavor light quarks (antiquarks)
                                      case 5: // scatterings between same-flavor light quarks (antiquarks) -- only use the first and fourth if's below
                                      case 7: // qqbar->qqbar same-flavor paritcle-anti-particle scattering -- only use the second and third if's below, some part of this channel should be moved to annihilation process in the future
                                      case 11: // scatterings between heavy and light quarks (antiquarks)    
                                          if(KATT1[i]>0 && KATT1[np0]>0) { // qq'->qq'
                                              maxColor++;
                                              color[i]=maxColor;
                                              anticolor[i]=0;
                                              color[np0]=parentColor;
                                              anticolor[np0]=0;
                                              color0[np0]=maxColor;
                                              anticolor0[np0]=0;
                                          } else if(KATT1[i]>0 && KATT1[np0]<0) { // qq'bar->qq'bar
                                              maxColor++;
                                              color[i]=maxColor;
                                              anticolor[i]=0;
                                              color[np0]=0;
                                              anticolor[np0]=maxColor;
                                              color0[np0]=0;
                                              anticolor0[np0]=parentColor;
                                          } else if(KATT1[i]<0 && KATT1[np0]>0) { // qbarq'->qbarq'
                                              maxColor++;
                                              color[i]=0;
                                              anticolor[i]=maxColor;
                                              color[np0]=maxColor;
                                              anticolor[np0]=0;
                                              color0[np0]=parentAntiColor;
                                              anticolor0[np0]=0;
                                          } else { // qbarq'bar->qbarq'bar
                                              maxColor++;
                                              color[i]=0;
                                              anticolor[i]=maxColor;
                                              color[np0]=0;
                                              anticolor[np0]=parentAntiColor;
                                              color0[np0]=0;
                                              anticolor0[np0]=maxColor;
                                          }
                                          break;
                                      case 6: // qqbar annihilation to q'q'bar 
                                          //if(KATT1[i]>0 && KATT1[np0]>0) { // +- -> +-
                                          if(parentColor>0 && KATT1[i]>0) { // +- -> +-
                                              maxColor++;
                                              color[i]=parentColor;
                                              anticolor[i]=0;
                                              color[np0]=0;
                                              anticolor[np0]=maxColor;
                                              color0[np0]=0;
                                              anticolor0[np0]=maxColor;
                                          //} else if(KATT1[i]>0 && KATT1[np0]<0) { // +- -> -+
                                          } else if(parentColor>0 && KATT1[i]<0) { // +- -> -+
                                              maxColor++;
                                              color[i]=0;
                                              anticolor[i]=maxColor;
                                              color[np0]=parentColor;
                                              anticolor[np0]=0;
                                              color0[np0]=0;
                                              anticolor0[np0]=maxColor;
                                          //} else if(KATT1[i]<0 && KATT1[np0]>0) { // -+ -> +-
                                          } else if(parentColor==0 && KATT1[i]>0) { // -+ -> +-
                                              maxColor++;
                                              color[i]=maxColor;
                                              anticolor[i]=0;
                                              color[np0]=0;
                                              anticolor[np0]=parentAntiColor;
                                              color0[np0]=maxColor;
                                              anticolor0[np0]=0;
                                          } else { // -+ -> -+
                                              maxColor++;
                                              color[i]=0;
                                              anticolor[i]=parentAntiColor;
                                              color[np0]=maxColor;
                                              anticolor[np0]=0;
                                              color0[np0]=maxColor;
                                              anticolor0[np0]=0;
                                          }
                                          break;
                                      case 8: // qqbar annihilation to gg                          
                                          //if(KATT1[i]>0) { // incoming jet parton is quark
                                          if(parentColor>0) { // incoming jet parton is quark
                                              maxColor++;
                                              color[i]=parentColor;
                                              anticolor[i]=maxColor;
                                              color[np0]=maxColor;
                                              maxColor++;
                                              anticolor[np0]=maxColor;
                                              color0[np0]=0;
                                              anticolor0[np0]=maxColor;
                                          } else { // incoming jet parton is antiquark
                                              maxColor++;
                                              color[i]=maxColor;
                                              anticolor[i]=parentAntiColor;
                                              anticolor[np0]=maxColor;
                                              maxColor++;
                                              color[np0]=maxColor;
                                              color0[np0]=maxColor;
                                              anticolor0[np0]=0;
                                          }
                                          break;
                                      default:
                                          cout << "Waning: unrecognized channel number for elastic scattering!" << endl;
                                          exit(EXIT_FAILURE);
                                  } // end switch for selecting scattering channel
                                  //cout << CT << "  " << parentColor << "  " << parentAntiColor << "  " << color[i] << "  " << anticolor[i] << "  " << color[np0] << "  " << anticolor[np0] << "  " << color0[np0] << "  " << anticolor0[np0] << endl;
                                  //cout << "newly produced: " << np0 << "  " << color0[np0] << "  " << anticolor0[np0] << endl;
                                  //if(fabs(ti-11.3) < epsilon) cout << CT << "  " << np0 << "  " << KATT1[i] << "  " << KATT1[np0] << "  " << color[np0] << "  " << anticolor[np0] << endl;
                              } // end if(colorFlag==1)


			  for(int j=0;j<=3;j++) {
//			    if(i<=nj) // always make sure the jet parton carries the largest momentum (maybe now all the parton includes recoiled for there could be multiple scattering in a single time step)
				
				P[j][i]=pc0[j];
				P[j][np0]=pc2[j];
				V[j][np0]=V[j][i];
				
//		        if(i<=nj || nsca!=1) // then for each recoiled parton there is a negative parton
//				{
				P0[j][np0]=pc3[j];
				V0[j][np0]=V[j][i];

                                Vfrozen[j][np0]=V[j][i];
                                Vfrozen0[j][np0]=V[j][i];
                                Tfrozen[np0]=temp0;
                                Tfrozen0[np0]=temp0;

                                if(j!=0) {                                 
                                    vcfrozen[j][np0]=vc0b[j];
                                    vcfrozen0[j][np0]=vc0b[j];
                                }
//				}
			  }

                          // pass initial pT and weight information from jet parton to recoil and negative parton
                          P[5][np0]=P[5][i];
                          P0[5][np0]=P[5][i];
                          WT[np0]=WT[i];
                          WT0[np0]=WT[i];
                          P[4][np0]=0.0;  // recoiled parton is always massless
                          P0[4][np0]=0.0;


//CAT!!!			  
/*	       				
*/
                          CAT[np0]=2;
                          } // done update particle list for valid elastic scattering
//.............................................................. collision between "negative" partons

//...background thermal parton scattering with the same thermal parton
//...the associated parton has scattered

//...negative parton can only get the rank i>nj and it can only act on the first collision in the current time step. 
//...This is not physical only because the recoiled parton and its negative partner appear in pairs 

//...if Ecut0
//

//.............................................................. collision between "negative" partons

//...background thermal parton scattering with the same thermal parton the associated parton has scattered
//...negative parton can only get the rank i>nj and it can only act on the first collision in the current time step. 
//...This is not physical only because the recoiled parton and its negative partner appear in pairs 

////...Scattering between P0[i] and the newest negative parton  // commented out by Shanshan, just let negative particles stream freely
//
//            if(CAT0[i] == 0)
//            {
////???
////                free0=1;
////???                
//                if(free0 == 0)
//                {
//		    if(i>nj && nsca==1)
//		    {
//
//                        temp00=temp0;
//
//		        qhat00=DebyeMass2(Kqhat0,alphas,temp00);
//
//                        pc00[1]=P0[1][i];
//                        pc00[2]=P0[2][i];
//                        pc00[3]=P0[3][i];
//                        pc00[0]=P0[0][i];
//				
//			pc30[1]=pc3[1];
//			pc30[2]=pc3[2];
//			pc30[3]=pc3[3];
//			pc30[0]=pc3[0];
//
//	                E=pc00[0];
//	                T=temp00;
//	                KATTC0=KATT10[i];
//			KATT30=KATT3;
//        
//			twflavor(CT,KATTC0,KATT30,E,T);		 			
//					   				   
//                        twcoll(CT,qhat00,vc0,pc00,pc30);		   
//				   				   				   
//	                KATT1[np0]=KATT2;
//	                KATT10[i]=KATTC0;
//	                KATT10[np0]=KATT30;
//
//                        for(int j=0;j<=3;j++)
//			{
//                            P0[j][i]=pc00[j];
//                            P0[j][np0]=pc30[j];
//                            V0[j][np0]=V[j][i];
//                        }
//               	        V0[0][np0]=-log(1.0-ran0(&NUM1));
//		    } //if(i>nj && nsca==1)
//                } // if(free0==0)
//            } //.............................................................. "negative" end
											
             
	          }//for(int nsca=1;nsca<=Nscatter;nsca++)			  			  
			  
//...............................................................................................the point where the elastic scattering end			  

                  if(qt<epsilon) { 
                      KINT=0;
                  } else KINT=1;

                  if(qt>epsilon)
                  {
                        qperp[nqperp]=qt;
                        nqperp=nqperp+1;
                  }
        
                  if(KINT0!=0 && KINT!=0) { // Switch on the radiation processes
                      for(int j=0;j<=3;j++) {
                          p0temp[j]=P[j][i];
                      }							
                	
                      px0=pc4[1];
                      py0=pc4[2];
                      pz0=pc4[3];
        
        	      Elab=pc4[0]; // p4-the initial 4-momentum of the jet parton in the lab frame
                      trans(vc0,pc4);
                      Ejp=pc4[0]; //  p4-the initial 4-momentum of the jet parton in the local rest frame
                      transback(vc0,pc4);
        
                      if(Ejp>2*sqrt(qhat0)) { // !radiation or not?
        
        		  Vtemp[1]=xcar;
        		  Vtemp[2]=ycar;
        		  Vtemp[3]=zcar;
        
                          if(KATT1[i] != 21) {
                              isp=1;
                              n_sp1 +=1;
                          } else {
                              isp=2;
                              n_sp2 +=1;
                          }
         
                          if(ran0(&NUM1)<probRad/probTot) { // radiation -- either heavy or light parton:
        
        	              for(int j=0;j<=3;j++) { // prepare for the collHQ23
        			  pc01[j]=pc4[j]; // 2->3 process
        			  pb[j]=pc4[j];
                      	      }
        
        //		      colljet23(temp0,qhat0,vc0,pc01,pc2,pc3,pc4,qt,icl23,tcar,tiscatter[i],tirad[i],Elab,Ejp);
        		      collHQ23(KATT1[i],temp0,qhat0,vc0,pc01,pc2,pc3,pc4,qt,icl23,Tdiff,Ejp,maxFncHQ,lim_low,lim_int);
                              n_coll23+=1;
        
        		      if(icl23!=1) {
        
                                  int ctGluon=1;
        
                                  ng_coll23+=1; 
        			  nrad=KPoisson(radng[i]);          
                                  int nrad0=nrad;
        
        //                          nrad=1; // only allow single gluon emission right now
                                  
                                  ng_nrad+=nrad;
        			  np0=np0+1;
        			  KATT1[np0]=21;
//                                  CAT[np0]=4;
                                  // SC:color assign color infomation to splitting
                                  if(colorFlag==1) {
                                      int parentColor=color[i];  // taken from final state of 22
                                      int parentAntiColor=anticolor[i];
                                      if(KATT1[i]==21) { // g->gg process
                                          maxColor++;
                                          color[i]=parentColor;
                                          anticolor[i]=maxColor;
                                          color[np0]=maxColor;    // radiated gluon
                                          anticolor[np0]=parentAntiColor;
                                      } else if(KATT1[i]>0){ // q->qg process
                                          maxColor++;
                                          color[i]=maxColor;
                                          anticolor[i]=0;
                                          color[np0]=parentColor; // radiated gluon
                                          anticolor[np0]=maxColor;
                                      } else { // qbar->qbarg
                                          maxColor++;
                                          color[i]=0;
                                          anticolor[i]=maxColor;
                                          color[np0]=maxColor;
                                          anticolor[np0]=parentAntiColor;
                                      }
                                  } 
     			  
                                  for(int j=0;j<=3;j++) {
            			      P[j][i]=pc01[j];            // heavy quark after colljet23
        			      P[j][np0-1]=pc2[j];         // replace the final state of 22
        			      V[j][np0-1]=V[j][i];
        			      P0[j][np0-1]=pc3[j];
        			      V0[j][np0-1]=V[j][i];
        			      P[j][np0]=pc4[j];           //radiated gluon from colljet23
        			      V[j][np0]=V[j][i];
        			      P0[j][np0]=0.0;
        			      V0[j][np0]=0.0;
        
                                      Vfrozen[j][np0]=V[j][i];
                                      Vfrozen0[j][np0]=0.0;
                                      if(j!=0) {                                 
                                          vcfrozen[j][np0]=vc0b[j];
                                          vcfrozen0[j][np0]=0.0;
                                      }
        			  }
        
                                  Tfrozen[np0]=temp0;
                                  Tfrozen0[np0]=0.0;
        
                                  // pass initial pT and weight information
                                  P[5][np0]=P[5][i];
                                  P0[5][np0]=0.0;
                                  WT[np0]=WT[i];
                                  WT0[np0]=0.0;   // doesn't exist
                                  fake0[np0]=1;   // SC:color this particle doesn't exist
                                  P[4][np0]=0.0;  // radiated gluons are massless
                                  P0[4][np0]=0.0;
        
                                  Tint_lrf[i]=0.0;  //reset radiation infomation for heavy quark
        
                                  eGluon=eGluon+pc4[0];
                                  nGluon=nGluon+1.0;
        
                       	          V[0][np0]=-log(1.0-ran0(&NUM1));
        
                                  // add multiple radiation for heavy quark
                                  while(nrad>1) {
        									
                                      for(int j=0;j<=3;j++) pc2[j]=pc4[j];
        
        //			      radiation(qhat0,vc0,pc01,pc2,pc4,pb,iclrad,tcar,tiscatter[i],tirad[i],Elab,Ejp);		   
        			      radiationHQ(KATT1[i],qhat0,vc0,pc2,pc01,pc4,pb,iclrad,Tdiff,Ejp,maxFncHQ,temp0,lim_low,lim_int);  
                                      n_radiation+=1;
        
        			      if(iclrad!=1) {
        			          np0=np0+1;
        			       	  ng_radiation+=1;
        			          KATT1[np0]=21;

                                          ctGluon++;

                                          // SC:color assign color information to additional splittings
                                          if(colorFlag==1){
                                              int parentColor=color[i];  // taken from jet parton from the previous splitting
                                              int parentAntiColor=anticolor[i];
                                              if(KATT1[i]==21) { // g->gg process
                                                  maxColor++;
                                                  color[i]=parentColor;
                                                  anticolor[i]=maxColor;
                                                  color[np0]=maxColor;    // radiated gluon
                                                  anticolor[np0]=parentAntiColor;
                                              } else if(KATT1[i]>0){ // q->qg process
                                                  maxColor++;
                                                  color[i]=maxColor;
                                                  anticolor[i]=0;
                                                  color[np0]=parentColor; // radiated gluon
                                                  anticolor[np0]=maxColor;
                                              } else { // qbar->qbarg
                                                  maxColor++;
                                                  color[i]=0;
                                                  anticolor[i]=maxColor;
                                                  color[np0]=maxColor;
                                                  anticolor[np0]=parentAntiColor;
                                              }
                                          } 
                                          for(int j=0;j<=3;j++) {
                    			      P[j][i]=pc01[j];            // heavy quark after one more gluon
                			      P[j][np0-1]=pc2[j];         // update the previous gluon
                			      P[j][np0]=pc4[j];           // record the new gluon
                			      V[j][np0]=V[j][i];
                			      P0[j][np0]=0.0;
                			      V0[j][np0]=0.0;
        
                                              Vfrozen[j][np0]=V[j][i];
                                              Vfrozen0[j][np0]=0.0;
                                              if(j!=0) {                                 
                                                  vcfrozen[j][np0]=vc0b[j];
                                                  vcfrozen0[j][np0]=0.0;
                                              }
                                          }
        
                                          Tfrozen[np0]=temp0;
                                          Tfrozen0[np0]=0.0;
         
                                          // pass initial pT and weight information
                                          P[5][np0]=P[5][i];
                                          P0[5][np0]=0.0;
                                          WT[np0]=WT[i];
                                          WT0[np0]=0.0;   // doesn't exist
                                          fake0[np0]=1;   // SC:color this particle doesn't exist
                                          P[4][np0]=0.0;  // radiated gluons are massless
                                          P0[4][np0]=0.0;
        
                                          eGluon=eGluon+pc4[0];
                                          nGluon=nGluon+1.0;
        
                              	          V[0][np0]=-log(1.0-ran0(&NUM1));
        
                                      } else { //end multiple radiation
                                          break;
                                      }  //if(icl!=1)radiation
        		   
                                      nrad=nrad-1;
        			  } //multiple radiation for heavy quark ends
        
        //                          cout<<"radiate! <Ng>: "<<nrad0<<"  real Ng H: "<< ctGluon << endl;
                              } // icl23 == 1, 1st gluon radiation from heavy quark
        
                          } // if(ran0(&NUM1)<probRad/probTot) 
        
        	      }  //if(Ejp>2*sqrt(qhat0))				  
        				  				  
        	  }  //if(KINT!=0)
        
        
        	  //...........collision end				  
        	  //...........determine the leading partons and set radng[i]
                      
        
                  //....tiscatter information
                  if(abs(KATT1[i])==4 || abs(KATT1[i])==5) radng[i]=0.0; // do it below
                  tiscatter[i]=tcar;
                  V[0][i]=-log(1.0-ran0(&NUM1));
        
        //          for(unsigned ip=np+1; ip<=np0; ++ip)
                  for(unsigned ip=nnpp+1; ip<=np0; ++ip)
                  {
                        tiscatter[ip]=tcar;
                    	V[0][ip]=-log(1.0-ran0(&NUM1));
                    	V0[0][ip]=-log(1.0-ran0(&NUM1));
                        tirad[ip]=tcar;
                        Tint_lrf[ip]=0.0;
                        radng[ip]=0.0;
                  }
        
        // SC: only do possible switch for g->ggg... here
                  if(KATT1[i]==21 && np0-nnpp>=2) {
                      idlead=nnpp+2;
                      p0temp1=P[0][idlead];
                      for(unsigned ip=nnpp+2; ip<=np0-1; ++ip) {
                          if(p0temp1 < P[0][ip+1]) {
                              p0temp1=P[0][ip+1];
                              idlead=ip+1;
                          }
                      }
        	      if(P[0][idlead]>P[0][i]) {
        	          KATTx=KATT1[i];
        	          KATTy=KATT1[idlead];
        	          KATT1[i]=KATTy;
        	          KATT1[idlead]=KATTx;					
                          for(int j=0;j<=3;j++) {
        	              pcx[j]=P[j][i];
        	              pcy[j]=P[j][idlead];
        	              P[j][i]=pcy[j];
        	              P[j][idlead]=pcx[j];
        	          }
                      }
                  }

	  //...........sorting block i np~np0(positive) np+1(negative)

	  //........................................................................................................

//CAT!!!			  
/*			  
*/			  
                  for(int m=nnpp+1;m<=np0;m++) {
                      if(CAT[i]==2) {
                          CAT[m]=2;             				 
                      }
                  }
/*                // SC: comment out, not necessary, keep all partons, they will do free streaming at the beginning of this function if their energies are low
                  if(P[0][nnpp+1]<Ecut) {
//                      cout << "cut recoil   " << nnpp << "    " << P[0][nnpp+1] << endl; 
                      CAT[nnpp+1]=1;
                      CAT0[nnpp+1]=1;
                      for(int j=0; j<=3; j++) {
                          P[j][nnpp+1]=0.0;
                          P0[j][nnpp+1]=0.0;
                      }
                  }

                  if(P[0][i]<Ecut) {
//                      cout << "cut leading" << endl;
                      CAT[i]=1;
                      for(int j=0; j<=3; j++) P[j][i]=0.0;
                  }

                  for(int m=nnpp+2; m<=np0; m++) {
                      if(P[0][m]<Ecut) {
//                          cout << "cut radiation   " << m << "    " << P[0][m] << endl;
                          CAT[m]=1;
                          for(int j=0; j<=3; j++) P[j][m]=0.0;
                      }
                  }
*/
			  
//CAT!!!			  
/*			  
			  if((i>nj) && (P[0][i]<Ecut))
				{
				  CAT[i]=1;
				  ncut=ncut+1;
//				  for(int j=0;j<=3;j++)
//					{
//					  PP[j][ncut]=P[j][i];
//					  VV[j][ncut]=V[j][i];
//					  P[j][i]=0.0;
//					}
				}			 			 
*/
								
				
//		     cout << i << " " << np << endl;
//             for(unsigned ip=1; ip<=np0; ++ip)
//             {
//                 cout << ip << " " << P[0][ip] << endl;
//	     }
//		     cout << endl;
 
			  //........................................................................................................sorting block for p2 and radiated gluons
              
//CAT!!!			  
/*				  
			  for(int m=np+1;m<=np0;m++)
				{
				  if(P[0][m]<=Ecut)
					{			
					  CAT[m]=1;             
					  ncut=ncut+1;
//					  for(int j=0;j<=3;j++)
//						{
//						  PP[j][ncut]=P[j][m];
//						  VV[j][ncut]=V[j][m];
//						  P[j][m]=0.0;			 
//						}				 
					}
				}
*/

			  //........................................................................................................dump all the negative particles
			  //..................................................................CAT0=0/1  stay in LBT/dumped into hydro (negative partons)
//			  CAT0[np0]=1;             
//			  ncut0=ncut0+1;
//			  for(int j=0;j<=3;j++)
//                {
//				  PP0[j][ncut0]=P0[j][np+1];
//				  VV0[j][ncut0]=V0[j][np+1];
//				  P0[j][np0]=0.0;				 
//                }				
			  //...........................................................................................................++++++++++++++++++++++++END			 				
              } //...!Yes, collision!			

// even if there is no scattering, one should reset last scatter time now in the new framework
              tiscatter[i]=tcar;
              radng[i]=0.0;

          } //....for if(free==0)
      } //..........if CAT[i]=0 end	

//      } //..........if tcar<tiform[i]

      nnpp=np0;  

  } //......for np end

  //........time step end, np: the number of particles at this point  		
  np=np0;
  
  if(Kprimary==1) {
     np=nj;
  }

  if(np>=499999) {
      cout << "np exceeds the grid size ... terminate " << endl;
      exit (EXIT_FAILURE); 
  }

}
		
		
////////////////////////////////////////////////////////////////////////////////////////////////////


