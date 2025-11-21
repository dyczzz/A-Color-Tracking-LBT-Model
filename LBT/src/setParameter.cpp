////////////////////////////////////////////////////////////////////////////////////////

  void setParameter(string fileName) {

       const char* cstring = fileName.c_str();

       ifstream input(cstring);

       if(!input) {
           cout << "Parameter file is missing!" << endl;
           exit(EXIT_FAILURE);
       }

       string line;

       while(getline(input,line)) {
          char str[1024];
          strncpy(str, line.c_str(), sizeof(str));
//          cout << str << endl;

          char* divideChar[5];
          divideChar[0] = strtok(str," ");

          int nChar=0;
          while(divideChar[nChar]!=NULL) {
              if(nChar==3) break;
//              cout << divideChar[nChar] << endl;
              nChar++;
              divideChar[nChar] = strtok(NULL," ");
          }

          if(nChar==3) {
//              cout << divideChar[0] << endl;
//              cout << divideChar[1] << endl;
//              cout << divideChar[2] << endl;
              string firstStr(divideChar[0]);
              string secondStr(divideChar[1]);
              string thirdStr(divideChar[2]);
//              cout << firstStr << "     " << secondStr << endl;
              if(secondStr=="=") {
                  if(firstStr=="flag:fixMomentum") fixMomentum=atoi(divideChar[2]);
                  else if(firstStr=="flag:colorFlag") colorFlag=atoi(divideChar[2]);   // SC:color
                  else if(firstStr=="flag:fixPosition") fixPosition=atoi(divideChar[2]);
                  else if(firstStr=="flag:initHardFlag") initHardFlag=atoi(divideChar[2]);
                  else if(firstStr=="flag:flagJetX") flagJetX=atoi(divideChar[2]);
                  else if(firstStr=="flag:vacORmed") vacORmed=atoi(divideChar[2]);
                  else if(firstStr=="flag:bulkType") bulkFlag=atoi(divideChar[2]);
                  else if(firstStr=="flag:Kprimary") Kprimary=atoi(divideChar[2]);
                  else if(firstStr=="flag:KINT0") KINT0=atoi(divideChar[2]);
                  else if(firstStr=="flag:outFormat") outFormat=atoi(divideChar[2]);
                  else if(firstStr=="flag:heavyOut") heavyOut=atoi(divideChar[2]);
                  else if(firstStr=="flag:lightOut") lightOut=atoi(divideChar[2]);
                  else if(firstStr=="para:nEvent") ncall=atoi(divideChar[2]);
                  else if(firstStr=="para:nParton") nj=atoi(divideChar[2]);
                  else if(firstStr=="para:parID") Kjet=atoi(divideChar[2]);
                  else if(firstStr=="para:tau0") tau0=atof(divideChar[2]);
                  else if(firstStr=="para:tauend") tauend=atof(divideChar[2]);
                  else if(firstStr=="para:dtau") dtau=atof(divideChar[2]);
                  else if(firstStr=="para:temperature") temp0=atof(divideChar[2]);
                  else if(firstStr=="para:alpha_s") fixAlphas=atof(divideChar[2]);
                  else if(firstStr=="para:hydro_Tc") hydro_Tc=atof(divideChar[2]);
                  else if(firstStr=="para:pT_min") ipTmin=atof(divideChar[2]);
                  else if(firstStr=="para:pT_max") ipTmax=atof(divideChar[2]);
                  else if(firstStr=="para:eta_cut") eta_cut=atof(divideChar[2]);
                  else if(firstStr=="para:cutRapidity") cutRapidity=atof(divideChar[2]);
                  else if(firstStr=="para:cutRapidity") cutRapidity=atof(divideChar[2]);
                  else if(firstStr=="para:ener") ener=atof(divideChar[2]);
                  else if(firstStr=="para:mass") amss=atof(divideChar[2]);
                  else if(firstStr=="para:KPamp") KPamp=atof(divideChar[2]);
                  else if(firstStr=="para:KPsig") KPsig=atof(divideChar[2]);
                  else if(firstStr=="para:KTamp") KTamp=atof(divideChar[2]);
                  else if(firstStr=="para:KTsig") KTsig=atof(divideChar[2]);

             }               


          }

       }

// calculated derived quantities
       if(vacORmed == 0) fixPosition = 1; // position is not important in vacuumm

       KTsig=KTsig*hydro_Tc;
       preKT=fixAlphas/0.3;

// write out parameters
       cout << endl;
       cout << "############################################################" << endl;
       cout << "Parameters for this run" << endl;
       cout << "flag:initHardFlag: " << initHardFlag << endl;
       cout << "flag:flagJetX: " << flagJetX << endl;
       cout << "flag:fixMomentum: " << fixMomentum << endl;
       cout << "flag:fixPosition: " << fixPosition << endl;
       cout << "flag:vacORmed: " << vacORmed << endl;
       cout << "flag:bulkType: " << bulkFlag << endl;
       cout << "flag:Kprimary: " << Kprimary << endl;
       cout << "flag:KINT0: " << KINT0 << endl;
       cout << "flag:colorFlag: " << colorFlag << endl;
       cout << "flag:outFormat: " << outFormat << endl;
       cout << "flag:heavyOut: " << heavyOut << endl;
       cout << "flag:lightOut: " << lightOut << endl;
       cout << "para:nEvent: " << ncall << endl;
       cout << "para:nParton: " << nj << endl;
       cout << "para:parID: " << Kjet << endl;
       cout << "para:tau0: " << tau0 << endl;
       cout << "para:tauend: " << tauend << endl;
       cout << "para:dtau: " << dtau << endl;
       cout << "para:cutRapidity: " << cutRapidity << endl;
       cout << "para:temperature: " << temp0 << endl;
       cout << "para:alpha_s: " << fixAlphas << endl;
       cout << "para:hydro_Tc: " << hydro_Tc << endl;
       cout << "para:pT_min: " << ipTmin << endl;
       cout << "para:pT_max: " << ipTmax << endl;
       cout << "para:eta_cut: " << eta_cut << endl;
       cout << "para:ener: " << ener << endl;
       cout << "para:mass: " << amss << endl;
       cout << "para:KPamp: " << KPamp << endl;
       cout << "para:KPsig: " << KPsig << endl;
       cout << "para:KTamp: " << KTamp << endl;
       cout << "para:KTsig: " << KTsig << endl;
       cout << "############################################################" << endl;
       cout << endl;

//       exit(EXIT_SUCCESS);

   
  }

////////////////////////////////////////////////////////////////////////////////////////

  int checkParameter(int nArg) {

      int ctErr=0;

      // better to make sure lightOut,heavyOut,fixPosition,fixMomentum,etc are either 0 or 1
      // not necessary at this moment

      if(colorFlag==1) { // SC:color colorFlag = 1 only compatible with initHardFlag = 2 and outFormat = 1 at this moment 
          if(initHardFlag!=2 || outFormat!=1) {
              cout<<"Color option only supports initHardFlag=2 and outFormat=1 at this moment." << endl;
              exit(EXIT_FAILURE);
          }
      }

      if(initHardFlag==1) { // initialize within LBT, no input particle list
          if(lightOut==1 && heavyOut==1) { // need both heavy and light output
              if(nArg!=5) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file HQ_output light_positive_output light_negative_output"<<endl;
              }
          } else if(heavyOut==1) { // only heavy output
              if(nArg!=3) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file  HQ_output"<<endl;
              }
          } else if(lightOut==1) { // only light output
              if(nArg!=4) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file light_positive_output light_negative_output"<<endl;
              }
          } else { // no output
               if(nArg!=2) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file"<<endl;
                  cout<<"No output specified, both heavyOut and lightOut are set to 0"<<endl;
               }
          }
      } else if(initHardFlag==2) { // need input particle list
          if(lightOut==1 && heavyOut==1) { // need both heavy and light output
              if(nArg!=6) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file input_parton_list HQ_output light_positive_output light_negative_output"<<endl;
              }
          } else if(heavyOut==1) { // only heavy output
              if(nArg!=4) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file input_parton_list HQ_output"<<endl;
              }
          } else if(lightOut==1) { // only light output
              if(nArg!=5) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file input_parton_list light_positive_output light_negative_output"<<endl;
              }
          } else { // no output
               if(nArg!=3) {
                  ctErr++;
                  cout<<"Wrong # of arguments for command line input"<<endl;
                  cout<<"./LBT parameter_file input_parton_list"<<endl;
                  cout<<"No output specified, both heavyOut and lightOut are set to 0"<<endl;
               }
          }
      } else {
          cout<<"Wrong input for initHardFlag!"<<endl;
          ctErr++;
      }

      return(ctErr);

  }



