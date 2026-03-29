#include <assert.h>
#include <time.h>
#include <sstream>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <vector>
#include "Pythia8/Pythia.h" 
#define PI 3.1415926
#define QQ0 1.0

using namespace Pythia8;
int searchmin(double*p, int *q,int*s,int len);
int searchmax(double*p, int *q,int len);
int searchmin2(double*p,int len);

char infiles[128];

//===========================================================================================
int main(int argv, char* argc[])
{
    string random_str = string(argc[1]);
    int n_event=atoi(argc[1]);
    double pT_cut=atoi(argc[2]);
    string path=string(argc[3]);
    string path2=string(argc[4]);
    string output_filename2;
    output_filename2 = path+"hadrons-posi.dat";// output files of final hadrons
    cout << output_filename2 << endl;
    string ramdomseed_str = "Random:seed = "+random_str;
    ofstream output2(output_filename2.c_str()); 
    if (!output2.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename2 << endl;
        return -1;
    }

    int jobid=1085153;
    std::string str = std::to_string(jobid);
    string ipput_filename2;
    ipput_filename2 = path2+"jet_shower_parton.dat";
    sprintf(infiles, ipput_filename2.c_str()); //input parton file 
    FILE* infile1;
    infile1 = fopen(infiles,"r");

    Pythia pythia;
    pythia.readFile(argc[5]);
    pythia.readFile(argc[6]);
    pythia.init();

    double hbarc = 0.19732, eweight=0;
    double c_px, c_py ,c_pz, c_e, c_m,c_x,c_y,c_z,c_t,c_c,c_ac,c_status;
    int c_id,tt;
    std::vector<int> acol_ip(5000, 0);
    int pp_collision=0;    //used for total cross section
    int  ccbar_num=0;
    double cmeson_px,cmeson_py,cmeson_pz, cmeson_P,tempp;
    double cbar_meson_px,cbar_meson_py,cbar_meson_pz, cbar_meson_P,pt_square,cbar_meson_energy;
    double pt,pm,amid,Qmid,weight,weightSum,ptHat;
    int mid,ie,II,col,acol,Npart,NN,Ntotal,PAosi,simble,mmaxindex;
    

    std::vector<int> idp(100000, 0);
    std::vector<int> idpo(100000, 0);
    std::vector<double> pxpo(100000, 0.0);
    std::vector<double> pypo(100000, 0.0);
    std::vector<double> pzpo(100000, 0.0);
    std::vector<double> epo(100000, 0.0);
    std::vector<double> ptpo(100000, 0.0);
    std::vector<double> xxpo(100000, 0.0);
    std::vector<double> yypo(100000, 0.0);
    std::vector<double> zzpo(100000, 0.0);
    std::vector<double> ttpo(100000, 0.0);
    std::vector<double> phio(100000, 0.0);
    std::vector<double> etao(100000, 0.0);
    std::vector<double> mass(100000, 0.0);
    std::vector<double> status(100000, 0.0);
    
    std::vector<double> pxp(100000, 0.0);
    std::vector<double> pyp(100000, 0.0);
    std::vector<double> pzp(100000, 0.0);
    std::vector<double> ep(100000, 0.0);
    std::vector<double> ptp(100000, 0.0);
    std::vector<double> xxp(100000, 0.0);
    std::vector<double> yyp(100000, 0.0);
    std::vector<double> zzp(100000, 0.0);
    std::vector<double> ttp(100000, 0.0);
    std::vector<double> phi(100000, 0.0);
    std::vector<double> distance(100000, 0.0);
    

    std::vector<std::vector<double> > dsting(1000, std::vector<double>(10000, 0.0));
    std::vector<std::vector<int> > strings(1000, std::vector<int>(10000, 0));
    
    std::vector<double> Qscale(100000, 0.0);
    std::vector<int> nncol(100000, 0);
    std::vector<int> aacol(100000, 0);
    std::vector<int> index(100000, 0);
    std::vector<int> qindex(100000, 0);
    std::vector<int> aqindex(100000, 0);
    std::vector<int> used(100000, 0);
    std::vector<int> aused(100000, 0);
    std::vector<int> uncol(10000, 0);
    std::vector<int> unacol(10000, 0);
    std::vector<int> gindex(10000, 0);
    std::vector<int> Snum(10000, 0);

    // event loop
    for (int iEvent=0; iEvent<n_event; iEvent++) {    
        fscanf(infile1,"%d %d\n",&mid,&Npart);
        if (Npart==0 ) {output2 << iEvent+1<<" "<<0 << endl;continue;} 
        int Nquark=0;
        int Naquark =0;
        int Ngluon=0;
        int Npair=0;
        

        if (Npart > idpo.size()) {
            idpo.resize(Npart + 1000);
            pxpo.resize(Npart + 1000);
            pypo.resize(Npart + 1000);
            pzpo.resize(Npart + 1000);
            epo.resize(Npart + 1000);
            ptpo.resize(Npart + 1000);
            xxpo.resize(Npart + 1000);
            yypo.resize(Npart + 1000);
            zzpo.resize(Npart + 1000);
            ttpo.resize(Npart + 1000);
            phio.resize(Npart + 1000);
            etao.resize(Npart + 1000);
            mass.resize(Npart + 1000);
            status.resize(Npart + 1000);
            nncol.resize(Npart + 1000);
            aacol.resize(Npart + 1000);
            Qscale.resize(Npart + 1000);
            index.resize(Npart + 1000);
            used.resize(Npart + 1000);
            aused.resize(Npart + 1000);
        }
        
        for (int ll=0;ll<Npart;ll++) {
            fscanf(infile1,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&mid,&c_id,&c_px,&c_py,&c_pz,&c_e,&c_m,&c_x, &c_y, &c_z,&c_t,&c_c,&c_ac);// format of input partons
            Qmid = 0.0; // The scale for the parton shower. 
            idpo[ll]=c_id;
            if(c_id==21){gindex[Ngluon]=ll;Ngluon++;}
            pxpo[ll]=c_px;
            pypo[ll]=c_py;
            pzpo[ll]=c_pz;
            ptpo[ll]=c_px*c_px+c_py*c_py;
            double pmg=sqrt(c_px*c_px+c_py*c_py+c_pz*c_pz);
            epo[ll]=c_e;
            mass[ll]=c_m;
            xxpo[ll]=c_x;
            yypo[ll]=c_y;
            zzpo[ll]=c_z;
            ttpo[ll]=c_t;
            nncol[ll]=c_c;
            aacol[ll]=c_ac;
            status[ll]=c_status;
            Qscale[ll]=sqrt(Qmid);
            double aamid=0.5*log((pmg+c_pz)/(pmg-c_pz));
            etao[ll]=aamid;
            phio[ll]=atan2(c_py,c_px);
            index[ll]=ll;
            used[ll]=0;
        }

// ****************** append the partons into pythia event *********************//

        double m_str=0.0, x_str=0.0,y_str=0.0,z_str=0.0,t_str=0.0;
        double x_hadron,y_hadron,z_hadron,t_hadron,hmt;
        pythia.event.clear();
        double maxQ0 = 2.0;//maxQ0>=QQ0; wenbin 
        double minQ0 = 0.4;//minum Q0;
        pythia.event.reset();
        for (int tt=0;tt<Npart;tt++){
            if(idpo[tt]==21){mass[tt]=0.0;}
            else{
                if(abs(idpo[tt])<=2)mass[tt]=0.330;
                if(abs(idpo[tt])==3)mass[tt]=0.50;
                if(abs(idpo[tt])==4)mass[tt]=1.50;
                if(abs(idpo[tt])==5)mass[tt]=4.80;
                epo[tt]=sqrt(abs(mass[tt]*mass[tt]+pzpo[tt]*pzpo[tt]+pypo[tt]*pypo[tt]+pxpo[tt]*pxpo[tt]));
            }
            pythia.event.append(idpo[tt],62,nncol[tt],aacol[tt],pxpo[tt],pypo[tt],pzpo[tt],epo[tt],mass[tt],pT_cut);
            if(maxQ0<Qscale[tt])Qscale[tt]=maxQ0;
            m_str=m_str+mass[tt];
            x_str=x_str+mass[tt]*xxpo[tt];
            y_str=y_str+mass[tt]*yypo[tt];
            z_str=z_str+mass[tt]*zzpo[tt];
            t_str=t_str+mass[tt]*ttpo[tt];
        }
        if(m_str==0)m_str=0.10;
        
        pythia.forceTimeShower(1,Npart,1000);//Continue the FSR to the defaulted scale



        pythia.forceHadronLevel();
        int simble = 0, pflag = 0;
        for(int i=0; i<pythia.event.size();i++) {
            if (pythia.event[i].isFinal() ) {
                c_id = pythia.event[i].id();
                if( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                    simble=simble+1;
                }
                if(abs(c_id)==21) pflag=pflag+1;
            }
        }

        if(simble==0){output2 << iEvent+1<<" "<<simble  << endl;}
        if(simble>0){
            output2 << iEvent+1<<" "<<simble << endl;
            for(int i=0; i<pythia.event.size();i++) {
                if (pythia.event[i].isFinal() ){
                    c_id = pythia.event[i].id();
                    if ( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                        cbar_meson_px = pythia.event[i].px();
                        cbar_meson_py = pythia.event[i].py();
                        cbar_meson_pz = pythia.event[i].pz();
                        cbar_meson_energy = pythia.event[i].e();
                        hmt=sqrt(pythia.event[i].m()*pythia.event[i].m()+cbar_meson_px*cbar_meson_px+cbar_meson_py*cbar_meson_py);
                        x_hadron=x_str/m_str+hbarc*cbar_meson_px/hmt;
                        y_hadron=y_str/m_str+hbarc*cbar_meson_py/hmt;
                        z_hadron=z_str/m_str+hbarc*cbar_meson_pz/hmt;
                        t_hadron=t_str/m_str+hbarc*cbar_meson_energy/hmt;
                        output2 << c_id << " "<<cbar_meson_px<<"  "<<cbar_meson_py<<"  "<<cbar_meson_pz <<" "<<cbar_meson_energy << "  "<< pythia.event[i].m()<<" "<< x_hadron<<" "<<y_hadron<<" "<<z_hadron<<" "<<t_hadron<<" "<<endl;
                    }
                }
            }
        }
        pythia.next();
    }
    output2.close();
    return 0;
}

