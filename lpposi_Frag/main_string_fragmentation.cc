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

    int jobid=1197175;
    std::string str = std::to_string(jobid);
    string ipput_filename2;

    ipput_filename2 = path2+"lp-posi.dat";
    sprintf(infiles, ipput_filename2.c_str()); //input parton file 
    FILE* infile1;
    infile1 = fopen(infiles,"r");

    char charName[1024];
    string ipput_filename3;
    ipput_filename3 = path2+"lp-nega.dat";
    sprintf(charName, ipput_filename3.c_str());
    ifstream f_nega(charName);

    Pythia pythia;
    pythia.readFile(argc[5]);
    pythia.readFile(argc[6]);
    pythia.init();

    double hbarc = 0.19732;
    double c_px, c_py ,c_pz, c_e, c_m,c_x,c_y,c_z,c_t,c_c,c_ac;
    double c_cat;
    int c_id,tt;
    int acol_ip[5000]={0};
    int pp_collision=0;    //used for total cross section
    int  ccbar_num=0;
    double cmeson_px,cmeson_py,cmeson_pz, cmeson_P,tempp;
    double cbar_meson_px,cbar_meson_py,cbar_meson_pz, cbar_meson_P,pt_square,cbar_meson_energy;
    double pt,amid,Qmid;
    int mid,ie,II,status,col,acol,Npart,NN,Ntotal,PAosi,simble,mmaxindex,eventid;
    int idp[100000]={0},idpo[100000]={0};
    double pxpo[100000]={0.0},pypo[100000]={0.0},pzpo[100000]={0.0},epo[100000]={0.0},ptpo[100000]={0.0},xxpo[100000]={0.0},yypo[100000]={0.0},zzpo[100000]={0.0},ttpo[100000]={0.0},phio[100000]={0.0},etao[100000]={0.0},mass[10000]={0.0};
    double pxp[100000]={0.0},pyp[100000]={0.0},pzp[100000]={0.0},ep[100000]={0.0},ptp[100000]={0.0},xxp[100000]={0.0},yyp[100000]={0.0},zzp[100000]={0.0},ttp[100000]={0.0},phi[100000]={0.0},distance[100000]={0.0},dsting[1000][10000]={0.0},Qscale[100000]={0.0};//,mass[100000]={0.0};
    int nncol[100000]={0},aacol[100000]={0},index[100000]={0},qindex[100000]={0},aqindex[100000]={0},used[100000]={0},pair[10000]={0},apair[10000]={0},gindex[10000]={0},strings[1000][10000]={0},Snum[10000]={0};//,nncol_mid[100000]={0},aacol_mid[100000]={0};
    // event loop
    int jetType = 0; //quark jet: 0    gluon jet: 1    quark + gluon jet : 2
    int partonid1 = 0, partonid2 = 0;
    double phi1 = 0, phi2 = 0, eweight = 0;


    for (int iEvent=0; iEvent<n_event; iEvent++) {  
        int Npart1=0,Npart2=0;  
        fscanf(infile1,"%d %d %d %d\n",&eventid,&Npart1,&partonid1,&partonid2);
//        cout << eventid << " " << Npart1 <<endl;
        if (Npart1==0 ) {output2 << iEvent+1<<" "<<0 << endl;continue;} 
        int Nquark=0;
        int Naquark =0;
        int Ngluon=0;
        int Npair=0;
        for (int ll=0;ll<Npart1;ll++) {
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
                if(ll<partonid1) Qscale[ll]=pT_cut; else Qscale[ll]=0.0;
                double aamid=0.5*log((pmg+c_pz)/(pmg-c_pz));
                etao[ll]=aamid;
                phio[ll]=atan2(c_py,c_px);
                index[ll]=ll;
                used[ll]=0;
        }
//////        fscanf(f_hq,"%d %d\n",&mid,&Npart);
        f_nega >> mid >> Npart2;
//        cout  << mid << " " << Npart2 <<endl;
        Npart=Npart1+Npart2;
        for (int ll=Npart1;ll<(Npart1+Npart2);ll++) {
//////                fscanf(f_hq,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",&mid,&c_id,&c_px,&c_py,&c_pz,&c_e,&c_x, &c_y, &c_z,&c_t);// format of input partons
                f_nega >> mid >> c_id >> c_px >> c_py >> c_pz >> c_e >> c_m >> c_x >> c_y >> c_z >> c_t >> c_c >> c_ac;
                Qmid = 0.0; // The scale for the parton shower.
                if(c_id==21) {
                    idpo[ll]=c_id;
                } else {
                    idpo[ll]=-c_id;
                } 
                if(c_id==21) {gindex[Ngluon]=ll;Ngluon++;}
                pxpo[ll]=0.1;
                pypo[ll]=0.1;
                pzpo[ll]=0.1;
                ptpo[ll]=c_px*c_px+c_py*c_py;
                double pmg=sqrt(c_px*c_px+c_py*c_py+c_pz*c_pz);
                epo[ll]=c_e;
                mass[ll]=c_m;
                xxpo[ll]=c_x;
                yypo[ll]=c_y;
                zzpo[ll]=c_z;
                ttpo[ll]=c_t;
                Qscale[ll]=sqrt(Qmid);
                double aamid=0.5*log((pmg+c_pz)/(pmg-c_pz));
                etao[ll]=aamid;
                phio[ll]=atan2(c_py,c_px);
                index[ll]=ll;
                used[ll]=0;
                nncol[ll]=c_ac;
                aacol[ll]=c_c;
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
            pythia.event.append(idpo[tt],62,nncol[tt],aacol[tt],pxpo[tt],pypo[tt],pzpo[tt],epo[tt],mass[tt],Qscale[tt]);
//            if(maxQ0<Qscale[tt])Qscale[tt]=maxQ0;
            //if(minQ0>Qscale[tt])Qscale[tt]=minQ0;
            //pythia.event[tt].scale(Qscale[tt]);//QQ0 the initial scale of input partons 
            // get the center of mass of the strings and corresponding posistion 
            m_str=m_str+mass[tt];
            x_str=x_str+mass[tt]*xxpo[tt];
            y_str=y_str+mass[tt]*yypo[tt];
            z_str=z_str+mass[tt]*zzpo[tt];
            t_str=t_str+mass[tt]*ttpo[tt];
        }
        if(m_str==0)m_str=0.10;
//****** fragment the remnant partons ***********
//        pythia.forceHadronLevel();
        pythia.forceTimeShower(1,Npart,1000);//Continue the FSR to the defaulted scale 
        pythia.forceHadronLevel();
        int simble = 0;
        for(int i=0; i<pythia.event.size();i++) {
            if (pythia.event[i].isFinal() ) {
                c_id = pythia.event[i].id();
//                if( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                    simble=simble+1;
//                }
            }
        }

        if(simble==0){output2 << iEvent+1<<" "<<simble<< endl;}
        if(simble>0){
            output2 << iEvent+1<<" "<<simble <<  endl;
//            cout << simble << endl;
            for(int i=0; i<pythia.event.size();i++)
                {
                if (pythia.event[i].isFinal() ){
                    c_id = pythia.event[i].id();
                    int color = pythia.event[i].col();
                    int acolor = pythia.event[i].acol();
//                    if ( (abs(c_id) !=22)&&(abs(c_id) !=11)) {
                        cbar_meson_px = pythia.event[i].px();
                        cbar_meson_py = pythia.event[i].py();
                        cbar_meson_pz = pythia.event[i].pz();
                        cbar_meson_energy = pythia.event[i].e();
                        // get the posistion of the final hadrons
                        hmt=sqrt(pythia.event[i].m()*pythia.event[i].m()+cbar_meson_px*cbar_meson_px+cbar_meson_py*cbar_meson_py);
                        x_hadron=x_str/m_str+hbarc*cbar_meson_px/hmt;
                        y_hadron=y_str/m_str+hbarc*cbar_meson_py/hmt;
                        z_hadron=z_str/m_str+hbarc*cbar_meson_pz/hmt;
                        t_hadron=t_str/m_str+hbarc*cbar_meson_energy/hmt;
                        output2 << c_id << " "<<cbar_meson_px<<"  "<<cbar_meson_py<<"  "<<cbar_meson_pz <<" "<<cbar_meson_energy << "  "<< pythia.event[i].m()<<" "<< x_hadron<<" "<<y_hadron<<" "<<z_hadron<<" "<<t_hadron<<endl;
//                    }
                }
            }
        }
        pythia.next();

}
output2.close();
f_nega.close();

  return 0;
}

int searchmin(double*p,int*q,int*s,int len)
{
    double m = 100000000.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if ((m > p[i])&&(s[i]==0))
        {
            m = p[i];
            k = i;
        }
    }
    return q[k];
} 


int searchmax(double*p,int*q,int len)
{
    double m = 0.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if (m < p[i])
        {
            m = p[i];
            k = i;
        }
    }
    return q[k];
} 

    
int searchmin2(double*p,int len)
{
    double m = 100000000.0;
    int k;
    for (int i = 0; i < len; ++i)
    {
        if (m > p[i])
        {
            m = p[i];
            k = i;
        }
    }
    return k;
} 
    
