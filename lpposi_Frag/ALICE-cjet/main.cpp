//#include"LBT.h"

#include "Pythia8/Pythia.h"    //....PYTHIA8 headers
#include "fastjet/ClusterSequence.hh"    //....FASTJET headers
//#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support

#include <math.h>


using namespace Pythia8;
using namespace fastjet;

#include "fastjet/contrib/Recluster.hh" // In external code, this should be fastjet/contrib/Recluster.hh
#include "fastjet/contrib/IteratedSoftDrop.hh" // In external code, this should be fastjet/contrib/IteratedSoftDrop.hh
#include "fastjet/contrib/RecursiveSoftDrop.hh" // In external code, this should be fastjet/contrib/RecursiveSoftDrop.hh
#include "fastjet/contrib/RecursiveSymmetryCutBase.hh"
int searchmax(double*p, int *q,int len);
//....Classes definition
//....FASTJET class dealing with negative particles
typedef fastjet::JetDefinition::Recombiner Recombiner;
/// Recombiner class that propagates the user index and arranges the
/// recombination accordingly
class NegativeEnergyRecombiner : public  Recombiner {
    public:
        NegativeEnergyRecombiner(const int ui) : _ui(ui) {}

        virtual std::string description() const {return "E-scheme Recombiner that checks a flag for a 'negative momentum' particle, and subtracts the 4-momentum in recombinations.";}

        /// recombine pa and pb and put result into pab
        virtual void recombine(const fastjet::PseudoJet & pa, 
                const fastjet::PseudoJet & pb, 
                fastjet::PseudoJet & pab) const {

            int ai=1,bi=1;

            // If a particle is flagged, restore its real negative energy. 
            // The -1 will flip the full 4-momentum, reversing the convention for 
            // negative energy particles.
            if (pa.user_index() < 0) { ai = -1;}
            if (pb.user_index() < 0) { bi = -1;}

            // recombine particles
            pab = ai*pa+bi*pb;

            // if the combination has negative energy, flip the whole 4-momentum and flag it, 
            // so that we have the same convention as for input particles
            if(pab.E() < 0) { 
                pab.set_user_index(_ui); 
                pab.reset_momentum(-pab.px(),
                        -pab.py(),
                        -pab.pz(),
                        -pab.E());
            } else { pab.set_user_index(0);}

        }

    private:
        const int _ui;  
};

class WTAPtNegativeEnergyRecombiner : public  Recombiner {
    public:
        WTAPtNegativeEnergyRecombiner(const int ui) : _ui(ui) {}

        virtual std::string description() const {return "WTA_pt-scheme Recombiner that checks a flag for a 'negative momentum' particle, and subtracts the 4-momentum in recombinations.";}

        /// recombine pa and pb and put result into pab
        virtual void recombine(const fastjet::PseudoJet & pa,
                const fastjet::PseudoJet & pb,
                fastjet::PseudoJet & pab) const {

            int ai=1,bi=1;

            // If a particle is flagged, restore its real negative energy. 
            // The -1 will flip the full 4-momentum, reversing the convention for 
            // negative energy particles.
            if (pa.user_index() < 0) { ai = -1;}
            if (pb.user_index() < 0) { bi = -1;}

            // recombine particles
            const fastjet::PseudoJet &phard = (pa.pt2() >= pb.pt2()) ? pa : pb;
            /// keep y,phi and m from the hardest, sum pt
            pab.reset_PtYPhiM(ai * pa.pt() + bi * pb.pt(), phard.rap(), phard.phi(), phard.m());

            // if the combination has negative energy, flip the whole 4-momentum and flag it, 
            // so that we have the same convention as for input particles
            if(pab.E() < 0) {
                pab.set_user_index(_ui);
                pab.reset_momentum(-pab.px(),
                        -pab.py(),
                        -pab.pz(),
                        -pab.E());
            } else { pab.set_user_index(0);}

        }

    private:
        const int _ui;
};

//////////////////////////////////////////
float ran33(long *idum);
void rotate(double px,double py,double pz,double pr[4],int icc);
double getDistance(double eta1, double eta2, double phi1, double phi2);
//////////////////////////////////////////

int main(int argc, char* argv[]){	  				
    //.... The program starts
    struct tm *local_start;
    time_t time_start;
    time_start = time(NULL);
    local_start = localtime(&time_start);

    char buf1[80];
    strftime(buf1, 80, "Current Time: %Y-%m-%d %H:%M:%S", local_start);
    cout << "the program starts at:" << endl;
    cout << buf1 << endl;
    //.... Time counts

    char charName[1024];

    sprintf(charName, "../hadrons-posi.dat");
    ifstream f_posi(charName);

    sprintf(charName, "../hadrons-nega.dat");
    ifstream f_nega(charName);

    sprintf(charName, "jet_pT.dat");
    ofstream f_pp(charName);

    sprintf(charName, "nega_pT.dat");
    ofstream f_test(charName);


    int numEvent = 0;
    int nEvent = atoi(argv[1]); // changed according to initial PYTHIA events
    int hadron1 = atoi(argv[2]);
    int hadron2 = atoi(argv[3]);
    int hadron3 = atoi(argv[4]);
    double sigma = 0.0;

    double R = 0.4;

    const double phoetacut = 0.0, phoptmin = 0.0, phoptmax = 0.0;
    const double jetetacut = 2.8, jetptmin = 40.0, jetptmax = 1000.0, leadingptcut = 0.0;

    double Njet = 0.0;
//    double alicebin[26] = {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0, 400.0, 500.0};
//    double numsum[25] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double atlasbin[20] = {40.0, 50.0, 60.0, 80.0, 100.0, 112.0, 125.0, 141.0, 158.0, 177.0, 199.0, 223.0, 251.0, 281.0, 316.0, 354.0, 398.0, 501.0, 630.0, 999.0};
    double numsum[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//    double atlasbin[20] = {0.0025, 0.0075, 0.0125, 0.0175, 0.0225, 0.0275, 0.0325, 0.0375, 0.0425, 0.0475, 0.0525, 0.0575, 0.0625, 0.0675, 0.0725, 0.0775, 0.0825, 0.0875, 0.0925, 0.//0975};
//    double numsum[19] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int enumber1=0,enumber2=0,enumber3=0,enumber4=0;
    double neE[19]={0.0}, nepx[19]={0.0}, nepy[19]={0.0}, nepz[19]={0.0};
    double poE[19]={0.0}, popx[19]={0.0}, popy[19]={0.0}, popz[19]={0.0};
    //....Begin event loop. Generate event. Skip if error.
    for (int n = 0; ; ++n) {
        if (numEvent >= nEvent) break;

        bool phoFind = false;

        double pm_flag = 0.0;

        int dummyInt;
        int numParton;
        int pid;
        double peta, ppt, pp; 
        double ppx, ppy, ppz, pe, pm;
        double dummyaa;
        int dummycc;
        int stat;
        int jetType = 0; //quark jet: 0    gluon jet: 1    quark + gluon jet : 2
        int partonid1 = 0, partonid2 = 0;
        double phi1 = 0, phi2 = 0;

        int njfast=-1;
        double E[10000]={0.0}, px[10000]={0.0}, py[10000]={0.0}, pz[10000]={0.0};
        int id[10000]={0};
        fastjet::PseudoJet pj0001;
        vector<fastjet::PseudoJet> input_particles0001;

//......fastjet vector initialize
        int ui = -123456;
        int wta_pt_ui = -123456;
//..................................fastjet loading ends

// **** start loading positive hadron
        f_posi >> dummyInt >> numParton; // >> jetType >> partonid1 >> phi1 >> partonid2 >> phi2;
//        cout << dummyInt << "  " << numParton << endl;
        enumber1 = enumber2;
        enumber2 = enumber3;
        enumber3 = enumber4;
        enumber4 = numParton;
        if (enumber1 == enumber2 && enumber2 == enumber3 && enumber3 == enumber4) break;
        for(int i=1; i<=numParton; i++) //negative
        {
            f_posi >> pid >> ppx >> ppy >> ppz >> pe >> pm >> dummyaa >> dummyaa >> dummyaa >> dummyaa;
            ppt=sqrt(ppx*ppx+ppy*ppy);
            pp=sqrt(ppx*ppx+ppy*ppy+ppz*ppz);
            peta=fabs(0.5*log((pp+ppz)/(pp-ppz)));


            if (fabs(pid) == hadron1 || fabs(pid) == hadron2 || fabs(pid) == hadron3 || fabs(pid) == 111 || fabs(pid) == 311 || fabs(pid) == 2112)
            {
                if (ppt > 0.15 && peta < 3.2)
                {
                    njfast=njfast+1;
                    E[njfast] = pe;
                    px[njfast] = ppx;
                    py[njfast] = ppy;
                    pz[njfast] = ppz;
                    id[njfast] = pid;
                    if ( E[njfast] < 0 ) 
                    { 
                         pj0001.reset_momentum(-px[njfast],-py[njfast],-pz[njfast],-E[njfast]);
                         pj0001.set_user_index(-njfast); 
                    }
                    else 
                    {
                         pj0001.reset_momentum(px[njfast],py[njfast],pz[njfast],E[njfast]);
                         pj0001.set_user_index(njfast);
                    }
                    input_particles0001.push_back(pj0001);			
                }
            }

		
        }//for(int i=1; i<=numParton; i++) //positive
						
// **** start loading positive hadron
        f_nega >> dummyInt >> numParton;
//        cout << dummyInt << "  " << numParton << endl;
        for(int i=1; i<=numParton; i++) //negative
        {
            f_nega >> pid >> ppx >> ppy >> ppz >> pe >> pm >> dummyaa >> dummyaa >> dummyaa >> dummyaa;
            ppx = -ppx;
            ppy = -ppy;
            ppz = -ppz;
            pe = -pe;
            ppt=sqrt(ppx*ppx+ppy*ppy);
            pp=sqrt(ppx*ppx+ppy*ppy+ppz*ppz);
            peta=fabs(0.5*log((pp+ppz)/(pp-ppz)));

            if (fabs(pid) == hadron1 || fabs(pid) == hadron2 || fabs(pid) == hadron3 || fabs(pid) == 111 || fabs(pid) == 311 || fabs(pid) == 2112)
            {
                if (ppt > 0.15 && peta < 3.2)
                {
                    njfast=njfast+1;
                    E[njfast] = pe;
                    px[njfast] = ppx;
                    py[njfast] = ppy;
                    pz[njfast] = ppz;
                    id[njfast] = pid;
                    if ( E[njfast] < 0 ) 
                    { 
                         pj0001.reset_momentum(-px[njfast],-py[njfast],-pz[njfast],-E[njfast]);
                         pj0001.set_user_index(-njfast); 
                    }
                    else 
                    {
                         pj0001.reset_momentum(px[njfast],py[njfast],pz[njfast],E[njfast]);
                         pj0001.set_user_index(njfast);
                    }
                    input_particles0001.push_back(pj0001);			
                }
            }

		
        }//for(int i=1; i<=numParton; i++) //positive
						
						
//..................................fastjet loading ends


//***** create a jet definition: 
//***** a jet algorithm with a given radius parameter	  
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
//***** create an instance of the negative energy recombiner, with a given flag ui
        NegativeEnergyRecombiner uir(ui);
//***** tell jet_def to use this new recombiner
        jet_def.set_recombiner(&uir);
//***** run the jet clustering with the above jet definition
//----------------------------------------------------------
        fastjet::ClusterSequence clust_seq0001(input_particles0001, jet_def);
//***** get the resulting jets ordered in pt
//----------------------------------------------------------
        vector<fastjet::PseudoJet> inclusive_jets0001 = sorted_by_pt(clust_seq0001.inclusive_jets(jetptmin));


//cc	cout << "D0 number is: " << "  " << nD0 << endl;
//cc	cout << "the constructed jet size is: " << "  " << inclusive_jets0001.size() << endl;


        for (int ij = 0; ij < inclusive_jets0001.size(); ij++) {
            int jetwanted = 0;
            double ldpx = 0, ldpy = 0, ldpz = 0, lde = 0, ldpt=0, ldeta=0;
            double ldphi = 0, R_axis = 0, value_pi = 3.1415926, ptsum = 0;
            double posiE = 0, posipx =0, posipy = 0, posipz = 0;
            double negaE = 0, negapx =0, negapy = 0, negapz = 0;
            ldpx = inclusive_jets0001[ij].px();
            ldpy = inclusive_jets0001[ij].py();
            ldpz = inclusive_jets0001[ij].pz();
            lde = inclusive_jets0001[ij].e();
            ldpt = inclusive_jets0001[ij].pt();
            ldeta = inclusive_jets0001[ij].eta();
//cc            ldphi = inclusive_jets0001[ij].phi_std();
            ldphi = atan2(ldpy,ldpx);
            if (ldpt > jetptmax) continue;
            if (fabs(ldeta) > jetetacut) continue;
            double idflag = 0.0;
            double LHptflag = 0.0;
            double ptmax = 0.0;
            if (inclusive_jets0001[ij].user_index()>=0)
            {
//                 jetwanted = jetwanted + 1;
                 vector<PseudoJet> constituents = inclusive_jets0001[ij].constituents();
                 for (int iD = 0; iD < constituents.size(); iD++)
                 { 
                      if(constituents[iD].user_index()>0)
                      {    
                      posiE += constituents[iD].e();
                      posipx += constituents[iD].px();
                      posipy += constituents[iD].py();
                      posipz += constituents[iD].pz();
                      int index = fabs(constituents[iD].user_index());
                      if(fabs(id[index]) == hadron1 || fabs(id[index]) == hadron2 || fabs(id[index]) == hadron3) {
                        if (constituents[iD].pt() > 4) ptsum = ptsum + constituents[iD].pt();
                      }
                      if (ptmax < constituents[iD].pt()) ptmax = constituents[iD].pt();
                 
                      if(constituents[iD].pt()>=LHptflag) LHptflag = constituents[iD].pt();

                      } else {
                        negaE += constituents[iD].e();
                        negapx += constituents[iD].px();
                        negapy += constituents[iD].py();
                        negapz += constituents[iD].pz();
                      }
                 }
            }
            if ((ptmax > leadingptcut)) ++jetwanted;

            if (ptsum < 8.0) continue;
//            if (jetwanted==0) continue;
//            if (LHptflag<3.0) continue;
/*
            fastjet::JetDefinition wta_jet_def(fastjet::cambridge_algorithm, R, fastjet::WTA_pt_scheme);
            WTAPtNegativeEnergyRecombiner wta_pt_uir(wta_pt_ui);
            wta_jet_def.set_recombiner(&wta_pt_uir);
            fastjet::Recluster wta_recluster(wta_jet_def);
            PseudoJet wta_jet = wta_recluster(inclusive_jets0001[ij]);

            double wtaeta = wta_jet.eta();
            double wtaphi = wta_jet.phi_std();

            phoFind = true;

            R_axis = getDistance(ldeta, wtaeta, ldphi, wtaphi);
        cout << R_axis << " " << wtaphi << " " << ldphi  << " " << wtaeta  << " " << ldeta << endl;
*/
//----------------------------------------------------------------
            int numflag = 0.0;
            if (inclusive_jets0001[ij].user_index()>=0)
            {
            for (int k=0; k<19; k++)
            {
                if ((ldpt>=atlasbin[k])&&(ldpt<atlasbin[k+1]))
                {
                     numflag = numflag + 1;
                     numsum[k] = numsum[k] + 1.0;
                     neE[k] += negaE;
                     nepx[k] += negapx;
                     nepy[k] += negapy;
                     nepz[k] += negapz;
                     poE[k] += posiE;
                     popx[k] += posipx;
                     popy[k] += posipy;
                     popz[k] += posipz;
                }
            }

            if (numflag >= 1) 
            {
                Njet = Njet + 1.0;
            }
            }



        }
//        if (!phoFind) continue;

        //  parton information
        numEvent++;

    } // event loop ends

    cout << numEvent  << "  " << Njet << endl;
    f_pp << "0" << "  " << numEvent << "  " << Njet << "  " <<  "0.0" << endl;
    f_test << "# pT    posiE    negaE    posipx    negapx    posipy    negapy    posipz    negapz" << endl;

    for (int k = 0; k < 19; k++)
    {
         double var = (atlasbin[k]+atlasbin[k+1])/2.0;
         double xerro = (atlasbin[k+1]-atlasbin[k])/2.0;
         double yvar = numsum[k] / (atlasbin[k+1]-atlasbin[k]);
         cout << "k: " << " " << k << " " << numsum[k] << endl;
         f_pp << var << "  " << yvar << "  " << xerro << "  " << xerro << endl;
         f_test << var << "    " << poE[k] << "    " << neE[k] << "    " << popx[k] << "    " << nepx[k] << "    " << popy[k] << "    " << nepy[k] << "    " << popz[k] << "    " << nepz[k] << endl;
    }


    f_posi.close();
    f_nega.close();
    f_pp.close();
    f_test.close();

    struct tm *local_end;
    time_t time_end;
    time_end = time(NULL);
    local_end = localtime(&time_end);

    char buf2[80];
    strftime(buf2, 80, "Current Time: %Y-%m-%d %H:%M:%S", local_end);
    cout << "the program ends at:" << endl;
    cout << buf2 << endl;

    int cost, nh, nm, ns;
    cost = difftime(time_end, time_start);

    nh = cost / 3600;
    nm = (cost % 3600) / 60;
    ns = (cost % 3600) % 60;

    cout << "the program costs:" << endl;
    cout << cost << "s:" << " " << nh << "h" << " " << nm << "m" << " " << ns << "s" << endl;
}


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

float ran33(long *idum)

{
    int j;
    long k;
    static long idum2=123456789;
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

    //  if(pt==0)
    if(pt<1e-6)
    {
        cosa=1;
        sina=0;
    } 
    else
    {
        cosa=px/pt;
        sina=py/pt;
    }

    if(E>1e-6) {

        cosb=pz/E;
        sinb=pt/E;

        if(icc==1) {
            wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
            wy1=-wx*sina+wy*cosa;
            wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
        } else {
            wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
            wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
            wz1=-wx*sinb+wz*cosb;
        }
        wx=wx1;
        wy=wy1;
        wz=wz1;
    } else {
        cout << "warning: small E in rotation" << endl;
    }

    pr[1]=wx;
    pr[2]=wy;
    pr[3]=wz;      

    //  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);

}


double getDistance(double eta1, double eta2, double phi1, double phi2){

       double value_pi = 3.1415926;
       if (phi1 < 0.0) phi1 = phi1 + 2.0*value_pi;
       if (phi2 < 0.0) phi2 = phi2 + 2.0*value_pi;
       double delta_eta = fabs(eta1 - eta2);
       double delta_phi = fabs(phi1 - phi2);
       if (delta_phi > value_pi) delta_phi = 2.0*value_pi - delta_phi;
       double dist = sqrt(pow(delta_eta,2)+pow(delta_phi,2));
       return dist;
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

