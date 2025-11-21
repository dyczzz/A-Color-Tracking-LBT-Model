//#include"LBT.h"

#include "Pythia8/Pythia.h"    //....PYTHIA8 headers
#include "fastjet/ClusterSequence.hh"    //....FASTJET headers
//#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support

#include <math.h>

//#include "TFile.h"    //....PYTHIA8 headers
//#include "TH1.h"    //....PYTHIA8 headers

using namespace Pythia8;
using namespace fastjet;

#include "fastjet/contrib/Recluster.hh" // In external code, this should be fastjet/contrib/Recluster.hh
#include "fastjet/contrib/IteratedSoftDrop.hh" // In external code, this should be fastjet/contrib/IteratedSoftDrop.hh
#include "fastjet/contrib/RecursiveSoftDrop.hh" // In external code, this should be fastjet/contrib/RecursiveSoftDrop.hh
#include "fastjet/contrib/RecursiveSymmetryCutBase.hh"

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


//////////////////////////////////////////
float ran33(long *idum);
void rotate(double px,double py,double pz,double pr[4],int icc);
double getDistance(double eta1, double eta2, double phi1, double phi2);
//////////////////////////////////////////

int main(int argc, char **argv){	  				
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

    //.... PYTHIA generator. Initialization. 
    Pythia pythia;
    Event& event = pythia.event;
    int hqid = atoi(argv[1]);
    pythia.readFile(argv[2]);
    int nEvent = pythia.mode("Main:numberOfEvents");
    pythia.readFile(argv[3]);
    pythia.init();


    char charName[1024];
    

    sprintf(charName, "jet_hadron.dat");
    ofstream f_hadron(charName);

    sprintf(charName, "jet_shower_parton.dat");
    ofstream f_pp(charName);

    sprintf(charName, "jet_sigmainfo.dat");
    ofstream f_info(charName);


    int numEvent = 0;
    int numHQ = 0;
    double sigma = 0.0;

    double R = 0.3;

    const double phoetacut = 0.0, phoptmin = 0.0, phoptmax = 0.0;
    const double trketacut = 0.6, trkptmin = 0.15, trkptmax = 300.0;
    const double jetetacut = 0.6, jetptmin = 5.0, jetptmax = 7000.0;

    double Njet = 0.0;


    //....Begin event loop. Generate event. Skip if error.
    for (int n = 1; ; ++n) {
        if (!pythia.next()) continue;
        if (numEvent == nEvent) break;
        bool phoFind = false;

	bool D0Find = false;
	int nD0 = 0;

        double pm_flag = 0.0;

        fastjet::PseudoJet pj0001;
        vector<fastjet::PseudoJet> input_particles0001;
        for (int i = 0; i < event.size(); ++i) {
            int pid = event[i].id();
            double peta = fabs(event[i].eta());
            double ppt = event[i].pT();
            double ppx = event[i].px();
            double ppy = event[i].py();
            double ppz = event[i].pz();
            double pe = event[i].e();

            if (event[i].isFinalPartonLevel()) {
                if (fabs(pid) == hqid) ++pm_flag;
              }


            if (event[i].isFinal()) {

            if (fabs(pid) == 2212 || fabs(pid) == 211 || fabs(pid) == 321 || fabs(pid) == 411 || fabs(pid) == 421 || fabs(pid) == 4122) {
//              if (fabs(pid) == 2212 || fabs(pid) == 321 || fabs(pid) == 211 || fabs(pid) == 111 || fabs(pid) == 311 || fabs(pid) == 2112 || fabs(pid) == 511 || fabs(pid) == 521 || fabs(pid) == 5122) {  
//              if (fabs(pid) == 2212 || fabs(pid) == 321 || fabs(pid) == 211 || fabs(pid) == 111 || fabs(pid) == 311 || fabs(pid) == 2112) {
//                   if (fabs(pid) == 321 || fabs(pid) == 211 || fabs(pid) == 2212) {
                   if (ppt > 0.15 && peta < 0.9)
                    {
                    pj0001.reset_momentum(ppx, ppy, ppz, pe);
                    pj0001.set_user_index(i);
                    input_particles0001.push_back(pj0001);
                    }
               }              

            }


        }


//        if (pm_flag == 0) continue;



//----------------------------------------------------------------
//------start with an example from anti-kt jets
//cc        cout << "--------------------------------------------------" << endl;

        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
        fastjet::ClusterSequence clust_seq0001(input_particles0001, jet_def);
        vector<fastjet::PseudoJet> inclusive_jets0001 = sorted_by_pt(clust_seq0001.inclusive_jets(jetptmin));
        if (inclusive_jets0001.size() == 0) continue;

        for (int ij = 0; ij < inclusive_jets0001.size(); ij++) {
            int jetwanted = 0;
            double ldpx = 0, ldpy = 0, ldpz = 0, lde = 0, ldpt=0, ldeta=0;
            double ldphi = 0;
            ldpx = inclusive_jets0001[ij].px();
            ldpy = inclusive_jets0001[ij].py();
            ldpz = inclusive_jets0001[ij].pz();
            lde = inclusive_jets0001[ij].e();
            ldpt = inclusive_jets0001[ij].pt();
            ldeta = fabs(inclusive_jets0001[ij].eta());
//cc            ldphi = inclusive_jets0001[ij].phi_std();
            ldphi = atan2(ldpy,ldpx);

/////////////            if ((ldpt<=jetptmin) || (ldpt>=jetptmax)) continue;
            if ((ldpt<=jetptmin)) continue;
            if (ldeta > jetetacut) continue; 
//            if ((ldpt>jetptmax)) continue;
   phoFind = true;
 
        }
        if (!phoFind) continue;


//*************************************************************************************************//
// calculate formation time for each shower parton
        int npar = 0;
        int irecord=1;
        int pidwant[5000]={0};
        double pxwant[5000]={0.0}, pywant[5000]={0.0};
        double pzwant[5000]={0.0}, pewant[5000]={0.0};
        double pmwant[5000]={0.0};
        double timewant[5000]={0.0};
        double xwant[5000]={0.0}, ywant[5000]={0.0};
        double zwant[5000]={0.0}, eweight=pythia.info.weight();
        int cwant[5000]={0};
        int acwant[5000]={0};
        int statuswant[5000]={0};

        int nhad = 0;
        int pidhadron[5000]={0};
        double pxhadron[5000]={0.0}, pyhadron[5000]={0.0};
        double pzhadron[5000]={0.0}, pehadron[5000]={0.0};
        double pmhadron[5000]={0.0};
        double thadron[5000]={0.0};
        double xhadron[5000]={0.0}, yhadron[5000]={0.0};
        double zhadron[5000]={0.0};

        for (int i = 0; i < event.size(); ++i) {
            int pid = event[i].id();
            double peta = fabs(event[i].eta());
            double ppt = event[i].pT();
            double ppx = event[i].px();
            double ppy = event[i].py();
            double ppz = event[i].pz();
            double pe = event[i].e();
            double pm = event[i].m();
            int color = event[i].col();
            int acolor = event[i].acol();
            int status = event[i].status();
            if (event[i].isHadron()) {
                pidhadron[nhad] = pid;
                pxhadron[nhad] = ppx;
                pyhadron[nhad] = ppy;
                pzhadron[nhad] = ppz;
                pehadron[nhad] = pe;
                pmhadron[nhad] = pm;
                thadron[nhad] = 0.0;
                xhadron[nhad] = 0.0;
                yhadron[nhad] = 0.0;
                zhadron[nhad] = 0.0;
                nhad += 1;
            }


            if (event[i].isFinalPartonLevel()) {
//                if (abs(pid) == 111 || fabs(pid) == 311 || fabs(pid) == 215 || fabs(pid) == 213 || fabs(pid) == 211 || fabs(pid) == 323 || fabs(pid) == 321 || fabs(pid) == 2212 || fabs(pid) == 325 || fabs(pid) == 329 || fabs(pid) == 3222 || fabs(pid) == 3312 || fabs(pid) == 3334 || fabs(pid) == 411 || fabs(pid) == 413 || fabs(pid) == 415 || fabs(pid) == 431 || fabs(pid) == 433 || fabs(pid) == 435 || fabs(pid) == 521 || fabs(pid) == 523 || fabs(pid) == 525 || fabs(pid) == 541 || fabs(pid) == 543 || fabs(pid) == 545 || fabs(pid) == 4122 || fabs(pid) == 4222 || fabs(pid) == 5112 || fabs(pid) == 5222) {
//             if (fabs(pid) == 2212 || fabs(pid) == 2112 || fabs(pid) == 321 || fabs(pid) == 111 || fabs(pid) == 311 || fabs(pid) == 211) {
//                  if (fabs(pid) == 1 || fabs(pid) == 2 || fabs(pid) == 3 || fabs(pid) == 21 || fabs(pid) == 4 || fabs(pid) == 5) {


//                if (event[i].status() > 0) {
                    npar += 1;
                    //......formation time +
                    int ishower=0;

                    double p0[4]={0.0};        
                    double p4[4]={0.0};

                    double qt;
                    double timeplus=0.0;
                    double timeplus_0=0.0;
                    double timeplus_0_1=0.0;

                    int IDmom1, IDmom2;
                    int timebreaker = 0;

                    int IDmom;
                    int IDmom0 = i;

                    //cout<<"------------------------------------"<<" "<<i<<endl;

                    while(timebreaker == 0){
                        int IDiii = IDmom0;
                        if (abs(event[IDiii].status())==23 || abs(event[IDiii].status())==21 || abs(event[IDiii].status())==12) {
                            ishower=42;
                            timebreaker=1;			
                        }

                        IDmom1=event[IDiii].mother1();
                        IDmom2=event[IDiii].mother2();

                        if (IDmom1==IDmom2 && IDmom1==0) {
                            //cout<<"IDmom1==IDmom2 && IDmom1==0 timebreak"<<endl;
                            ishower=40;
                            timebreaker=1;
				}//if(IDmom1==IDmom2 && IDmom1==0)

				if(IDmom1==IDmom2 && IDmom1>0) {
				    IDmom0=IDmom1;
				}//if(IDmom1==IDmom2 && IDmom1>0)

				if(IDmom1>0 && IDmom2==0) {
				    IDmom0=IDmom1;
				    //........................................			  
				    double IDdaughter1=event[IDmom0].daughter1();
				    double IDdaughter2=event[IDmom0].daughter2();

				    if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){
					p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
					p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
					p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
					p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();

					double x_split=event[IDiii].e()/p4[0];

					//double x_split=event[IDiii].e()/event[IDmom0].e();

					if (x_split>1) x_split=1.0/x_split;

					p0[0]=event[IDiii].e();
					p0[1]=event[IDiii].px();
					p0[2]=event[IDiii].py();
					p0[3]=event[IDiii].pz(); 

					//p4[0]=event[IDmom0].e();
					//p4[1]=event[IDmom0].px();
					//p4[2]=event[IDmom0].py();
					//p4[3]=event[IDmom0].pz();

					double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
					double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));			  

					//cout<<"daughter"<<" "<<p0[1]<<" "<<p0[2]<<" "<<p0[3]<<" "<<p0[0]<<" "<<endl;
					//cout<<"mother"<<" "<<p4[1]<<" "<<p4[2]<<" "<<p4[3]<<" "<<p4[0]<<" "<<endl;

					rotate(p4[1],p4[2],p4[3],p0,1);

					qt=sqrt(pow(p0[1],2)+pow(p0[2],2));

					rotate(p4[1],p4[2],p4[3],p0,-1);

					double kt_daughter=qt;

					double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

					if(x_split<0.5){
					    //if(kt_daughter > 0.0001 && Q2 > 1.0)
					    if(kt_daughter > 0.0001) {
						//timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
						timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
					    }
					}
				    }
				}//if(IDmom1>0 && IDmom2==0)

				if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0) {
				    if(event[IDmom1].e()>event[IDmom2].e()) {
					IDmom0=IDmom1;			
				    }
				    if(event[IDmom1].e()<=event[IDmom2].e()) {
					IDmom0=IDmom2;			
				    }

				    double IDdaughter1=event[IDmom0].daughter1();
				    double IDdaughter2=event[IDmom0].daughter2();

				    if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){

					p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
					p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
					p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
					p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();

					double x_split=event[IDiii].e()/p4[0];

					//double x_split=event[IDiii].e()/event[IDmom0].e();

					if(x_split>1) x_split=1.0/x_split;

					p0[0]=event[IDiii].e();
					p0[1]=event[IDiii].px();
					p0[2]=event[IDiii].py();
					p0[3]=event[IDiii].pz(); 

					//p4[0]=event[IDmom0].e();
					//p4[1]=event[IDmom0].px();
					//p4[2]=event[IDmom0].py();
					//p4[3]=event[IDmom0].pz();

					rotate(p4[1],p4[2],p4[3],p0,1);

					qt=sqrt(pow(p0[1],2)+pow(p0[2],2));

					rotate(p4[1],p4[2],p4[3],p0,-1);

					double kt_daughter=qt;

					double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

					if(x_split<0.5) {
					    //if(kt_daughter > 0.0001 && Q2 > 1.0)
					    if(kt_daughter > 0.0001) {
						//timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
						timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
					    }
					}
				    }
				}//if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0)

				if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0) {
				    if(event[IDmom1].e()>event[IDmom2].e()) {
					IDmom0=IDmom1;			
				    }
				    if(event[IDmom1].e()<=event[IDmom2].e()) {
					IDmom0=IDmom2;			
				    }			

				    double IDdaughter1=event[IDmom0].daughter1();
				    double IDdaughter2=event[IDmom0].daughter2();

				    if(IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0){

					p4[0]=event[IDdaughter1].e()+event[IDdaughter2].e();
					p4[1]=event[IDdaughter1].px()+event[IDdaughter2].px();
					p4[2]=event[IDdaughter1].py()+event[IDdaughter2].py();
					p4[3]=event[IDdaughter1].pz()+event[IDdaughter2].pz();

					double x_split=event[IDiii].e()/p4[0];

					//double x_split=event[IDiii].e()/event[IDmom0].e();

					if(x_split>1) x_split=1.0/x_split;

					p0[0]=event[IDiii].e();
					p0[1]=event[IDiii].px();
					p0[2]=event[IDiii].py();
					p0[3]=event[IDiii].pz(); 

					//p4[0]=event[IDmom0].e();
					//p4[1]=event[IDmom0].px();
					//p4[2]=event[IDmom0].py();
					//p4[3]=event[IDmom0].pz();

					rotate(p4[1],p4[2],p4[3],p0,1);

					qt=sqrt(pow(p0[1],2)+pow(p0[2],2));

					rotate(p4[1],p4[2],p4[3],p0,-1);

					double kt_daughter=qt;

					double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

					if(x_split<0.5) {
					    //if(kt_daughter > 0.0001 && Q2 > 1.0)
					    if(kt_daughter > 0.0001) {						
						//timeplus=timeplus+2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
						timeplus=timeplus+2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;					
					    }
					}
				    }
				}//if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0)		


				if(abs(event[IDmom0].status())==23 || abs(event[IDmom0].status())==21 || abs(event[IDmom0].status())==12) {

				    // cout<<"event[IDmom0].status()"<<" "<<event[IDmom0].status()<<endl;

				    ishower=abs(event[IDmom0].status());

				    if(abs(event[IDmom0].status()) == 23)
				    {
					//ishower=23;
					//cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++23"<<endl;  
				    }

				    double x_split=event[IDiii].e()/event[IDmom0].e();

				    if(x_split>1) x_split=0;

				    p0[0]=event[IDiii].e();
				    p0[1]=event[IDiii].px();
				    p0[2]=event[IDiii].py();
				    p0[3]=event[IDiii].pz(); 

				    p4[0]=event[IDmom0].e();
				    p4[1]=event[IDmom0].px();
				    p4[2]=event[IDmom0].py();
				    p4[3]=event[IDmom0].pz();

				    rotate(p4[1],p4[2],p4[3],p0,1);

				    qt=sqrt(pow(p0[1],2)+pow(p0[2],2));

				    rotate(p4[1],p4[2],p4[3],p0,-1);

				    double kt_daughter=qt;

				    double Q2=1.0/(x_split*(1-x_split)/pow(kt_daughter,2));

				    if(kt_daughter > 0.0001 && Q2 > 1.0)
				    {

					timeplus_0=2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;

				    }

				    timeplus_0_1=2.0*event[IDmom0].e()*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;

				    timebreaker=1;			
				}

			    }//while(timebreaker == 0)
/*			    
			       mother1 = mother2 = 0: for lines 0 - 2, where line 0 represents the event as a whole, and 1 and 2 the two incoming beam particles;
			       mother1 = mother2 > 0: the particle is a "carbon copy" of its mother, but with changed momentum as a "recoil" effect, e.g. in a shower;
			       mother1 > 0, mother2 = 0: the "normal" mother case, where it is meaningful to speak of one single mother to several products, in a shower or decay;
			       mother1 < mother2, both > 0, for abs(status) = 81 - 86: primary hadrons produced from the fragmentation of a string spanning the range from mother1 to mother2, so that all partons in this range should be considered mothers; and analogously for abs(status) = 101 - 106, the formation of R-hadrons;
			       mother1 < mother2, both > 0, except case 4: particles with two truly different mothers, in particular the particles emerging from a hard 2 → n interaction.
			       mother2 < mother1, both > 0: particles with two truly different mothers, notably for the special case that two nearby partons are joined together into a status 73 or 74 new parton, in the g + q → q case the q is made first mother to simplify flavour tracing.
			    */   


			    //if(abs(ishower) != 23) timeplus=0.0;


			    //......formation time -				
			    //          cout << pid << "  " << pe << " " << ppx << " " << ppy << "  " << ppz << " " << peta << endl;
			    pidwant[irecord] = pid;
			    pxwant[irecord] = ppx;
			    pywant[irecord] = ppy;
			    pzwant[irecord] = ppz;
			    pewant[irecord] = pe;
			    pmwant[irecord] = pm;
			    xwant[irecord] = 0.0;
			    ywant[irecord] = 0.0;
			    zwant[irecord] = 0.0;
			    timewant[irecord] = timeplus;
                            cwant[irecord] = color;
                            acwant[irecord] = acolor;
                            statuswant[irecord] = status;
			    irecord = irecord + 1;
	                }
//			} 
//		    } // if (isparton)
		} // for (i in event)


                // mark the colour which is paired
                int unindex=0, unaindex=0, used[10000]={0}, aused[10000]={0}, uncol[1000]={0}, unacol[1000]={0};
                for (int ll=1;ll<=npar;++ll) {
                    for (int tt=1;tt<=npar;++tt) {
                        if ((cwant[ll]!=0) && (acwant[tt]!=0)) {
                            if (cwant[ll]==acwant[tt]) used[ll]=used[ll]+1, aused[tt]=aused[tt]+1;
            
                        }// if if ((nncol[ll]!=0) && (aacol[tt]!=0)
                    } // for (int tt=0;tt<Npart;++tt)      
                } // for (int ll=0;ll<Npart;++ll)

                // pick up the unpaird colour, put them in uncol and unacol
                for (int iD=1;iD<=npar;++iD) {
                    if ((used[iD]==0)&&(cwant[iD]!=0)) uncol[unindex]=cwant[iD], ++unindex;
                    if ((aused[iD]==0)&&(acwant[iD]!=0)) unacol[unaindex]=acwant[iD], ++unaindex;
                }


        for (int i=0;i<unindex;++i) {
//            cout << uncol[i] << endl;

            pidwant[npar+1]=-3;pxwant[npar+1]=0.10;pywant[npar+1]=0.20;pzwant[npar+1]=-0.10;
            pewant[npar+1]=sqrt(pzwant[npar+1]*pzwant[npar+1]+pywant[npar+1]*pywant[npar+1]+pxwant[npar+1]*pxwant[npar+1]+0.5*0.5);
            xwant[npar+1]=0.0;ywant[npar+1]=0.0;zwant[npar+1]=0.0;timewant[npar+1]=100.0;
            used[npar+1]=1;
            cwant[npar+1]=0;
            acwant[npar+1]=uncol[i];
            pmwant[npar+1]=0.5;
            statuswant[npar+1]=62;
            npar=npar+1;
        }
//        cout << "the unpaired anti-colour:" << endl;
        for (int i=0;i<unaindex;++i) {
//            cout << unacol[i] << endl;
            pidwant[npar+1]=3;pxwant[npar+1]=0.10;pywant[npar+1]=0.20;pzwant[npar+1]=0.10;           
            pewant[npar+1]=sqrt(pzwant[npar+1]*pzwant[npar+1]+pywant[npar+1]*pywant[npar+1]+pxwant[npar+1]*pxwant[npar+1]+0.5*0.5);            
            xwant[npar+1]=0.0;ywant[npar+1]=0.0;zwant[npar+1]=0.0;timewant[npar+1]=100.0;
            aused[npar+1]=1;
            cwant[npar+1]=unacol[i];
            acwant[npar+1]=0;
            pmwant[npar+1]=0.5;
            statuswant[npar+1]=62;
            npar=npar+1;
        }

                f_hadron << numEvent+1 << "  " << nhad << endl;
                for(int iD=0; iD<nhad; iD++)
                {                       f_hadron << pidhadron[iD] << "  " << pxhadron[iD] << "  " << pyhadron[iD] << "  " << pzhadron[iD] << "  " << pehadron[iD] << "  " << pmhadron[iD] << "  " << xhadron[iD]  << "  " << yhadron[iD] << "  " << zhadron[iD] << "  " << thadron[iD] << endl;
                }

		f_pp << numEvent+1 << "  " << npar << endl;
		for(int iD=1; iD<=npar; iD++) 
		{
                    f_pp << iD << "  " << pidwant[iD] << "  " << pxwant[iD] << "  " << pywant[iD] << "  " << pzwant[iD] << "  " << pewant[iD] << "  " << pmwant[iD] << "  " << xwant[iD]  << "  " << ywant[iD] << "  " << zwant[iD] << "  " << timewant[iD] << "  " << cwant[iD] << "  " << acwant[iD] << endl;
		}
	//*************************************************************************************************/


		//  parton information
		sigma += pythia.info.sigmaGen();
		numEvent++;
		if (numEvent == nEvent) 
		{
		    cout << n << "  " << numEvent  << "  " << sigma << endl;
                    f_info << "# N_total" << "  " << "N_selected" << "  " << "Simga_seclected" << endl;
		    f_info << n << "  " << numEvent << "  " << sigma << endl;
                    f_info << "# N_total: Total number of events generated by PYTHIA."<<endl;
                    f_info << "# N_selected: Number of events that pass the selection criteria."<< endl;
                    f_info << "# Simga_seclected: Sum of the cross sections of the selected events."<< endl;
		}

    } // event loop ends


    f_hadron.close();
    f_pp.close();
    f_info.close();

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
