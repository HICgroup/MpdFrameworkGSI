#define _MAX_TRACKS 200000
#define _N2_ZDC 225
#define _N_ZDC 15
#define _N_MODULES_TOTAL 90
#define _N_ARM 2
#define _N_HARM 2
#define _N_METHOD 2
#define _N_QCOMP 2
#define UNDEFINED_CENTRALITY -9999
#define UNDEFINED_PID -9999

#include <iostream>

#include <TMath.h>
#include <TSystem.h>
#include <TVector3.h>
//#include "TRoot.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "FairMCEventHeader.h"
//#include "MCHeader.h"
#include "MpdEvent.h"
//#include "MPDEvent.h"
#include "MpdZdcDigi.h"
#include "MpdPid.h"
//#include "MpdPid_AZ.h"
//#include "ZDCHit.h"
#include "TClonesArray.h"
#include "MpdTrack.h"
#include "FairMCTrack.h"
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom1.h>
#include <MpdKalmanTrack.h> 
#include <MpdVertex.h>

using std::cout;
using std::endl;
using TMath::ATan2;

#define	NmultiplicityBins 100

const Float_t pt_bins[]={0.,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.};
const Int_t   n_pt_bin = 12;

const Float_t eta_bins[]={-1.5,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.5};
const Int_t   n_eta_bin = 14;
const Int_t   n_proj    = 3;

const Int_t   n_hits_cut= 32;

const Double_t PIDsigM 	 	= 4.0;
const Double_t PIDsigE 	 	= 4.0;
const Double_t PIDenergy 	= 11.;
const Double_t PIDkoeff	 	= 1.;
const TString  PIDgenerator = "URQMD";
const TString  PIDtracking  = "CF";
const TString  PIDparticles = "pikapr";

class reducedTreeCreator
{
	public:
		reducedTreeCreator(TString inFileHistName, TString inFileTreeName, TString outFileName , TString dcaFileName);
		void CreateReducedTree();
		Int_t GetCentrality(Int_t multiplicity);
		Float_t integrate(TH1F *h, Int_t max_bin, Float_t sum);
		bool FillTrack(FairMCTrack *mctrack, long int j, int &m);

	private:
	
		TFile *inFile;
		TTree *inTree;
	
		TFile  *outFile;
		TTree  *outTree;
		
		TFile *inFileHist;

		long int mc_side[_MAX_TRACKS][10];
		long int mpd_side[_MAX_TRACKS];

		MpdPid *pid;

		TFile* dcaFile;
		//TF1*   f_dca[n_proj][n_pt_bin][n_eta_bin];
		TF1* f_pt_fit[n_proj][n_eta_bin];

		Float_t b_mc;
		Float_t phiEP_mc;
		Float_t x_vertex_mc;
		Float_t y_vertex_mc;
		Float_t z_vertex_mc;
		Float_t x_vertex_mpd;
		Float_t y_vertex_mpd;
		Float_t z_vertex_mpd;
		Long_t n_tracks_mc;
		Float_t eta_mc[_MAX_TRACKS];
		Float_t pt_mc[_MAX_TRACKS];
		Int_t mother_ID_mc[_MAX_TRACKS];
		Int_t PDG_code_mc[_MAX_TRACKS];
		Float_t px_mc[_MAX_TRACKS];
		Float_t py_mc[_MAX_TRACKS];
		Float_t pz_mc[_MAX_TRACKS];
		Float_t start_x_mc[_MAX_TRACKS];
		Float_t start_y_mc[_MAX_TRACKS];
		Float_t start_z_mc[_MAX_TRACKS];
		Float_t mass_mc[_MAX_TRACKS];
		Float_t energy_mc[_MAX_TRACKS];
		
		Long_t n_tracks_mpd;
		Long_t k_tracks_mpd;
		Float_t eta_mpd[_MAX_TRACKS];
		Float_t phi_mpd[_MAX_TRACKS];
		Float_t theta_mpd[_MAX_TRACKS];
		Int_t TOF_flag_mpd[_MAX_TRACKS];
		Float_t ZDC_energy_mpd[_N_MODULES_TOTAL];
		Float_t pid_tpc_prob_electron_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_pion_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_kaon_mpd[_MAX_TRACKS];
		Float_t pid_tpc_prob_proton_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_electron_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_pion_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_kaon_mpd[_MAX_TRACKS];
		Float_t pid_tof_prob_proton_mpd[_MAX_TRACKS];
		Float_t tof_beta_mpd[_MAX_TRACKS];
		Float_t tof_mass2_mpd[_MAX_TRACKS];
		Float_t dEdx_tpc_mpd[_MAX_TRACKS];
		Float_t chi2_mpd[_MAX_TRACKS];
		Float_t chi2_vertex[_MAX_TRACKS];
		Float_t pt_error_mpd[_MAX_TRACKS];
		Float_t theta_error_mpd[_MAX_TRACKS];
		Float_t phi_error_mpd[_MAX_TRACKS];
		Float_t DCA_x_mpd[_MAX_TRACKS];
		Float_t DCA_y_mpd[_MAX_TRACKS];
		Float_t DCA_z_mpd[_MAX_TRACKS];
		Int_t n_hits_mpd[_MAX_TRACKS];
		Int_t n_hits_poss_mpd[_MAX_TRACKS];
		Float_t signed_pt_mpd[_MAX_TRACKS];
		Float_t p_mpd[_MAX_TRACKS];
		Int_t centrality_tpc_mpd;

		
		FairMCEventHeader *MCHeader;
		TClonesArray *MCTracks;
		MpdEvent *MPDEvent;
		TClonesArray *ZDCHits;
		TClonesArray *MpdGlobalTracks;
		MpdZdcDigi* ZDCHit;
		TClonesArray *mpdKalmanTracks;
		TClonesArray *vertexes;

		TRandom *RNG = new TRandom();
		Float_t multiplicity_bins[NmultiplicityBins+1];
		TH1F *h_multiplicity_before;

		ClassDef(reducedTreeCreator,1);
};
