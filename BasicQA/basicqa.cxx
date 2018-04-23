#define basicqa_cxx
#include "basicqa.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <TRandom3.h>
#include <TDatime.h>

#include <iostream>
#include <stdlib.h> 

std::pair<Double_t, Double_t> GetSmearedPair(std::pair<Double_t, Double_t> Center){
  std::pair<Double_t, Double_t> pair;
  TDatime *  t = new TDatime();
  TRandom3 * r = new TRandom3(t->GetDate()*t->GetTime());
  
  pair.first  = (r->Rndm()-0.5)*15 + Center.first;
  pair.second = (r->Rndm()-0.5)*15 + Center.second;
  
  return pair;
}

std::pair<Double_t, Double_t> GetModuleCoordinates(Int_t Nmodule){
  // 0-44 - right FHCal; 45 -90 - left FHCal;
  std::pair<Double_t, Double_t> pXY;
  
  if      (Nmodule>=0 && Nmodule<=4)    { pXY.first =  15*(Nmodule-2);          pXY.second = 45.; }
  else if (Nmodule>=5 && Nmodule<=39)   { pXY.first =  15*(3-(Nmodule+2)%7);    pXY.second = 15*(3-(Nmodule+2)/7); }
  else if (Nmodule>=40 && Nmodule<=44)  { pXY.first =  15*(Nmodule-42);         pXY.second = -45.; }
  else if (Nmodule>=45 && Nmodule<=49)  { pXY.first = -15*(Nmodule-45-2);       pXY.second = 45.; }
  else if (Nmodule>=50 && Nmodule<=84)  { pXY.first = -15*(3-(Nmodule-45+2)%7); pXY.second = 15*(3-(Nmodule-45+2)/7); }
  else if (Nmodule>=85 && Nmodule<=89)  { pXY.first = -15*(Nmodule-45-42);      pXY.second = -45.; }
  
  return pXY;
}

void basicqa::LoopPID()
{
   std::cout << "\nPID QA started." << std::endl;
   
   if (fChain == 0) return;

   const Double_t PIDptBins  []  	 = {0.,0.2,0.5,1.,1.5,2.,3.};
   const Int_t    PIDNpt        	 = 6;
   const Double_t PIDetaBins [] 	 = {-1.5,-1.,0.,1.,1.5};
   const Int_t    PIDNeta       	 = 4;
   		 TString  PIDparticleBins [] = {TString("pion"), TString("kaon"), TString("proton")}; 
   const Int_t    PIDNparticles  	 = 3;
   const Int_t    PIDNhist		     = 300;
   const Double_t PIDMass2 []		 = {0.1396*0.1396,0.4937*0.4937,0.9383*0.9383};
   const Double_t PIDMass2Width []	 = {PIDMass2[0]*3.,PIDMass2[1]*3.,PIDMass2[2]*3.};
   const Int_t    PIDpdg [] 		 = {211, 321, 2212};
   const Int_t    PIDTofFlags [] 	 = {0,2,4,6};
   const Int_t    PIDNTofFlag  		 = 4; 
   const Int_t    Nmodules		 = 90;

   TH1F* hPIDm2[PIDNparticles][PIDNpt][PIDNeta];
   TH1F* hPIDm2Before[PIDNparticles][PIDNpt][PIDNeta];
   TH1F* hPIDm2PDGpion[PIDNparticles][PIDNpt][PIDNeta];
   TH1F* hPIDm2PDGkaon[PIDNparticles][PIDNpt][PIDNeta];
   TH1F* hPIDm2PDGproton[PIDNparticles][PIDNpt][PIDNeta];
   TH1F* hPIDm2All[PIDNpt][PIDNeta];
   TH1F* hPIDtofFlagDiff[PIDNparticles][PIDNpt][PIDNeta];
   TH2F* hPIDtofFlagEtaDiff[PIDNparticles][PIDNpt];
   TH2F* hPIDtofFlagPtDiff[PIDNparticles][PIDNeta];
   TH1F* hPIDtofFlag = new TH1F("hPIDtofFlag","TOF flag;TOFflag;N_{count}",10,0,10);
   TH2F* hPIDtofFlagEta = new TH2F("hPIDtofFlagEta","TOF flag vs #eta (0 < p_{T} < 3 [GeV/c]);#eta;TOF flag"	   ,100,-1.5,1.5,10,0,10);
   TH2F* hPIDtofFlagPt  = new TH2F("hPIDtofFlagPt","TOF flag vs p_{T} (-1.5 < #eta < 1.5);p_{T}, [GeV/c];TOF flag",100,0.,3.,10,0,10);

   TH1F* hEFFptBefore   = new TH1F("hEFFptBefore","pt before TOF flag selection",100,0.,3.); hEFFptBefore -> Sumw2();
   TH1F* hEFFptAfter    = new TH1F("hEFFptAfter" ,"pt after TOF flag selection",100,0.,3.); hEFFptAfter -> Sumw2();

   TH1F* hPIDProbability[PIDNparticles];

   TH2F* hPIDAcceptance[PIDNparticles][PIDNTofFlag];
   TH2F* hPIDAcceptancePDG[PIDNparticles][PIDNTofFlag];
   TH1F* hFHCalEtot[Nmodules];
   TH2F* hFHCalEtotPosR = new TH2F("hFHCalEtotPosR","E_{tot} in R FHCal",7,-52.5,52.5,7,-52.5,52.5);
   TH2F* hFHCalEtotPosL = new TH2F("hFHCalEtotPosL","E_{tot} in L FHCal",7,-52.5,52.5,7,-52.5,52.5);

   TH2F* hPIDm2MomTotal = new TH2F("hPIDm2MomTotal","m^{2} vs momentum",220,0.,11.,400,-2.,2.);
   TH2F* hPIDm2Mom[PIDNparticles];
   
   for (int ModuleIterator=0;ModuleIterator<Nmodules;ModuleIterator++){
      std::string nameEtot = "hFHCalEtot" + std::to_string(ModuleIterator);
      char titleEtot[200];
      sprintf(titleEtot, "Energy deposited in the %i module;E_{tot}, [a.u.];N_{count}", ModuleIterator);
      hFHCalEtot[ModuleIterator] = new TH1F(nameEtot.c_str(),titleEtot,300,0.,30.);
   }

   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){

   std::string nameMom = "hPIDm2Mom" + std::to_string(ParticleIterator);
	   char titleMom[200];
	   sprintf(titleMom, "m^{2} vs momentum of %s (realistic PID);p, [GeV/c];m^{2}, [GeV/c^{2}]^{2}", PIDparticleBins[ParticleIterator].Data());
	   hPIDm2Mom[ParticleIterator] = new TH2F(nameMom.c_str(),titleMom, 220, 0., 11., 400, -2., 2.);

	std::string nameProb = "hPIDProbability" + std::to_string(ParticleIterator);
	   char titleProb[200];
	   sprintf(titleProb, "Probability of %s (with TOF hits);Probablility;N_{count}", PIDparticleBins[ParticleIterator].Data());
	   hPIDProbability[ParticleIterator] = new TH1F(nameProb.c_str(),titleProb, 100, 0., 1.);

	 for (int TofFlagIterator=0; TofFlagIterator < PIDNTofFlag; TofFlagIterator++){
	   std::string nameAcc = "hPIDAcceptance" + std::to_string(ParticleIterator) + std::to_string(TofFlagIterator);
	   char titleAcc[200];
	   sprintf(titleAcc, "p_{T} vs #eta of %s (real PID) at TOF Flag %i;#eta;p_{T}, [GeV/c]", PIDparticleBins[ParticleIterator].Data(), PIDTofFlags[TofFlagIterator]);
	   hPIDAcceptance[ParticleIterator][TofFlagIterator] = new TH2F(nameAcc.c_str(),titleAcc, 100, -1.5, 1.5, 100, 0., 3.);

	   std::string nameAccPDG = "hPIDAcceptancePDG" + std::to_string(ParticleIterator) + std::to_string(TofFlagIterator);
	   char titleAccPDG[200];
	   sprintf(titleAccPDG, "p_{T} vs #eta of %s (PDG) at TOF Flag %i;#eta;p_{T}, [GeV/c]", PIDparticleBins[ParticleIterator].Data(), PIDTofFlags[TofFlagIterator]);
	   hPIDAcceptancePDG[ParticleIterator][TofFlagIterator] = new TH2F(nameAccPDG.c_str(),titleAccPDG, 100, -1.5, 1.5, 100, 0., 3.);
	 }

	 for (int PtIterator=0; PtIterator < PIDNpt; PtIterator++){
	   for (int EtaIterator=0;EtaIterator < PIDNeta; EtaIterator++){
		 std::string name = "hPIDmSqr" + std::to_string(ParticleIterator) + std::to_string(PtIterator) + std::to_string(EtaIterator);
		 char title[200];
		 sprintf(title, "m^{2} of %s at (%1.1f < p_{T} < %1.1f [GeV/c]) (%1.1f < #eta < %1.1f);m^{2}, [GeV/c]^{2};N_{count}", PIDparticleBins[ParticleIterator].Data(),
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1],PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 hPIDm2[ParticleIterator][PtIterator][EtaIterator] = new TH1F(name.c_str(),title,PIDNhist,PIDMass2[ParticleIterator]-PIDMass2Width[ParticleIterator],
																	  PIDMass2[ParticleIterator]+PIDMass2Width[ParticleIterator]);

		 std::string nameBefore = "hPIDmSqrBeforeCut" + std::to_string(ParticleIterator) + std::to_string(PtIterator) + std::to_string(EtaIterator);
		 char titleBefore[200];
		 sprintf(titleBefore, "m^{2} before PID at (%1.1f < p_{T} < %1.1f [GeV/c]) (%1.1f < #eta < %1.1f);m^{2}, [GeV/c]^{2};N_{count}",
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1],PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 hPIDm2Before[ParticleIterator][PtIterator][EtaIterator] = new TH1F(nameBefore.c_str(),titleBefore,PIDNhist,PIDMass2[ParticleIterator]-PIDMass2Width[ParticleIterator],
																	  PIDMass2[ParticleIterator]+PIDMass2Width[ParticleIterator]);

		 std::string namePDGpi = "hPIDmSqrPDGpion" + std::to_string(ParticleIterator) + std::to_string(PtIterator) + std::to_string(EtaIterator);
		 char titlePDGpi[200];
		 sprintf(titlePDGpi, "m^{2} (pion PDG) of %s at (%1.1f < p_{T} < %1.1f [GeV/c]) (%1.1f < #eta < %1.1f);m^{2}, [GeV/c]^{2};N_{count}", PIDparticleBins[ParticleIterator].Data(),
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1],PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 hPIDm2PDGpion[ParticleIterator][PtIterator][EtaIterator] = new TH1F(namePDGpi.c_str(),titlePDGpi,PIDNhist,PIDMass2[ParticleIterator]-PIDMass2Width[ParticleIterator],
																	  PIDMass2[ParticleIterator]+PIDMass2Width[ParticleIterator]);

		 std::string namePDGka = "hPIDmSqrPDGkaon" + std::to_string(ParticleIterator) + std::to_string(PtIterator) + std::to_string(EtaIterator);
		 char titlePDGka[200];
		 sprintf(titlePDGka, "m^{2} (kaon PDG) of %s at (%1.1f < p_{T} < %1.1f [GeV/c]) (%1.1f < #eta < %1.1f);m^{2}, [GeV/c]^{2};N_{count}", PIDparticleBins[ParticleIterator].Data(),
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1],PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 hPIDm2PDGkaon[ParticleIterator][PtIterator][EtaIterator] = new TH1F(namePDGka.c_str(),titlePDGka,PIDNhist,PIDMass2[ParticleIterator]-PIDMass2Width[ParticleIterator],
																	  PIDMass2[ParticleIterator]+PIDMass2Width[ParticleIterator]);

		 std::string namePDGpr = "hPIDmSqrPDGproton" + std::to_string(ParticleIterator) + std::to_string(PtIterator) + std::to_string(EtaIterator);
		 char titlePDGpr[200];
		 sprintf(titlePDGpr, "m^{2} (proton PDG) of %s at (%1.1f < p_{T} < %1.1f [GeV/c]) (%1.1f < #eta < %1.1f);m^{2}, [GeV/c]^{2};N_{count}", PIDparticleBins[ParticleIterator].Data(),
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1],PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 hPIDm2PDGproton[ParticleIterator][PtIterator][EtaIterator] = new TH1F(namePDGpr.c_str(),titlePDGpr,PIDNhist,PIDMass2[ParticleIterator]-PIDMass2Width[ParticleIterator],
																	  PIDMass2[ParticleIterator]+PIDMass2Width[ParticleIterator]);

		 std::string nameTOF = "hPIDtofFlagDiff" + std::to_string(ParticleIterator) + std::to_string(PtIterator) + std::to_string(EtaIterator);
		 char titleTOF[200];
		 sprintf(titleTOF, "TOF flag of %s at (%1.1f < p_{T} < %1.1f [GeV/c]) (%1.1f < #eta < %1.1f);TOF flag;N_{count}", PIDparticleBins[ParticleIterator].Data(),
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1],PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 hPIDtofFlagDiff[ParticleIterator][PtIterator][EtaIterator] = new TH1F(nameTOF.c_str(),titleTOF,10,0,10);

		 if (EtaIterator == 0){
			std::string nameTOFEta = "hPIDtofFlagEtaDiff" + std::to_string(ParticleIterator) + std::to_string(PtIterator);
		 	char titleTOFEta[200];
		 	sprintf(titleTOFEta, "TOF flag vs #eta of %s at (%1.1f < p_{T} < %1.1f [GeV/c]);TOF flag;p_{T}, [GeV/c]", PIDparticleBins[ParticleIterator].Data(),
				 PIDptBins[PtIterator],PIDptBins[PtIterator+1]);
		 	hPIDtofFlagEtaDiff[ParticleIterator][PtIterator] = new TH2F(nameTOFEta.c_str(),titleTOFEta,100,-1.5,1.5,10,0,10);
		 }

		 if (PtIterator == 0){
			std::string nameTOFPt = "hPIDtofFlagPtDiff" + std::to_string(ParticleIterator) + std::to_string(EtaIterator);
		 	char titleTOFPt[200];
		 	sprintf(titleTOFPt, "TOF flag vs #eta of %s at (%1.1f < #eta < %1.1f);TOF flag;#eta", PIDparticleBins[ParticleIterator].Data(),
				 PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		 	hPIDtofFlagPtDiff[ParticleIterator][EtaIterator] = new TH2F(nameTOFPt.c_str(),titleTOFPt,100,0.,3.,10,0,10);
		 }

		 if (ParticleIterator == 0){
		   std::string name = "hPIDmSqrTotal" + std::to_string(PtIterator) + std::to_string(EtaIterator);
		   char title[200];
		   sprintf(title, "Total m^{2} at (%1.1f < p_{T} < %1.1f) (%1.1f < #eta < %1.1f);m^{2}, [GeV/c]^{2};N_{count}", PIDptBins[PtIterator],PIDptBins[PtIterator+1],
				   PIDetaBins[EtaIterator],PIDetaBins[EtaIterator+1]);
		   hPIDm2All[PtIterator][EtaIterator] = new TH1F(name.c_str(),title,PIDNhist*10,-1.,10.);
		 }
	   }
	 }
   }

   Long64_t nentries = fChain->GetEntriesFast();
   Int_t Npt, Neta, Npid, Ntof;
   Bool_t isPion, isKaon, isProton;
   std::pair<Double_t, Double_t> FHcalXY;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	  if (jentry % 100 == 0) std::cout << "Event N: " << jentry << std::endl;
	  
	  for (int ModuleIterator = 0; ModuleIterator<Nmodules;ModuleIterator++){
	    hFHCalEtot[ModuleIterator] -> Fill(ZDC_energy_mpd[ModuleIterator]);
	    FHcalXY = GetModuleCoordinates(ModuleIterator);
	    if (ModuleIterator < 45) hFHCalEtotPosR -> Fill(FHcalXY.first, FHcalXY.second, ZDC_energy_mpd[ModuleIterator]);
	    else 		     hFHCalEtotPosL -> Fill(FHcalXY.first, FHcalXY.second, ZDC_energy_mpd[ModuleIterator]);
	  }

	  for (Long64_t jtrack=0; jtrack<n_tracks_mpd; jtrack++){
		Npt = -1; Neta = -1; Npid = -1; Ntof = -1;
		isPion = kFALSE; isKaon = kFALSE; isProton = kFALSE;
		for (int PtIterator=0; PtIterator < PIDNpt; PtIterator++){
		  if (TMath::Abs(signed_pt_mpd[jtrack]) >= PIDptBins[PtIterator] &&
			  TMath::Abs(signed_pt_mpd[jtrack]) <  PIDptBins[PtIterator+1]){
			Npt = PtIterator;
		  }
		}
		for (int EtaIterator=0;EtaIterator < PIDNeta; EtaIterator++){
		  if (eta_mpd[jtrack] >= PIDetaBins[EtaIterator] && eta_mpd[jtrack] < PIDetaBins[EtaIterator+1]){
			Neta = EtaIterator;
		  }
		}

		isPion = (pid_tpc_prob_pion_mpd[jtrack] > 0.9 && pid_tpc_prob_kaon_mpd[jtrack] < 0.1 &&
				  pid_tpc_prob_proton_mpd[jtrack] < 0.1);
		isKaon = (pid_tpc_prob_pion_mpd[jtrack] < 0.1 && pid_tpc_prob_kaon_mpd[jtrack] > 0.9 &&
				  pid_tpc_prob_proton_mpd[jtrack] < 0.1);
		isProton = (pid_tpc_prob_pion_mpd[jtrack] < 0.1 && pid_tpc_prob_kaon_mpd[jtrack] < 0.1 &&
				  pid_tpc_prob_proton_mpd[jtrack] > 0.9);

		for (int TofFlagIterator=0; TofFlagIterator < PIDNTofFlag; TofFlagIterator++){
			if (TOF_flag_mpd[jtrack] == PIDTofFlags[TofFlagIterator]) Ntof = TofFlagIterator;
		}

		if (Npt == -1 || Neta == -1) continue;

		if (TMath::Abs(eta_mpd[jtrack]) < 1.4){
			hEFFptBefore -> Fill(TMath::Abs(signed_pt_mpd[jtrack]));
			if (TOF_flag_mpd[jtrack] == 2 || TOF_flag_mpd[jtrack] == 6) hEFFptAfter -> Fill(TMath::Abs(signed_pt_mpd[jtrack]));
		}

		if (TOF_flag_mpd[jtrack] == 2 || TOF_flag_mpd[jtrack] == 6){
			hPIDProbability[0] -> Fill(pid_tpc_prob_pion_mpd[jtrack]);
			hPIDProbability[1] -> Fill(pid_tpc_prob_kaon_mpd[jtrack]);
			hPIDProbability[2] -> Fill(pid_tpc_prob_proton_mpd[jtrack]);
		}
		if (isPion) Npid = 0; if (isKaon) Npid = 1; if (isProton) Npid = 2;
		if (signed_pt_mpd[jtrack] > 0) continue;
		if (n_hits_mpd[jtrack] < 32) continue;

		if (TOF_flag_mpd[jtrack] == 2 || TOF_flag_mpd[jtrack] == 6) hPIDm2MomTotal -> Fill(p_mpd[jtrack], tof_mass2_mpd[jtrack]);
		
		for (int i=0;i<3;i++){
			if (TOF_flag_mpd[jtrack] == 2 || TOF_flag_mpd[jtrack] == 6)
				hPIDm2Before[i][Npt][Neta] -> Fill(tof_mass2_mpd[jtrack]);
		}

		if (Npid == -1) continue;

		if (TOF_flag_mpd[jtrack] == 2 || TOF_flag_mpd[jtrack] == 6) hPIDm2Mom[Npid] -> Fill(p_mpd[jtrack], tof_mass2_mpd[jtrack]);

		hPIDtofFlag -> Fill(TOF_flag_mpd[jtrack]);
		hPIDtofFlagEta -> Fill(eta_mpd[jtrack], TOF_flag_mpd[jtrack]);
		hPIDtofFlagPt  -> Fill(-signed_pt_mpd[jtrack], TOF_flag_mpd[jtrack]);
		hPIDtofFlagDiff[Npid][Npt][Neta] -> Fill(TOF_flag_mpd[jtrack]);

		hPIDtofFlagEtaDiff[Npid][Npt] -> Fill(eta_mpd[jtrack],TOF_flag_mpd[jtrack]);
		hPIDtofFlagPtDiff[Npid][Neta] -> Fill(-signed_pt_mpd[jtrack],TOF_flag_mpd[jtrack]);

		hPIDAcceptance[Npid][Ntof]    -> Fill(eta_mpd[jtrack],-signed_pt_mpd[jtrack]);
		if (PIDpdg[Npid] == PDG_code_mc[id_from_mc_mpd[jtrack]]) 
			hPIDAcceptancePDG[Npid][Ntof] -> Fill(eta_mpd[jtrack],-signed_pt_mpd[jtrack]);

		if (TOF_flag_mpd[jtrack] != 2 && TOF_flag_mpd[jtrack] != 6) continue;

		// std::cout << "Pt (" << Npt << "): " << TMath::Abs(signed_pt_mpd[jtrack]) << std::endl;
		// std::cout << "PID(" << Npid << "): " << pid_tpc_prob_pion_mpd[jtrack] << "|" << 
		//			 pid_tpc_prob_kaon_mpd[jtrack] << "|" << pid_tpc_prob_proton_mpd[jtrack] << std::endl;
		// std::cout << "Eta(" << Neta << "): " << eta_mpd[jtrack] << std::endl;

		hPIDm2[Npid][Npt][Neta] -> Fill(tof_mass2_mpd[jtrack]);
		hPIDm2All[Npt][Neta]    -> Fill(tof_mass2_mpd[jtrack]);
		
		if (PIDpdg[0] == PDG_code_mc[id_from_mc_mpd[jtrack]]) 
		  hPIDm2PDGpion[Npid][Npt][Neta] -> Fill(tof_mass2_mpd[jtrack]);
		if (PIDpdg[1] == PDG_code_mc[id_from_mc_mpd[jtrack]]) 
		  hPIDm2PDGkaon[Npid][Npt][Neta] -> Fill(tof_mass2_mpd[jtrack]);
		if (PIDpdg[2] == PDG_code_mc[id_from_mc_mpd[jtrack]]) 
		  hPIDm2PDGproton[Npid][Npt][Neta] -> Fill(tof_mass2_mpd[jtrack]);
		

	  }// end reco track loop

   }// end event loop

   std::cout << "Writing output." << std::endl;

   TFile* outFile = new TFile(outputFileName.Data(),"recreate");
   outFile -> mkdir("EnergyDeposition");
   outFile -> mkdir("Probability");
   outFile -> mkdir("Acceptance");
   outFile -> mkdir("tofFlag");
   outFile -> mkdir("mass2");

   outFile->cd();

   hEFFptBefore -> Write();
   hEFFptAfter 	-> Write();

   hPIDm2MomTotal -> Write();
   for (int ParticleIterator=0; ParticleIterator<PIDNparticles;ParticleIterator++){
		hPIDm2Mom[ParticleIterator] -> Write();
   }
   
   outFile -> cd("EnergyDeposition");
   hFHCalEtotPosR -> Write();
   hFHCalEtotPosL -> Write();
   for (int ModuleIterator=0; ModuleIterator<Nmodules; ModuleIterator++){
     hFHCalEtot[ModuleIterator] -> Write();
   }

   outFile -> cd("Probability");

   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){
	 hPIDProbability[ParticleIterator] -> Write();
   }

   outFile -> cd("Acceptance");

   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){
	 for (int TofFlagIterator=0; TofFlagIterator < PIDNTofFlag; TofFlagIterator++){
		hPIDAcceptance[ParticleIterator][TofFlagIterator]    -> Write();
		hPIDAcceptancePDG[ParticleIterator][TofFlagIterator] -> Write();
	 }
   }

   outFile -> cd ("tofFlag");

   hPIDtofFlag -> Write();
   hPIDtofFlagEta -> Write();
   hPIDtofFlagPt  -> Write();
   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){
	 for (int PtIterator=0; PtIterator < PIDNpt; PtIterator++){
	   for (int EtaIterator=0;EtaIterator < PIDNeta; EtaIterator++){
		 hPIDtofFlagDiff[ParticleIterator][PtIterator][EtaIterator] -> Write();
	   }
	 }
   }

   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){
	 for (int PtIterator=0; PtIterator < PIDNpt; PtIterator++){
		hPIDtofFlagEtaDiff[ParticleIterator][PtIterator] -> Write();
	 }
   }

   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){
	 for (int EtaIterator=0;EtaIterator < PIDNeta; EtaIterator++){
		hPIDtofFlagPtDiff[ParticleIterator][EtaIterator] -> Write();
	 }
   }

   outFile -> cd("mass2");
   
   for (int ParticleIterator=0; ParticleIterator < PIDNparticles; ParticleIterator++){
	 for (int PtIterator=0; PtIterator < PIDNpt; PtIterator++){
	   for (int EtaIterator=0;EtaIterator < PIDNeta; EtaIterator++){
		 hPIDm2[ParticleIterator][PtIterator][EtaIterator] -> Write();
		 hPIDm2Before[ParticleIterator][PtIterator][EtaIterator] -> Write();
		 hPIDm2PDGpion[ParticleIterator][PtIterator][EtaIterator] -> Write();
		 hPIDm2PDGkaon[ParticleIterator][PtIterator][EtaIterator] -> Write();
		 hPIDm2PDGproton[ParticleIterator][PtIterator][EtaIterator] -> Write();
	   }
	 }
   }
   for (int PtIterator=0; PtIterator < PIDNpt; PtIterator++){
	 for (int EtaIterator=0;EtaIterator < PIDNeta; EtaIterator++){
	   hPIDm2All[PtIterator][EtaIterator] -> Write();
	 }
   }
   outFile->Close();
}

void basicqa::LoopDCA()
{
   std::cout << "\nDCA QA started." << std::endl;
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	  if (jentry % 1000 == 0) std::cout << "Event N: " << jentry << std::endl;
   }
}

int main(int argc, char** argv){
  TString iFileName, oFileName;
  Bool_t  isPID, isDCA;

  isPID = kFALSE; isDCA = kFALSE;
  if (argc < 3){
	std::cerr << "./basicqa -i inputfile -o outputfile" << std::endl;
	return 1;
  }

  for (int i=1; i<argc; i++){
	if(std::string(argv[i]) != "-i" &&
	  std::string(argv[i]) != "-o" &&
	  std::string(argv[i]) != "--pid" &&
	  std::string(argv[i]) != "--dca"){
		std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " <<  argv[i] << std::endl;
		return 1;
	} else {
	  if(std::string(argv[i]) == "-i" && i!=argc-1) {
		iFileName = argv[++i];
		continue;
	  }
	  if(std::string(argv[i]) == "-i" && i==argc-1) {
		std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
		return 1;
	  }
	  if(std::string(argv[i]) == "-o" && i!=argc-1) {
		oFileName = argv[++i];
		continue;
	  }
	  if(std::string(argv[i]) == "-o" && i==argc-1){
		std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
		return 1;
	  }
	  if(std::string(argv[i]) == "--pid") {
		isPID = kTRUE;
		continue;
	  }
	  if(std::string(argv[i]) == "--dca") {
		isDCA = kTRUE;
		continue;
	  }
    }
  }

  if (iFileName == "" || oFileName == ""){
	std::cerr << "\n[ERROR]: Some of the input information is absent " << std::endl;
	return 1;
  }

  if ((isPID && isDCA) || (!isPID && !isDCA)){
	std::cerr << "\n[ERROR]: Need to choose one QA option (--pid, --dca) " << std::endl;
	return 0;
  }
  basicqa * qa = new basicqa(iFileName, oFileName);

  if (isPID) qa->LoopPID();
  if (isDCA) qa->LoopDCA();
  return 0;
}
