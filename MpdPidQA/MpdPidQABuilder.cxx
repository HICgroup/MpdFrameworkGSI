#define MpdPidQABuilder_cxx
#include "MpdPidQABuilder.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <fstream>
#include <stdlib.h>

ClassImp(MpdPidQABuilder);

void MpdPidQABuilder::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   pidQA = new MpdPidQA(sigM, sigE, energy, koef, generator, tracking, nSigPart);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	 if (jentry % 100 == 0) std::cout << "Event N: " << jentry << std::endl;

	 for (Long64_t jtrack=0; jtrack<n_tracks_mpd; jtrack++){

		if (n_hits_mpd[jtrack] < 32) 				continue;
		if (TMath::Abs(signed_pt_mpd[jtrack]) > 3.) continue;
		if (TMath::Abs(eta_mpd[jtrack]) > 1.5)		continue;
		if (TOF_flag_mpd[jtrack] == 0)				continue;
		if (TOF_flag_mpd[jtrack] == 4) 				continue;

		if (signed_pt_mpd[jtrack] <= 0)	charge = 1;
		else							charge = -1;

		pidQA -> FillDedxHists(p_mpd[jtrack],dEdx_tpc_mpd[jtrack],PDG_code_mc[id_from_mc_mpd[jtrack]]);
		pidQA -> Fillm2Hists(p_mpd[jtrack],tof_mass2_mpd[jtrack],PDG_code_mc[id_from_mc_mpd[jtrack]]);
		pidQA -> FillAmplHists(p_mpd[jtrack],PDG_code_mc[id_from_mc_mpd[jtrack]]);
		pidQA -> FillEffContHists(p_mpd[jtrack],dEdx_tpc_mpd[jtrack],tof_mass2_mpd[jtrack],charge,PDG_code_mc[id_from_mc_mpd[jtrack]],0.9);
	 }
   }

	pidQA -> GetDedxQA("TMP/dEdx/");
	pidQA -> Getm2QA("TMP/m2/");
	pidQA -> GetAmplQA("TMP/Ampl/");
	pidQA -> GetEffContQA("TMP/Eff/","pikapr");
}

int main(int argc, char** argv)
{
  TString iFileName, iFileList;
  Bool_t  isFile, isFileList, isMaxDefined;
  Int_t   NmaxLine;

  isFile = kFALSE; isFileList = kFALSE; isMaxDefined = kFALSE;

  if (argc < 3){
	std::cerr << "./MpdPidQABuilder -i inputfile" << std::endl;
    std::cerr << "or " << std::endl;
	std::cerr << "./MpdPidQABuilder -f inputfilelist [OPTIONS]" << std::endl;
	std::cerr << "OPTIONS: " << std::endl;
	std::cerr << "-N [NUMBER] - read [NUMBER] files from the filelist (if -N isn't set, the MpdPidQABuilder will take all files from filelist) " << std::endl;

	return 1;
  }

  for (int i=1; i<argc; i++){
	if(std::string(argv[i]) != "-i" &&
	   std::string(argv[i]) != "-f" &&
	   std::string(argv[i]) != "-N"){
		std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " <<  argv[i] << std::endl;
		return 1;
	} else {
	  if(std::string(argv[i]) == "-i" && i!=argc-1) {
		iFileName = argv[++i];
		isFile = kTRUE;
		continue;
	  }
	  if(std::string(argv[i]) == "-i" && i==argc-1) {
		std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
		return 1;
	  }
	  if(std::string(argv[i]) == "-f" && i!=argc-1) {
		iFileList = argv[++i];
		isFileList = kTRUE;
		continue;
	  }
	  if(std::string(argv[i]) == "-f" && i==argc-1) {
		std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
		return 1;
	  }
	  if(std::string(argv[i]) == "-N" && i!=argc-1) {
		NmaxLine = atoi(argv[++i]);
		isMaxDefined = kTRUE;
		continue;
	  }
	  if(std::string(argv[i]) == "-N" && i==argc-1) {
		std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
		return 1;
	  }
    }
  }

	if (isFile){
		MpdPidQABuilder * QA = new MpdPidQABuilder(iFileName);
		QA -> Loop();
	}
	if (isFileList){
		int Ncount = 0;
		std::ifstream infile(iFileList.Data());
		TChain * chain = new TChain("cbmsim_reduced");
		TString sLine = "";
		while (infile >> sLine)
		{
			if (isMaxDefined && (NmaxLine <= Ncount) ) break;
			chain -> Add(sLine.Data());
			Ncount++;
		}
		TTree * inTree = (TTree*) chain;
		MpdPidQABuilder * QA = new MpdPidQABuilder(inTree);
		QA -> Loop();
	}

	return 0;
}
