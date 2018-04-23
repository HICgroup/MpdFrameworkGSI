#include "reducedTreeCreator.h"

reducedTreeCreator::reducedTreeCreator(TString inFileHistName, TString inFileTreeName, TString outFileName , TString dcaFileName)
{

	pid = new MpdPid(PIDsigM, PIDsigE, PIDenergy, PIDkoeff, PIDgenerator, PIDtracking, PIDparticles);

	inFile = new TFile(inFileTreeName.Data(),"READ");
	inTree = (TTree*) inFile->Get("cbmsim");

	outFile = new TFile(outFileName.Data(),"RECREATE");
	outTree = new TTree("cbmsim_reduced","cbmsim_reduced");
	
	inFileHist = new TFile(inFileHistName.Data(),"READ");
	dcaFile = new TFile(dcaFileName.Data(),"READ");
	
	for (Int_t i_proj=0;i_proj<n_proj;i_proj++){
		//for (Int_t i_pt=0;i_pt<n_pt_bin;i_pt++){
	    		for (Int_t i_eta=0;i_eta<n_eta_bin;i_eta++){
				//f_dca[i_proj][i_pt][i_eta] = (TF1*) dcaFile->Get(Form("dca_fit[%i][%i][%i]",i_proj,i_pt,i_eta));
				f_pt_fit[i_proj][i_eta] = (TF1*) dcaFile->Get(Form("sigma_pt_fit%i%i",i_proj,i_eta));
	    		}
		//}
    	}

	//RNG = new TRandom();
	
	outTree->Branch("b_mc",&b_mc,"b_mc/F");
	outTree->Branch("phiEP_mc",&phiEP_mc,"phiEP_mc/F");
	outTree->Branch("x_vertex_mc",&x_vertex_mc,"x_vertex_mc/F");
	outTree->Branch("y_vertex_mc",&y_vertex_mc,"y_vertex_mc/F");
	outTree->Branch("z_vertex_mc",&z_vertex_mc,"z_vertex_mc/F");
	outTree->Branch("x_vertex_mpd",&x_vertex_mpd,"x_vertex_mpd/F");
	outTree->Branch("y_vertex_mpd",&y_vertex_mpd,"y_vertex_mpd/F");
	outTree->Branch("z_vertex_mpd",&z_vertex_mpd,"z_vertex_mpd/F");
	outTree->Branch("n_tracks_mc",&n_tracks_mc,"n_tracks_mc/L");
	outTree->Branch("eta_mc",eta_mc,"eta_mc[n_tracks_mc]/F");
	outTree->Branch("pt_mc",pt_mc,"pt_mc[n_tracks_mc]/F");
	outTree->Branch("mother_ID_mc",mother_ID_mc,"mother_ID_mc[n_tracks_mc]/I");
	outTree->Branch("PDG_code_mc",PDG_code_mc,"PDG_code_mc[n_tracks_mc]/I");
	outTree->Branch("px_mc",px_mc,"px_mc[n_tracks_mc]");
	outTree->Branch("py_mc",py_mc,"py_mc[n_tracks_mc]");
	outTree->Branch("pz_mc",pz_mc,"pz_mc[n_tracks_mc]");
	outTree->Branch("start_x_mc",start_x_mc,"start_x_mc[n_tracks_mc]");
	outTree->Branch("start_y_mc",start_y_mc,"start_y_mc[n_tracks_mc]");
	outTree->Branch("start_z_mc",start_z_mc,"start_z_mc[n_tracks_mc]");
	outTree->Branch("mass_mc",mass_mc,"mass_mc[n_tracks_mc]");
	outTree->Branch("energy_mc",energy_mc,"energy_mc[n_tracks_mc]");
	outTree->Branch("n_tracks_mpd",&n_tracks_mpd,"n_tracks_mpd/L");
	outTree->Branch("k_tracks_mpd",&k_tracks_mpd,"k_tracks_mpd/L");
	outTree->Branch("eta_mpd",eta_mpd,"eta_mpd[n_tracks_mpd]/F");
	outTree->Branch("phi_mpd",phi_mpd,"phi_mpd[n_tracks_mpd]/F");
	outTree->Branch("theta_mpd",theta_mpd,"theta_mpd[n_tracks_mpd]/F");
	outTree->Branch("TOF_flag_mpd",TOF_flag_mpd,"TOF_flag_mpd[n_tracks_mpd]/I");
	outTree->Branch("ZDC_energy_mpd",ZDC_energy_mpd,"ZDC_energy_mpd[90]/F"); /////////////////////////////////
	outTree->Branch("pid_tpc_prob_electron_mpd",pid_tpc_prob_electron_mpd,"pid_tpc_prob_electron_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tpc_prob_pion_mpd",pid_tpc_prob_pion_mpd,"pid_tpc_prob_pion_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tpc_prob_kaon_mpd",pid_tpc_prob_kaon_mpd,"pid_tpc_prob_kaon_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tpc_prob_proton_mpd",pid_tpc_prob_proton_mpd,"pid_tpc_prob_proton_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_electron_mpd",pid_tof_prob_electron_mpd,"pid_tof_prob_electron_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_pion_mpd",pid_tof_prob_pion_mpd,"pid_tof_prob_pion_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_kaon_mpd",pid_tof_prob_kaon_mpd,"pid_tof_prob_kaon_mpd[n_tracks_mpd]/F");
	outTree->Branch("pid_tof_prob_proton_mpd",pid_tof_prob_proton_mpd,"pid_tof_prob_proton_mpd[n_tracks_mpd]/F");
	outTree->Branch("tof_beta_mpd",tof_beta_mpd,"tof_beta_mpd[n_tracks_mpd]/F");
	outTree->Branch("tof_mass2_mpd",tof_mass2_mpd,"tof_mass2_mpd[n_tracks_mpd]/F");
	outTree->Branch("dEdx_tpc_mpd",dEdx_tpc_mpd,"dEdx_tpc_mpd[n_tracks_mpd]/F");
	outTree->Branch("chi2_mpd",chi2_mpd,"chi2_mpd[n_tracks_mpd]/F");
	outTree->Branch("chi2_vertex",chi2_vertex,"chi2_vertex[n_tracks_mpd]/F");
	outTree->Branch("pt_error_mpd",pt_error_mpd,"pt_error_mpd[n_tracks_mpd]/F");
	outTree->Branch("theta_error_mpd",theta_error_mpd,"theta_error_mpd[n_tracks_mpd]/F");
	outTree->Branch("phi_error_mpd",phi_error_mpd,"phi_error_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_x_mpd",DCA_x_mpd,"DCA_x_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_y_mpd",DCA_y_mpd,"DCA_y_mpd[n_tracks_mpd]/F");
	outTree->Branch("DCA_z_mpd",DCA_z_mpd,"DCA_z_mpd[n_tracks_mpd]/F");
	outTree->Branch("n_hits_mpd",n_hits_mpd,"n_hits_mpd[n_tracks_mpd]/I");
	outTree->Branch("n_hits_poss_mpd",n_hits_poss_mpd,"n_hits_poss_mpd[n_tracks_mpd]/I");
	outTree->Branch("signed_pt_mpd",signed_pt_mpd,"signed_pt_mpd[n_tracks_mpd]/F");
	outTree->Branch("centrality_tpc_mpd",&centrality_tpc_mpd,"centrality_tpc_mpd/I");
	outTree->Branch("id_from_mc_mpd",mpd_side,"id_from_mc_mpd[n_tracks_mpd]/L");
	outTree->Branch("p_mpd",p_mpd,"p_mpd[n_tracks_mpd]/F");
	
	MCHeader = 0;
	MCTracks = 0;
	MPDEvent = 0;
	ZDCHits = 0;
	mpdKalmanTracks = (TClonesArray*) inFile->FindObjectAny("TpcKalmanTrack");
	vertexes = (TClonesArray*) inFile->FindObjectAny("Vertex");
	
	inTree->SetBranchAddress("MCEventHeader.", &MCHeader);
	inTree->SetBranchAddress("MCTrack", &MCTracks);
	inTree->SetBranchAddress("MPDEvent.", &MPDEvent);
	inTree->SetBranchAddress("ZdcDigi",&ZDCHits);
	inTree->SetBranchAddress("TpcKalmanTrack",&mpdKalmanTracks);
	inTree->SetBranchAddress("Vertex",&vertexes);
	
	for (int i = 0; i <= NmultiplicityBins; ++i) multiplicity_bins[i] = 1;
	h_multiplicity_before = (TH1F*)inFileHist->Get("h_multiplicity");
	//h_multiplicity_before = (TH1F*)inFileHist->Get("h_mult_old");
	//Float_t histo_int = h_multiplicity_before->Integral("WIDTH");
	//Float_t one_tenth = histo_int / NmultiplicityBins; //the fraction of events equal to on centrality bin
	Float_t one_tenth = h_multiplicity_before->Integral("WIDTH") / NmultiplicityBins; //the fraction of events equal to on centrality bin
	
	Int_t n_mult_bins = h_multiplicity_before->GetNbinsX(); //total number of bins in multiplicity histogram
	multiplicity_bins[0] = 0.; // the first bin is always zero
	for (int i = 1; i <= NmultiplicityBins; ++i)
	{
		Float_t sum = one_tenth*i;
		multiplicity_bins[i] = integrate(h_multiplicity_before, n_mult_bins, sum);
	}
	//the last bin is always the last
	multiplicity_bins[NmultiplicityBins] = h_multiplicity_before->GetBinLowEdge(n_mult_bins) + h_multiplicity_before->GetBinWidth(1);
	
	for (int i = 0; i <= NmultiplicityBins; ++i) cout << "multiplicity bin = " << multiplicity_bins[i] << endl;
}

void reducedTreeCreator::CreateReducedTree()
{
	Int_t n_pt;
	Int_t n_eta;
	
	Bool_t isGoodPID;
	
	Int_t n_entries = inTree->GetEntries();
    //Int_t n_entries = 50;

    //dcaFile = new TFile("/lustre/nyx/hades/user/parfenov/dca_out_file_TDR.root","READ");

    
    for (int i = 0; i < n_entries; ++i)
    {
		for (long int j = 0; j < _MAX_TRACKS; ++j)
        {
			for (int k = 0; k < 10; ++k)
			{
				mc_side[j][k] = -1;
			}
	    }
	    for (long int j = 0; j < _MAX_TRACKS; ++j)
		{
			mpd_side[j] = -1;
		}
		
		for (long int j = 0; j < 90; ++j)
        {
			ZDC_energy_mpd[j] = 0;
		}
	    
		cout << "EVENT N "<< i <<endl;
		inTree->GetEntry(i);
	    MpdGlobalTracks = MPDEvent->GetGlobalTracks();
	    
	    Int_t number_of_zdchits = ZDCHits->GetEntries();
		for (Int_t zdchit_number = 0; zdchit_number < number_of_zdchits; ++zdchit_number)
		{
			ZDCHit = (MpdZdcDigi*) ZDCHits->At(zdchit_number);
		    Int_t detector_ID = ZDCHit->GetDetectorID();//1,2
	        Int_t module_ID = ZDCHit->GetModuleID();//1-45
			Double_t energy_deposit_per_hit = ZDCHit->GetELoss();
			
			ZDC_energy_mpd[ (detector_ID - 1)*45 + module_ID - 1] += energy_deposit_per_hit;
		}
        
        n_tracks_mpd = MpdGlobalTracks->GetEntriesFast();
        b_mc = MCHeader->GetB();
        //        phiEP_mc = MCHeader->GetEP();
		phiEP_mc = MCHeader->GetRotZ();
	    //phi_event_mc = MCHeader->GetPhi();
	    x_vertex_mc = MCHeader->GetX();
	    y_vertex_mc = MCHeader->GetY();
	    z_vertex_mc = MCHeader->GetZ();
	    
	    MpdVertex *vertex = (MpdVertex*) vertexes->First();
	    TVector3 primaryVertex;
	    vertex->Position(primaryVertex);
	    
	    x_vertex_mpd = primaryVertex.X();
	    y_vertex_mpd = primaryVertex.Y();
	    z_vertex_mpd = primaryVertex.Z();

	    Long_t k_check = 0;

	    for (Int_t track = 0;track<MpdGlobalTracks->GetEntriesFast();track++){
		MpdTrack* mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(track);
		
		n_pt=-1;
		n_eta=-1;
		if (mpdtrack->GetNofHits()<n_hits_cut) continue;

		for (Int_t i_pt=0;i_pt<n_pt_bin;i_pt++)
		    if (TMath::Abs(mpdtrack->GetPt())>pt_bins[i_pt] && TMath::Abs(mpdtrack->GetPt())<=pt_bins[i_pt+1]) n_pt=i_pt;
		for (Int_t i_eta=0;i_eta<n_eta_bin;i_eta++)
		    if (mpdtrack->GetEta()>eta_bins[i_eta] && mpdtrack->GetEta()<=eta_bins[i_eta+1]) n_eta=i_eta;

		if (n_pt==-1) continue;
		if (n_eta==-1) continue;

		TF1 sigma_fit_X = *f_pt_fit[0][n_eta];
		TF1 sigma_fit_Y = *f_pt_fit[1][n_eta];
		TF1 sigma_fit_Z = *f_pt_fit[2][n_eta];
		

		if (TMath::Abs(mpdtrack->GetDCAX()) >= sigma_fit_X(TMath::Abs(mpdtrack->GetPt()))*2) continue;
		if (TMath::Abs(mpdtrack->GetDCAY()) >= sigma_fit_Y(TMath::Abs(mpdtrack->GetPt()))*2) continue;
		if (TMath::Abs(mpdtrack->GetDCAZ()) >= sigma_fit_Z(TMath::Abs(mpdtrack->GetPt()))*2) continue;

		//if (TMath::Abs(mpdtrack->GetEta())>1.5) continue;
		//if (TMath::Abs(mpdtrack->GetPt())<0 && TMath::Abs(mpdtrack->GetPt())>3) continue;
		
		k_check++;
	    }
		k_tracks_mpd = k_check;
        
	    if (k_tracks_mpd>0){
		centrality_tpc_mpd = GetCentrality(k_tracks_mpd);
	    }else centrality_tpc_mpd = -9999; 
	   
	    for (long int j = 0; j < n_tracks_mpd; ++j)
	    { 
			MpdTrack *mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(j);
			MpdKalmanTrack *kalmanTrack = (MpdKalmanTrack*) mpdKalmanTracks->UncheckedAt(j);
			//if (mctrack->GetMotherId() > -1)  continue; ///////////////////////////////////////Mother ID cut is here
			for (int k = 0; k < 10; ++k) 
			{
				if ((mc_side[mpdtrack->GetID()][k]) == -1) { mc_side[mpdtrack->GetID()][k] = j; break;} 
		    }

			int charge_mpd;
			if (mpdtrack->GetPt() < 0) charge_mpd = 1;
			if (mpdtrack->GetPt() >= 0) charge_mpd = -1;

		
			if (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6)
			 	isGoodPID = pid -> FillProbs(TMath::Abs(mpdtrack->GetPt())*TMath::CosH(mpdtrack->GetEta()),
											 mpdtrack->GetdEdXTPC(), mpdtrack->GetTofMass2(), charge_mpd);
			else isGoodPID = pid -> FillProbs(TMath::Abs(mpdtrack->GetPt())*TMath::CosH(mpdtrack->GetEta()),
											 mpdtrack->GetdEdXTPC(), charge_mpd);

			eta_mpd[j] = mpdtrack->GetEta();
			n_hits_mpd[j] = mpdtrack->GetNofHits();
			n_hits_poss_mpd[j] = mpdtrack->GetNofHitsPossTpc();
			pid_tpc_prob_electron_mpd[j] = mpdtrack->GetTPCPidProbElectron();
			if (isGoodPID){
				pid_tpc_prob_pion_mpd[j] = pid -> GetProbPi();//mpdtrack->GetTPCPidProbPion();
				pid_tpc_prob_kaon_mpd[j] = pid -> GetProbKa();//mpdtrack->GetTPCPidProbKaon();
				pid_tpc_prob_proton_mpd[j] = pid -> GetProbPr();//mpdtrack->GetTPCPidProbProton();
			} else {
				pid_tpc_prob_pion_mpd[j] = -999.;
				pid_tpc_prob_kaon_mpd[j] = -999.;
				pid_tpc_prob_proton_mpd[j] = -999.;
			}
			pid_tof_prob_electron_mpd[j] = mpdtrack->GetTOFPidProbElectron();
			pid_tof_prob_pion_mpd[j] = mpdtrack->GetTOFPidProbPion();
			pid_tof_prob_kaon_mpd[j] = mpdtrack->GetTOFPidProbKaon();
			pid_tof_prob_proton_mpd[j] = mpdtrack->GetTOFPidProbProton();
			tof_beta_mpd[j] = mpdtrack->GetTofBeta();
			tof_mass2_mpd[j] = mpdtrack->GetTofMass2();
			dEdx_tpc_mpd[j] = mpdtrack->GetdEdXTPC();
			chi2_mpd[j] = mpdtrack->GetChi2();
			chi2_vertex[j] = kalmanTrack->GetChi2Vertex();
			pt_error_mpd[j] = mpdtrack->GetPtError();
			theta_error_mpd[j] = mpdtrack->GetThetaError();
			phi_error_mpd[j] = mpdtrack->GetPhiError();
			DCA_x_mpd[j] = mpdtrack->GetDCAX();
			DCA_y_mpd[j] = mpdtrack->GetDCAY();
			DCA_z_mpd[j] = mpdtrack->GetDCAZ();
			signed_pt_mpd[j] = mpdtrack->GetPt();
			phi_mpd[j] = mpdtrack->GetPhi();
            theta_mpd[j] = mpdtrack->GetTheta();
            p_mpd[j] = TMath::Abs(signed_pt_mpd[j])/TMath::Sin(theta_mpd[j]);
            TOF_flag_mpd[j] = mpdtrack->GetTofFlag();

	    }

		
		
		n_tracks_mc = MCTracks->GetEntries();
		int m = 0; //ID OF MCTRACK IN THE FINAL TREE ARRAY (AFTER REDUCTION)
		for(long int j = 0; j < n_tracks_mc; ++j)
	    {
			FairMCTrack *mctrack = (FairMCTrack*)MCTracks->At(j);
			if (FillTrack(mctrack,j,m)) continue; //TRACK IS ADDED, GO TO NEXT ONE
			
			if (mctrack->GetMotherId() == -1) // IF MCTRACK HAS NO CORRESPONDING MPDTRACK BUT IT'S PRIMARY, WE WRITE IT ANYWAY
			{	
				
				Float_t theta_mc = TMath::ATan2(mctrack->GetPt(),mctrack->GetPz());
				eta_mc[m] = -TMath::Log(TMath::Tan(theta_mc/2.));
				pt_mc[m] = mctrack->GetPt();
				mother_ID_mc[m] = mctrack->GetMotherId();
				PDG_code_mc[m] = mctrack->GetPdgCode();
				px_mc[m] = mctrack->GetPx();
				py_mc[m] = mctrack->GetPy();
				pz_mc[m] = mctrack->GetPz();
				start_x_mc[m] = mctrack->GetStartX();
				start_y_mc[m] = mctrack->GetStartY();
				start_z_mc[m] = mctrack->GetStartZ();
				mass_mc[m] = mctrack->GetMass();
				energy_mc[m] = mctrack->GetEnergy();  
				
				++m;
			}
			
			//WE DO NOT WRITE SECONDARY TRACKS NOT CORRESPONDING TO ANY MPDTRACKS
		}
		n_tracks_mc = m;
		
		outTree->Fill();
    }
    
    outFile->cd();
    outTree->Write();
    outFile->Close();
}

Float_t reducedTreeCreator::integrate(TH1F *h, Int_t max_bin, Float_t sum)
{
	for (Int_t mult_bin = 1; mult_bin <= max_bin; ++mult_bin)
	{
		if (h->Integral(1,mult_bin,"WIDTH") >= sum)
		{
			//RETURN FLOATING-POINT MULTIPLICITY VALUE
			return h->GetBinLowEdge(mult_bin) + 1 - ( (h->Integral(1,mult_bin,"WIDTH") - sum) / h->GetBinContent(mult_bin));
		}
	}
	return max_bin;
}

bool reducedTreeCreator::FillTrack(FairMCTrack *mctrack, long int j, int &m)
{
	for (int k = 0; k < 10; ++k) //CHECK WHETHER MC TRACK IS ASSOCIATED WITH ANY MPDTRACKS AND TAKE FIRST IF IT IS
	{
		if ( mc_side[j][k] != -1) 
		{
			mpd_side[mc_side[j][k]] = m;
			
			Float_t theta_mc = TMath::ATan2(mctrack->GetPt(),mctrack->GetPz());
			eta_mc[m] = -TMath::Log(TMath::Tan(theta_mc/2.));
			pt_mc[m] = mctrack->GetPt();
			mother_ID_mc[m] = mctrack->GetMotherId();
			PDG_code_mc[m] = mctrack->GetPdgCode();
			px_mc[m] = mctrack->GetPx();
			py_mc[m] = mctrack->GetPy();
			pz_mc[m] = mctrack->GetPz();
			start_x_mc[m] = mctrack->GetStartX();
			start_y_mc[m] = mctrack->GetStartY();
			start_z_mc[m] = mctrack->GetStartZ();
			mass_mc[m] = mctrack->GetMass();
			energy_mc[m] = mctrack->GetEnergy();  
			
			++m;
			return true; //MC TRACK IS ASSOCIATED WITH ONE OR MORE MPDTRACKS -> TRACK IS WRITTEN INTO FINAL TREE -> TRUE 
		}
	}
	
	return false; //MC TRACK IS NOT ASSOCIATED WITH ANY MPDTRACKS -> TRACK IS NOT WRITTEN INTO FINAL TREE -> FALSE
}

Int_t reducedTreeCreator::GetCentrality(Int_t multiplicity)
{
	int centrality_bin = -1;
	for (int multiplicityBin = 0; multiplicityBin < NmultiplicityBins; ++multiplicityBin)
	{
		if ((multiplicity > (Int_t)multiplicity_bins[multiplicityBin]) && (multiplicity <= (Int_t)multiplicity_bins[multiplicityBin+1])) //inclusive borders
		{
			if ((multiplicity == (Int_t)multiplicity_bins[multiplicityBin+1]))
			{
				if ((multiplicity == (Int_t)multiplicity_bins[multiplicityBin+2]))
				{
					
					Double_t random_number = RNG->Uniform();
					if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+1]) &&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+2]))
						return multiplicityBin+1;
					//else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+2]))
					else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+2]) &&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+3]))
						return multiplicityBin+2;
					else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+3]) &&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+4]))
						return multiplicityBin+3;
					else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+4]) &&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+5]))
						return multiplicityBin+4;
					else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+5]) &&((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) < multiplicity_bins[multiplicityBin+6]))
						return multiplicityBin+5;
					else if (((random_number +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+6]))
						return multiplicityBin + 6;
					else 
						return multiplicityBin;
				}
				if ((RNG->Uniform() +(Int_t)multiplicity_bins[multiplicityBin+1]) > multiplicity_bins[multiplicityBin+1]) 
					return multiplicityBin+1; 
				else 
					return multiplicityBin;
				//~ Int_t n = 0;
				//~ while (multiplicity == (Int_t)multiplicity_bins[multiplicityBin + 1 + n]) ++n;
				//~ cout << "multiplicity = "<< multiplicity <<"; n = " << n << endl;
				//~ Float_t *array = new Float_t[n+2];
				//~ array[0] = (Int_t)multiplicity_bins[multiplicityBin + 1];
				//~ for (int k = 1; k < n+2; ++k) array[k] = multiplicity_bins[k];
				//~ Double_t random_number = RNG->Uniform();
				//~ for (int j = 0; j < n + 1; ++j) if ( (((Int_t)multiplicity_bins[multiplicityBin+1] + random_number) > array[j]) &&
					//~ (((Int_t)multiplicity_bins[multiplicityBin+1] + random_number) < array[j+1]) ) return  multiplicityBin + j;
				//~ delete[] array;
			}
			return multiplicityBin;
		}
	}
	return centrality_bin;
}

int main(int argc, char** argv)
{
  TString iFileName, oFileName, dcaFileName, centFileName;
  if (argc < 3){
	std::cerr << "./reducedTreeCreator -i inputfile -o outputfile -dca dcafile -centrality centralityfile" << std::endl;
	return 1;
  }

  for (int i=1; i<argc; i++){
	if(std::string(argv[i]) != "-i" &&
	  std::string(argv[i]) != "-dca" &&
	  std::string(argv[i]) != "-centrality" &&
	  std::string(argv[i]) != "-o"){
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
	  if(std::string(argv[i]) == "-dca" && i!=argc-1) {
		dcaFileName = argv[++i];
		continue;
	  }
	  if(std::string(argv[i]) == "-dca" && i==argc-1) {
		std::cerr << "\n[ERROR]: DCA file name was not specified " << std::endl;
		return 1;
	  }
	  if(std::string(argv[i]) == "-centrality" && i!=argc-1) {
		centFileName = argv[++i];
		continue;
	  }
	  if(std::string(argv[i]) == "-centrality" && i==argc-1) {
		std::cerr << "\n[ERROR]: Centrality file name was not specified " << std::endl;
		return 1;
	  }	   
    }
  }
  reducedTreeCreator* rtc = new reducedTreeCreator(centFileName, iFileName, oFileName, dcaFileName);
  rtc->CreateReducedTree();
  return 0;
}

ClassImp(reducedTreeCreator);
