#ifndef RESTORE_DCA_H
#define RESTORE_DCA_H

//Seting up undefind value. To check if something was written here one can make if(value!=UNDEFINED_DCA)
#define UNDEFINED_DCA -9999

#include "MpdFillDstTask.h"
#include "MpdKalmanTrack.h"
#include "MpdEctKalmanTrack.h"
#include "MpdTofMatching.h"
#include "MpdEtofMatching.h"
#include "MpdVertex.h"
#include "MpdHelix.h"
#include "MpdParticleIdentification.h"

#include "FairMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>

#include <iostream>
#include <TString.h>
using namespace std;

const TString DCA_INPUT_TREE_NAME = "cbmsim";
const TString MPD_KALMAN_TRACKS_BRANCH_NAME = "TpcKalmanTrack";
const TString MC_TRACK_BRANCH_NAME = "MCTrack";
const TString MPD_EVENT_TRACKS_BRANCH_NAME = "MPDEvent.";
const TString VERTEXES_BRANCH_NAME = "Vertex";

class MpdGlobalTracks;
class MpdHelix;
class MpdKalmanTrack;
class TTree;
class TClonesArray;

//Function that returns correct 3D helicity in MpdHelix format
MpdHelix MakeHelix(const MpdKalmanTrack *tr);

//Filling values=UNDEFINED_DCA
void fill_zeroes(TClonesArray *MpdGlobalTracks);

//Main function
//
//USAGE:
//root -l rootlogon.C
//  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C")
//  mpdloadlibs(kTRUE,kTRUE)
//  .L restore_dca.c+
//  restore_dca("inFileName","outFileName")
//
//ARGUMENTS:
// inFileName - input cbmsim tree with uncorrect dca values
// outFileName - output cbmsim tree with correct dca values. Everything else is the same as it was in inFileName
//
void restore_dca(const TString &inFileName, const TString &outFileName);

void restore_dca(TTree *inTree, TTree *outTree,
                 TClonesArray *mpdKalmanTracks, TClonesArray *MCTracks,
                 MpdEvent *MPDEvent, TClonesArray *vertexes);
                 
void restore_dca(const TString &inFileName, const TString &inFileTreeName, const TString &outFileName, const TString &mpdKalmanTracksBranchName,
                 const TString &MCTracksBranchName, const TString &MPDEventName, const TString &vertexesBranchName);

#endif // RESTORE_DCA_H
