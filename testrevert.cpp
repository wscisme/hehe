// C++
#include <iostream>
#include <stdio.h>  
#include <vector>
#include <math.h>
#include <iomanip>

// ROOT
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLegend.h"
#include "TPad.h"
#include "TH2F.h"
#include "TH1.h"
#include "THStack.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"

// CMS2
#include "CMS2.h"
#include "muonSelections.h"
#include "susySelections.h"
#include "eventSelections.h"
#include "MT2/MT2.h"
#include "MT2/MT2Utility.h"

// Versions here!
#include "Version_No.2"


// Good run list
#include "/home/users/jgran/CMSSW_5_3_2_patch4_V05-03-23/src/CMS2/NtupleMacros/Tools/goodrun.cc"

// My includes
#include "myheader.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVec;

using namespace std;
using namespace tas;
using namespace ROOT::Math;

/***************************************************
 * In this file we want to calculate the tt~ cross section and use those events to plot the mt2 and compare the plot with the MC plots

 *
 ***************************************************/


int ScanChain( TChain* chain, char* suffix = "", int nEvtFdr = -1, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // Set Good Run List
  set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    
    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms2.Init(tree);
      
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
    
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      cms2.GetEntry(event);
      ++nEventsTotal;
    
      // Progress	
      CMS2::progress( nEventsTotal, nEventsChain );

      // mt2_bisect::mt2 mt2_event;

      // Select Good Runs

      if(evt_isRealData() && !goodrun( evt_run(), evt_lumiBlock() ) ) {
	++notGoodRun;
	continue;
      }
      else
	++goodRun;

      if(evt_isRealData()){
        DorkyEventIdentifier id = { evt_run(), evt_event(), evt_lumiBlock() };
        if ( is_duplicate(id) ){
          continue;
        }
      }
      isRealData = evt_isRealData(); 
      scale1fb = evt_scale1fb();
      // Analysis Code for MT2

    } //loop over events in the current file

    file_count++;

    // Clean Up
    delete tree;
    file->Close();
    delete file;
  } //file loop
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }


  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << "Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:   " << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:   " << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return 0;
 }

