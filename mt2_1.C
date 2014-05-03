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

  const int FIRTBIN = 0;
  const int LASTBIN = 150;
  const int BINNUM = 60;
  const int INVMLASTBIN = 300;
  //const int METLASTBIN = 150;

////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *MT2_ttbar = new TH1F("MT2","MT2", BINNUM, FIRTBIN, LASTBIN+1);
  TH1F *h2_nozero = new TH1F("hnz","MT2 without 0 bin", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_emuCut50  = new TH1F("hem50","MT2 without 0 bin for emu with met cut 50", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_emu  = new TH1F("hem","MT2 without 0 bin for emu  events", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_ee   = new TH1F("hee","MT2 without 0 bin for ee   events", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_mumu = new TH1F("hmm","MT2 without 0 bin for mumu events", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_emj2 = new TH1F("hemj2","MT2 ex 0 for emu  events with 2 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_eej2 = new TH1F("heej2","MT2 ex 0 for ee   events with 2 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_mmj2 = new TH1F("hmmj2","MT2 ex 0 for mumu events with 2 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_emj3 = new TH1F("hemj3","MT2 ex 0 for emu  events with 3 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_eej3 = new TH1F("heej3","MT2 ex 0 for ee   events with 3 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_mmj3 = new TH1F("hmmj3","MT2 ex 0 for mumu events with 3 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_emj4 = new TH1F("hemj4","MT2 ex 0 for emu  events with 4 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_eej4 = new TH1F("heej4","MT2 ex 0 for ee   events with 4 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_mmj4 = new TH1F("hmmj4","MT2 ex 0 for mumu events with 4 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_emj5 = new TH1F("hemj5","MT2 ex 0 for emu  events with >= 5 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_eej5 = new TH1F("heej5","MT2 ex 0 for ee   events with >= 5 jets", BINNUM, 1, LASTBIN+1);
  TH1F *MT2_mmj5 = new TH1F("hmmj5","MT2 ex 0 for mumu events with >= 5 jets", BINNUM, 1, LASTBIN+1);


  TH1F *InvM_lep  = new TH1F("h5","Dileptons InvMass all",	     70, 0, INVMLASTBIN+1);
  TH1F *InvM_ee   = new TH1F("h6","Dilepton InvMass of ee events",   70, 0, INVMLASTBIN+1);
  TH1F *InvM_emu  = new TH1F("h7","Dilepton InvMass of emu events",  70, 0, INVMLASTBIN+1);
  TH1F *InvM_mumu = new TH1F("h8","Dilepton InvMass of mumu events", 70, 0, INVMLASTBIN+1);
  TH1F *JetMult_a = new TH1F("h9","Jet mutiplicity all", 7, 0, 7);
  TH1F *JetMult_b = new TH1F("h10","b-Jet mutiplicity", 7, 0, 7);

///////////////////////////////////////////////////////////////////////////////////////////

  int file_count = 0;

  int notGoodRun = 0;
  int goodRun = 0;
  int reject_at_selection = 0;
  int no_good_pair = 0;
  int less_jets = 0;
  int nttbarEvents = 0;
  int metCut30 = 0;
  int metCut50 = 0;
  int emuZeroBin = 0;
  int noGoodVtx = 0;
  bool isRealData = true;  
  float nttbarScaled = 0;

  float bTagDiscriminant = 0.244; 
  float scale1fb = -1;

  // Loop over events to Analyze
  // Int_t allFileEntriesTotal = 27137253;     // for mcdy
  // Int_t allFileEntriesTotal = 11794288;     // for mctt 
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
      int n_good = -1;	
      int n_nextGood = -1;	
      Int_t n_jets = 0;
      float pt_max = 0;
      float pt_nextMax = 0;
      float lumi = 20.0;
      Int_t n_bTag = 0;
      float M_ll = -9;


      // Analysis Code for MT2
      for(unsigned int i=0; i<hyp_p4().size(); i++){
       // hyp: hypothetical 2 lep events;
       // ll: lep1, lt: lep2; describe the criteria strength 

        // int hyp_type_good = hyp_type().at(i);  
        int lt_id  = cms2.hyp_lt_id().at(i);
        int ll_id  = cms2.hyp_ll_id().at(i);
        int lt_idx = cms2.hyp_lt_index().at(i); 
        int ll_idx = cms2.hyp_ll_index().at(i);
	
	// select only e mu os event to eliminate containmination from Z
	// if(lt_id*ll_id != -143) continue;
	// throw the events near Z channel

       // select opposite charged dilepton
	if(hyp_ll_charge().at(i)*hyp_lt_charge().at(i) > 0)  continue;
       // pt cut
        if(hyp_ll_p4().at(i).pt()  < 20.0)   continue;
        if(hyp_lt_p4().at(i).pt()  < 20.0)   continue;        

       // eta cut
        if(abs(ll_id) == 11 && fabs(hyp_ll_p4().at(i).eta()) > 2.5)   continue; // but 2.4 for the baby
        if(abs(ll_id) == 13 && fabs(hyp_ll_p4().at(i).eta()) > 2.1)   continue;
        if(abs(lt_id) == 11 && fabs(hyp_lt_p4().at(i).eta()) > 2.5)   continue;
        if(abs(lt_id) == 13 && fabs(hyp_lt_p4().at(i).eta()) > 2.1)   continue;

	// To avoid no good vertex warning
	if ( firstGoodVertex() < 0) { noGoodVtx++; continue;}

       // confirm that they pass tests
        if( (abs(ll_id) == 11) && (!passElectronSelection_ZMet2012_v3_Iso(ll_idx)) )  continue; 
        if( (abs(lt_id) == 11) && (!passElectronSelection_ZMet2012_v3_Iso(lt_idx)) )  continue;
        if( (abs(ll_id) == 13) && (!muonId(ll_idx, ZMet2012_v1)) ) continue;
        if( (abs(lt_id) == 13) && (!muonId(lt_idx, ZMet2012_v1)) ) continue;

	float invM_ll = (hyp_ll_p4().at(i) + hyp_lt_p4().at(i)).M();
       // depress D-Y process by requiring InvMas cut
	if(invM_ll < 60)						continue; 
       // Z mass veto
	if(lt_id*ll_id != -143 && (fabs(invM_ll - 91.2)<15) )		continue;

       // find the greatest scalar sum of dileps at position i
	float pt_sum = hyp_ll_p4().at(i).pt() + hyp_lt_p4().at(i).pt();
	if (pt_sum > pt_max) {
	  pt_max = pt_sum;
	  n_good = i;
	  M_ll = invM_ll;
	}
      }// loop hyp_p4().size() 
       // hyp:all possible combination of lep pairs, choose the pair has highest scalar sum

      if(n_good == -1) { no_good_pair++; continue;}
      // for later usage
      int ltid = hyp_lt_id().at(n_good);
      int llid = hyp_ll_id().at(n_good);

      //// tt~ event selection ////
      if(pfjets_p4().size() < 2) {less_jets++; continue;}

      for(unsigned int c = 0; c < pfjets_p4().size(); c++) {
      	float _jetPt = pfjets_p4().at(c).pt() * pfjets_corL1FastL2L3().at(c);
      	// jet pt times correction to jets
      	float dr_lt = VectorUtil::DeltaR(pfjets_p4().at(c), hyp_lt_p4().at(n_good));
      	float dr_ll = VectorUtil::DeltaR(pfjets_p4().at(c), hyp_ll_p4().at(n_good));
      	// delta R is the distance in the eta-phi plane
      	if(_jetPt < 30) continue;       // disgard those with small pt

      	if(dr_lt < 0.4) continue;
      	if(dr_ll < 0.4) continue;      // disgard small delta R 

	if(fabs(pfjets_p4().at(c).eta()) > 2.5) continue;
	if(fabs(pfjets_p4().at(c).eta()) > 2.5) continue; // May here cause less tt~ events as dumping higher eta jets

	float _bTag = pfjets_combinedSecondaryVertexBJetTag().at(c);
	if ( _bTag > bTagDiscriminant)  n_bTag++;	 // mark as the has qualified b-jet, _bTag should be 0~1

      	// float _jetdp = abs(pfjets_p4().at(c).pt() - _jetPt);  
      	// should be diff in correction, tested equal to pt*abs(1-correc)
      	n_jets++;					 // as qualified jet multiplicity

      }//pfjets_p4().size()

      if(n_jets < 2 || n_bTag == 0) {less_jets++; continue;}

      float metcor = evt_pfmet_type1cor();
     // MET correction = evt_particleFloat_type1correction

      // MET not cut!
      if(llid*ltid != -143 && metcor < 30) {metCut30++; }//continue;}

      nttbarEvents++;
     // All selection finished
      nttbarScaled += evt_scale1fb();

      if(metcor < 50) {metCut50++; }//continue;}
/*
      if(llid*ltid == -121) {
	if(metcor > METLASTBIN)  MetDis_ee->Fill((float) METLASTBIN);
	else MetDis_ee->Fill(metcor);
      }
      else if(llid*ltid == -143) {
	if(metcor > METLASTBIN)  MetDis_em->Fill((float) METLASTBIN);
	else MetDis_em->Fill(metcor);
      }
      else if(llid*ltid == -169) {
	if(metcor > METLASTBIN)  MetDis_mm->Fill((float) METLASTBIN);
	else MetDis_mm->Fill(metcor);
      }
*/

      // Fill the jet mutiplicity to the histogram
      if(n_jets > 7)  JetMult_a->Fill((float) 7);
      else JetMult_a->Fill(n_jets);
      if(n_bTag > 7)  JetMult_b->Fill((float) 7);
      else JetMult_b->Fill(n_bTag);

      float mt2 = MT2(metcor, evt_pfmetPhi_type1cor(), hyp_ll_p4().at(n_good), hyp_lt_p4().at(n_good));

      if(ltid*llid == -121) {
	InvM_ee->Fill(M_ll);

	MT2_ee->Fill(mt2);
	if(n_jets == 2)        MT2_eej2->Fill(mt2);
	else if(n_jets == 3)   MT2_eej3->Fill(mt2);
	else if(n_jets == 4)   MT2_eej4->Fill(mt2);
	else if(n_jets > 4)    MT2_eej5->Fill(mt2);
	else cout << "Bug for jet multi !!\n";
      }
      else if(ltid*llid == -143) {
	if(M_ll > (float)INVMLASTBIN)  InvM_emu->Fill((float)INVMLASTBIN);
	else InvM_emu->Fill(M_ll);

	MT2_emu->Fill(mt2);
	if(n_jets == 2)        MT2_emj2->Fill(mt2);
	else if(n_jets == 3)   MT2_emj3->Fill(mt2);
	else if(n_jets == 4)   MT2_emj4->Fill(mt2);
	else if(n_jets > 4)    MT2_emj5->Fill(mt2);
	else cout << "Bug for jet multi !!\n";

	if(evt_pfmet() > 50) MT2_emuCut50->Fill(mt2);
      }
      else if(ltid*llid == -169) {
	InvM_mumu->Fill(M_ll);

	MT2_mumu->Fill(mt2);
	if(n_jets == 2)        MT2_mmj2->Fill(mt2);
	else if(n_jets == 3)   MT2_mmj3->Fill(mt2);
	else if(n_jets == 4)   MT2_mmj4->Fill(mt2);
	else if(n_jets > 4)    MT2_mmj5->Fill(mt2);
	else cout << "Bug for jet multi !!\n";
      }
      else cout << "Uncounted event: " << ltid << ' '<<llid <<endl;

      if(M_ll > (float)INVMLASTBIN)  InvM_lep->Fill((float)INVMLASTBIN);
      else InvM_lep->Fill(M_ll);

      if(mt2 > (float)LASTBIN){
	MT2_ttbar->Fill((float)LASTBIN);
	h2_nozero->Fill((float)LASTBIN);
      }
      else{
	MT2_ttbar->Fill(mt2);
	h2_nozero->Fill(mt2);
      }

      // GenMet->Fill(gen_met());
      // RecMet->Fill(metcor);

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

  cout << "\nNumber of total Events: " << nEventsTotal << "  identified tt~ events: " << nttbarEvents 
       //<< " at 20GeV cut, " << nttbarEvents-metCut30 
       << " raw 0 bin for emu is: "<< emuZeroBin
       <<endl <<endl;
  cout << "For the data samples: "
       << "\n# notGoodRun: " << notGoodRun 
       << "\n# goodRun: " << goodRun
       << "\n# noGoodVtx: " << noGoodVtx
       << endl;

  // scaling the hists to the whole folder for evt_scale1fb to take effect
  // Double_t scale = nEvtFdr/nEventsChain;
  // if(! isRealData)    scale *= scale1fb;

  // MT2_emu->Scale(scale);
  // MT2_emuCut50->Scale(scale);
  // MT2_ttbar->Scale(scale);
  // h2_nozero->Scale(scale);

  // InvM_ee->Scale(scale);
  // InvM_emu->Scale(scale);
  // InvM_mumu->Scale(scale);
  // InvM_lep->Scale(scale);
  // JetMult_a->Scale(scale);
  // JetMult_b->Scale(scale);

  // MT2_ee  ->Scale(scale);
  // MT2_mumu->Scale(scale);
  // MT2_emj2->Scale(scale);
  // MT2_eej2->Scale(scale);
  // MT2_mmj2->Scale(scale);
  // MT2_emj3->Scale(scale);
  // MT2_eej3->Scale(scale);
  // MT2_mmj3->Scale(scale);
  // MT2_emj4->Scale(scale);
  // MT2_eej4->Scale(scale);
  // MT2_mmj4->Scale(scale);
  // MT2_emj5->Scale(scale);
  // MT2_eej5->Scale(scale);
  // MT2_mmj5->Scale(scale);
  
  // nttbarScaled *= scale; // no need for running on baby

  // char* suffix = "_1_mcdy";
  TFile* fout = new TFile(Form("./hists/mt2_hists%s.root",suffix),"RECREATE");
  MT2_emuCut50->Write();
  MT2_emu->Write();
  MT2_ttbar->Write();
  h2_nozero->Write();
  InvM_ee->Write();
  InvM_emu->Write();
  InvM_mumu->Write();
  InvM_lep->Write();
  JetMult_a->Write();
  JetMult_b->Write();

  MT2_ee  ->Write();
  MT2_mumu->Write();
  MT2_emj2->Write();
  MT2_eej2->Write();
  MT2_mmj2->Write();
  MT2_emj3->Write();
  MT2_eej3->Write();
  MT2_mmj3->Write();
  MT2_emj4->Write();
  MT2_eej4->Write();
  MT2_mmj4->Write();
  MT2_emj5->Write();
  MT2_eej5->Write();
  MT2_mmj5->Write();  

  // MetDis_ee->Write();
  // MetDis_em->Write();
  // MetDis_mm->Write();
  fout->Close();

  TCanvas* c1 = new TCanvas;
  // TLegend* l1 = new TLegend(0.4,0.1,1,0.4,"Legend");
  // l1->AddEntry(MT2_hist,"MT2","f");

  // MT2_ttbar->Draw();
  // // c1->BuildLegend(0.2, 0.1, 0.3, 0.2);
  // c1->SaveAs(Form("./hists/mt2%s.png", suffix));
  // h2_nozero->Draw();
  // c1->SaveAs(Form("./hists/mt2_nonzero%s.png", suffix));

  //TCanvas* c2 = new TCanvas;		   
					   
  TCanvas* c3 = new TCanvas;		   
  MT2_emu->Draw();			   
  MT2_emuCut50->Draw("same");
  MT2_emuCut50->SetLineColor(kRed+1);
  // c3->SaveAs(Form("./hists/mt2_emuonly%s.png", suffix));
					   
  //  c1->SetLogy();			   
  // c1->Close();			   
  					   

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

