#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVARegGui.h"
#include "TMVA/Tools.h"

void PlotClassifiers() {

  auto infile = std::unique_ptr<TFile>(TFile::Open("../data/pass6.root"));
  auto _tempTree = std::unique_ptr<TTree>(
      (TTree *)infile->Get("VFTMVABuildTree/TrainingTree"));

  auto reader_b2 = std::unique_ptr<TMVA::Reader>(new TMVA::Reader("V:Color"));
  auto reader_b3 = std::unique_ptr<TMVA::Reader>(new TMVA::Reader("V:Color"));

  int Category;
  float Beta;
  float Rigidity;
  float Mass;
  float RichBDT;

  float TrEdepL2On;
  float TrEdepL2Off;
  float TrEdepInnOff0;
  float TrEdepInnOff1;
  float TrChiSqChY;
  float TrChiSqChX;
  float TrChiSqKaY;
  float TrChiSqKaX;
  float TrHalfRUpCh;
  float TrHalfRDwCh;
  float TrHalfRUpKa;
  float TrHalfRDwKa;
  float L2ChRes_x;
  float L2ChRes_y;
  float L2KaRes_x;
  float L2KaRes_y;
  float L34Scat_x;
  float L34Scat_y;
  float L56Scat_x;
  float L56Scat_y;
  float L56ScatBR_y;
  float L56FeetDist;
  float L2FeetDist;
  float TrQMin;
  float TrQAsymm;

  _tempTree->SetBranchAddress("Category", &Category);
  _tempTree->SetBranchAddress("Mass", &Mass);
  _tempTree->SetBranchAddress("Beta", &Beta);
  _tempTree->SetBranchAddress("Rigidity", &Rigidity);
  _tempTree->SetBranchAddress("RichBDT", &RichBDT);

  _tempTree->SetBranchAddress("TrEdepL2On", &TrEdepL2On);
  _tempTree->SetBranchAddress("TrEdepL2Off", &TrEdepL2Off);
  _tempTree->SetBranchAddress("TrEdepInnOff0", &TrEdepInnOff0);
  _tempTree->SetBranchAddress("TrEdepInnOff1", &TrEdepInnOff1);
  _tempTree->SetBranchAddress("TrChiSqChY", &TrChiSqChY);
  _tempTree->SetBranchAddress("TrChiSqChX", &TrChiSqChX);
  _tempTree->SetBranchAddress("TrChiSqKaY", &TrChiSqKaY);
  _tempTree->SetBranchAddress("TrChiSqKaX", &TrChiSqKaX);
  _tempTree->SetBranchAddress("TrHalfRUpCh", &TrHalfRUpCh);
  _tempTree->SetBranchAddress("TrHalfRDwCh", &TrHalfRDwCh);
  _tempTree->SetBranchAddress("TrHalfRUpKa", &TrHalfRUpKa);
  _tempTree->SetBranchAddress("TrHalfRDwKa", &TrHalfRDwKa);
  _tempTree->SetBranchAddress("L2ChRes_x", &L2ChRes_x);
  _tempTree->SetBranchAddress("L2ChRes_y", &L2ChRes_y);
  _tempTree->SetBranchAddress("L2KaRes_x", &L2KaRes_x);
  _tempTree->SetBranchAddress("L2KaRes_y", &L2KaRes_y);
  _tempTree->SetBranchAddress("L34Scat_x", &L34Scat_x);
  _tempTree->SetBranchAddress("L34Scat_y", &L34Scat_y);
  _tempTree->SetBranchAddress("L56Scat_x", &L56Scat_x);
  _tempTree->SetBranchAddress("L56Scat_y", &L56Scat_y);
  _tempTree->SetBranchAddress("L56FeetDist", &L56FeetDist);
  _tempTree->SetBranchAddress("L2FeetDist", &L2FeetDist);
  _tempTree->SetBranchAddress("TrQMin", &TrQMin);
  _tempTree->SetBranchAddress("TrQAsymm", &TrQAsymm);

  reader_b2->AddSpectator("Category", &Category);
  reader_b2->AddSpectator("Mass", &Mass);
  reader_b2->AddSpectator("Beta", &Beta);
  reader_b2->AddSpectator("Rigidity", &Rigidity);
  reader_b2->AddSpectator("RichBDT", &RichBDT);
  reader_b2->AddVariable("TrEdepL2On", &TrEdepL2On);
  reader_b2->AddVariable("TrEdepL2Off", &TrEdepL2Off);
  reader_b2->AddVariable("TrEdepInnOff0", &TrEdepInnOff0);
  reader_b2->AddVariable("TrEdepInnOff1", &TrEdepInnOff1);
  reader_b2->AddVariable("TrChiSqChY", &TrChiSqChY);
  reader_b2->AddVariable("TrChiSqChX", &TrChiSqChX);
  reader_b2->AddVariable("TrChiSqKaY", &TrChiSqKaY);
  reader_b2->AddVariable("TrChiSqKaX", &TrChiSqKaX);
  reader_b2->AddVariable("TrHalfRUpCh", &TrHalfRUpCh);
  reader_b2->AddVariable("TrHalfRDwCh", &TrHalfRDwCh);
  reader_b2->AddVariable("TrHalfRUpKa", &TrHalfRUpKa);
  reader_b2->AddVariable("TrHalfRDwKa", &TrHalfRDwKa);
  reader_b2->AddVariable("L2ChRes_x", &L2ChRes_x);
  reader_b2->AddVariable("L2ChRes_y", &L2ChRes_y);
  reader_b2->AddVariable("L2KaRes_x", &L2KaRes_x);
  reader_b2->AddVariable("L2KaRes_y", &L2KaRes_y);
  reader_b2->AddVariable("L34Scat_x", &L34Scat_x);
  reader_b2->AddVariable("L34Scat_y", &L34Scat_y);
  reader_b2->AddVariable("L56Scat_x", &L56Scat_x);
  reader_b2->AddVariable("L56FeetDist", &L56FeetDist);
  reader_b2->AddVariable("L2FeetDist", &L2FeetDist);
  reader_b2->AddVariable("TrQMin", &TrQMin);
  reader_b2->AddVariable("TrQAsymm", &TrQAsymm);
  reader_b2->AddVariable("L56ScatBR_y := L56Scat_y*(Beta*Rigidity)",
                         &L56ScatBR_y);

  reader_b3->AddSpectator("Category", &Category);
  reader_b3->AddSpectator("Mass", &Mass);
  reader_b3->AddSpectator("Beta", &Beta);
  reader_b3->AddSpectator("Rigidity", &Rigidity);
  reader_b3->AddSpectator("RichBDT", &RichBDT);
  reader_b3->AddVariable("TrEdepL2On", &TrEdepL2On);
  reader_b3->AddVariable("TrEdepL2Off", &TrEdepL2Off);
  reader_b3->AddVariable("TrEdepInnOff0", &TrEdepInnOff0);
  reader_b3->AddVariable("TrEdepInnOff1", &TrEdepInnOff1);
  reader_b3->AddVariable("TrChiSqChY", &TrChiSqChY);
  reader_b3->AddVariable("TrChiSqChX", &TrChiSqChX);
  reader_b3->AddVariable("TrChiSqKaY", &TrChiSqKaY);
  reader_b3->AddVariable("TrChiSqKaX", &TrChiSqKaX);
  reader_b3->AddVariable("TrHalfRUpCh", &TrHalfRUpCh);
  reader_b3->AddVariable("TrHalfRDwCh", &TrHalfRDwCh);
  reader_b3->AddVariable("TrHalfRUpKa", &TrHalfRUpKa);
  reader_b3->AddVariable("TrHalfRDwKa", &TrHalfRDwKa);
  reader_b3->AddVariable("L2ChRes_x", &L2ChRes_x);
  reader_b3->AddVariable("L2ChRes_y", &L2ChRes_y);
  reader_b3->AddVariable("L2KaRes_x", &L2KaRes_x);
  reader_b3->AddVariable("L2KaRes_y", &L2KaRes_y);
  reader_b3->AddVariable("L34Scat_x", &L34Scat_x);
  reader_b3->AddVariable("L34Scat_y", &L34Scat_y);
  reader_b3->AddVariable("L56Scat_x", &L56Scat_x);
  reader_b3->AddVariable("L56FeetDist", &L56FeetDist);
  reader_b3->AddVariable("L2FeetDist", &L2FeetDist);
  reader_b3->AddVariable("TrQMin", &TrQMin);
  reader_b3->AddVariable("TrQAsymm", &TrQAsymm);
  reader_b3->AddVariable("L56ScatBR_y := L56Scat_y*(Beta*Rigidity)",
                         &L56ScatBR_y);

  reader_b2->BookMVA(
      "BDTDG", "TMVAClass_Agl_2/weights/TMVAClassification_BDTDG.weights.xml");
  reader_b3->BookMVA(
      "BDTDG", "TMVAClass_Agl_3/weights/TMVAClassification_BDTDG.weights.xml");

  TH3D *class_vs_mass =
      new TH3D("class_vs_mass", ";Mass;#Lambda_{2};#Lambda_{3}", 500, -6, 6,
               100, -1, 1, 100, -1, 1);

  float class2, class3;
  int perc=0, prev_perc=0;

  Long64_t nEv = _tempTree->GetEntries();

  // nEv = 1e6;

  for(Long64_t iEv=0; iEv<nEv; iEv++){
  // for (Long64_t iEv = 0; iEv < 100; iEv++) {

    _tempTree->GetEntry(iEv);
    perc = TMath::Floor(100*(float)iEv/nEv);

    // cout << perc << endl;

    if(perc != prev_perc){
      prev_perc = perc;
      cout << perc << "\% - " << iEv << "/" << nEv << "\r";
    }

    if (TMath::Floor(Category / 10) != 3)
      continue;

    L56ScatBR_y = L56Scat_y * Beta * Rigidity;

    class2 = reader_b2->EvaluateMVA("BDTDG");
    class3 = reader_b3->EvaluateMVA("BDTDG");

    class_vs_mass->Fill(Mass, class2, class3);
  }

  auto outfile = std::unique_ptr<TFile>(new TFile("classout.root", "RECREATE"));
  outfile->WriteTObject(class_vs_mass);

}
