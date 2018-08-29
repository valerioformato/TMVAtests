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

#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/TMVARegGui.h"
#include "TMVA/Tools.h"

#include "parser.h"

using namespace TMVA;

static std::vector<TString> trainingVars = {
    "TrEdepL2On",  "TrEdepL2Off", "TrEdepInnOff0", "TrEdepInnOff1",
    "TrChiSqChY",  "TrChiSqChX",  "TrChiSqKaY",    "TrChiSqKaX",
    "TrHalfRUpCh", "TrHalfRDwCh", "TrHalfRUpKa",   "TrHalfRDwKa",
    "L2ChRes_x",   "L2ChRes_y",   "L2KaRes_x",     "L2KaRes_y",
    "L34Scat_x",   "L34Scat_y",   "L56Scat_x",     "L56FeetDist",
    "L2FeetDist",  "TrQMin",      "TrQAsymm",
};
static std::vector<TString> spectatorVars = {"Category", "Mass", "Beta",
                                             "Rigidity", "RichBDT"};

int main(int argc, char **argv) {

  struct arguments arguments;
  arguments.flag_overwrite = false;
  arguments.flag_debug = false;
  arguments.flag_input = false;
  arguments.flag_output = false;

  int parse_status = argp_parse(&argp, argc, argv, 0, 0, &arguments);

  if (parse_status)
    throw 42;

  if (!arguments.flag_input) {
    std::cerr << "ERROR: Missing input file" << std::endl;
    return 1;
  }
  TString infilename = arguments.infilename;

  if (!arguments.flag_output) {
    std::cerr << "ERROR: Missing output file" << std::endl;
    return 1;
  }
  TString outfilename = arguments.outfilename;

  int nTrees = arguments.nTrees;
  int maxDepth = arguments.maxDepth;
  int nCuts = arguments.nCuts;

  int bkgType = 2;
  int detCat = 3;

  TMVA::Tools::Instance();
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1e4;

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  TFile *outputFile = TFile::Open(outfilename + ".root", "RECREATE");

  auto factory = std::unique_ptr<TMVA::Factory>(
      new TMVA::Factory("TMVAClassification", outputFile,
                        "V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;"
                        "G,D;G_Signal,D_Signal:AnalysisType=Classification"));

  auto dataloader =
      std::unique_ptr<TMVA::DataLoader>(new TMVA::DataLoader(outfilename));

  for (auto var : trainingVars) {
    dataloader->AddVariable(var, var, "", 'F');
  }
  dataloader->AddVariable("L56ScatBR_y := L56Scat_y*(Beta*Rigidity)",
                          "L56ScatBR_y", "", 'F');

  TFile *infile = TFile::Open(infilename);
  TTree *regTree = (TTree *)infile->Get("TrainingTree");
  Double_t signalWeight = 1.0;
  Double_t backgroundWeight = 1.0;
  dataloader->AddSignalTree(regTree, signalWeight);
  dataloader->AddBackgroundTree(regTree, backgroundWeight);

  TString bkgString = (TString) "(Category%10)==";
  bkgString += bkgType;

  TCut detCatCut = Form("TMath::Floor(Category/10)==%i", detCat);
  TCut mycuts = detCatCut + "(Category%10)==1";
  TCut mycutb = detCatCut + (TCut)bkgString;
  TString evSepString = "";
  dataloader->PrepareTrainingAndTestTree(
      mycuts, mycutb,
      "SplitMode=Random:nTrain_Signal=150000:nTest_Signal=200000:NormMode="
      "NumEvents:V");

  TString parString =
      Form("!H:V:NTrees=%i:MinNodeSize=5%%:MaxDepth=%i:BoostType=AdaBoost:"
           "SeparationType=GiniIndex:nCuts=%i:VarTransform=G,D",
           nTrees, maxDepth, nCuts);
  factory->BookMethod(dataloader.get(), TMVA::Types::kBDT, "BDTDG", parString);

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  return 0;
}
