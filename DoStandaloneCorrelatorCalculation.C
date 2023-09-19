// ----------------------------------------------------------------------------
// 'DoStandaloneCorrelatorCalculation.C'
// Derek Anderson
// 02.02.2023
//
// Use this to run the SEnergyCorrelator
// class in standalone mode.
// ----------------------------------------------------------------------------

#ifndef DOSTANDALONECORRELATORCALCULATION_C
#define DOSTANDALONECORRELATORCALCULATION_C

// standard c includes
#include <string>
#include <vector>
#include <cstdlib>
#include <utility>
// user includes
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelator.h"

// load libraries
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libsenergycorrelator.so)

using namespace std;



void DoStandaloneCorrelatorCalculation() {

  // io parameters
  const vector<string> inFile = {
    "../SCorrelatorJetTree/output/condor/final_merge/correlatorJetTree.pa200hijing500bkg010jet20run6_trks.d17m8y2023.root",
    "../SCorrelatorJetTree/output/condor/final_merge/correlatorJetTree.pa200hijing500bkg010jet20run6_trks.d17m8y2023.root"
  };
  const vector<string> inTree = {
    "RecoJetTree",
    "TruthJetTree"
  };
  const vector<string> outFile = {
    "pa200hijing50khzrun6jet20.reco.d22m8y2023.root",
    "pa200hijing50khzrun6jet20.true.d22m8y2023.root"
  };

  // module parameters
  const vector<string> moduleName = {
    "SRecoEnergyCorrelator",
    "STrueEnergyCorrelator"
  };
  const vector<bool> isTruth = {
    false,
    true
  };

  // correlator parameters
  const uint32_t             nPointCorr = 2;
  const uint64_t             nBinsDr    = 75;
  const pair<double, double> binRangeDr = {1e-5, 1.};

  // jet/cst parameters
  const pair<double, double> etaJetRange = {-1., 1.};
  const pair<double, double> momCstRange = {0.,  100.};
  const pair<double, double> drCstRange  = {0.,  5.};

  // jet pT bins
  const vector<pair<double, double>> ptJetBins = {{5., 10.}, {10., 15.}, {15., 20.}, {20., 30.}, {30., 50.}};

  // embedding-specific options
  const int  subEvtOpt = 0;
  const bool isEmbed   = true;

  // misc parameters
  const int  verbosity = 0;
  const bool isComplex = false;
  const bool doDebug   = false;
  const bool inBatch   = false;

  // do correlator calculation on reco jets
  SEnergyCorrelator* recoCorrelator = new SEnergyCorrelator(moduleName[0], isComplex, doDebug, inBatch);
  recoCorrelator -> SetVerbosity(verbosity);
  recoCorrelator -> SetInputFile(inFile[0]);
  recoCorrelator -> SetInputTree(inTree[0], isTruth[0]);
  recoCorrelator -> SetOutputFile(outFile[0]);
  recoCorrelator -> SetJetParameters(ptJetBins, etaJetRange.first, etaJetRange.second);
  recoCorrelator -> SetConstituentParameters(momCstRange.first, momCstRange.second, drCstRange.first, drCstRange.second);
  recoCorrelator -> SetCorrelatorParameters(nPointCorr, nBinsDr, binRangeDr.first, binRangeDr.second);
  recoCorrelator -> Init();
  recoCorrelator -> Analyze();
  recoCorrelator -> End();

  // do correlator calculation on truth jets
  SEnergyCorrelator* trueCorrelator = new SEnergyCorrelator(moduleName[1], isComplex, doDebug, inBatch);
  trueCorrelator -> SetVerbosity(verbosity);
  trueCorrelator -> SetInputFile(inFile[1]);
  trueCorrelator -> SetInputTree(inTree[1], isTruth[1], isEmbed);
  trueCorrelator -> SetOutputFile(outFile[1]);
  trueCorrelator -> SetJetParameters(ptJetBins, etaJetRange.first, etaJetRange.second);
  trueCorrelator -> SetConstituentParameters(momCstRange.first, momCstRange.second, drCstRange.first, drCstRange.second);
  trueCorrelator -> SetCorrelatorParameters(nPointCorr, nBinsDr, binRangeDr.first, binRangeDr.second);
  if (isEmbed) {
    trueCorrelator -> SetSubEventsToUse(subEvtOpt);
  }
  trueCorrelator -> Init();
  trueCorrelator -> Analyze();
  trueCorrelator -> End();
  return;

}

#endif

// end ------------------------------------------------------------------------
