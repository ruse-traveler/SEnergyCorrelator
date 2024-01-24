// ----------------------------------------------------------------------------
// 'DoStandaloneCorrelatorCalculation.cxx'
// Derek Anderson
// 02.02.2023
//
// Use this to run the SEnergyCorrelator
// class in standalone mode.
//
// NOTE: subEvtOpt controls which sub-event
// to use in MC
//   0 -- use everything
//   1 -- use only signal subevent
//   2 -- use only background subevents
//        (pileup and underlying event)
//   3 -- use only primary background
//        (i.e. underlying event)
//   4 -- use only pileup subevents
//   5 -- use specific set of subevents
// ----------------------------------------------------------------------------

#ifndef DOSTANDALONECORRELATORCALCULATION_C
#define DOSTANDALONECORRELATORCALCULATION_C

// c++ utilities
#include <string>
#include <vector>
#include <cstdlib>
#include <utility>
// analysis utilities
#include "/sphenix/user/danderson/eec/SCorrelatorUtilities/CstTools.h"
#include "/sphenix/user/danderson/eec/SCorrelatorUtilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelator.h"
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelatorConfig.h"

using namespace std;
using namespace SColdQcdCorrelatorAnalysis;
using namespace SColdQcdCorrelatorAnalysis::ScorrelatorUtilities;

// load libraries
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libsenergycorrelator.so)

// default options
static const bool DebugDefault = false;



// macro body -----------------------------------------------------------------

void DoStandaloneCorrelatorCalculation() {

  // options ------------------------------------------------------------------

  // input/output files
  const vector<string> inFiles = {
    "../SCorrelatorJetTree/output/condor/final_merge/correlatorJetTree.pp200py8jet10run8_trksWithRoughCuts.d26m9y2023.root"
  };
  const pair<string, string> outFiles = {
    "twoPoint.pp200py8jet10run8.refactor_afterConfigStruct_reco.d18m1y2024.root",
    "twoPoint.pp200py8jet10run8.refactor_afterConfigStruct_true.d18m1y2024.root"
  };

  // correlator parameters
  const vector<uint32_t>     nPoints    = {2};
  const pair<double, double> binRangeDr = {1e-5, 1.};
  const uint64_t             nBinsDr    = 75;

  // jet pT bins
  const vector<pair<double, double>> ptJetBins = {
    {5., 10.},
    {10., 15.},
    {15., 20.},
    {20., 30.},
    {30., 50.}
  };

  // misc options
  const int  subEvtOpt = 0;
  const bool isEmbed   = false;
  const bool doCstCuts = true;

  // jet acceptance
  pair<JetInfo, JetInfo> cfg_jetAccept;
  cfg_jetAccept.first.eta  = -0.7;
  cfg_jetAccept.second.eta = 0.7;

  // cst acceptance
  pair<CstInfo, CstInfo> cfg_cstAccept;
  cfg_cstAccept.first.dr  = 0.;
  cfg_cstAccept.first.pt  = 0.1;
  cfg_cstAccept.second.dr = 100.;
  cfg_cstAccept.second.pt = 100.;

  // reco configuration
  SEnergyCorrelatorConfig cfg_reco {
    .isDebugOn    = DebugDefault,
    .isEmbed      = isEmbed,
    .sInTreeTruth = false,
    .moduleName   = "SRecoEnergyCorrelator",
    .inTreeName   = "RecoJetTree",
    .outFileName  = outFiles.first,
    .applyCstCuts = doCstCuts,
    .subEvtOpt    = subEvtOpt,
    .inFileNames  = inFiles,
    .nPoints      = nPoints,
    .nBinsDr      = nBinsDr,
    .drBinRange   = binRangeDr,
    .ptJetBins    = ptJetBins,
    .inFileNames  = inFiles,
    .jetAccept    = cfg_jetAccept,
    .cstAccept    = cfg_cstAccept
  };

  // truth configuration
  SEnergyCorrelatorConfig cfg_true {
    .isDebugOn    = DebugDefault,
    .isEmbed      = isEmbed,
    .sInTreeTruth = true,
    .moduleName   = "TrueEnergyCorrelator",
    .inTreeName   = "TruthJetTree",
    .outFileName  = outFiles.second,
    .applyCstCuts = doCstCuts,
    .subEvtOpt    = subEvtOpt,
    .inFileNames  = inFiles,
    .nPoints      = nPoints,
    .nBinsDr      = nBinsDr,
    .drBinRange   = binRangeDr,
    .ptJetBins    = ptJetBins,
    .jetAccept    = cfg_jetAccept,
    .cstAccept    = cfg_cstAccept
  };

  // run calculations ---------------------------------------------------------

  // do correlator calculation on reco jets
  SEnergyCorrelator* recoCorrelator = new SEnergyCorrelator(cfg_reco);
  recoCorrelator -> Init();
  recoCorrelator -> Analyze();
  recoCorrelator -> End();

  // do correlator calculation on truth jets
  SEnergyCorrelator* trueCorrelator = new SEnergyCorrelator(cfg_true);
  trueCorrelator -> Init();
  trueCorrelator -> Analyze();
  trueCorrelator -> End();
  return;

}

#endif

// end ------------------------------------------------------------------------
