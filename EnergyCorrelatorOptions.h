// ----------------------------------------------------------------------------
// 'EnergyCorrelatorOptions.h'
// Derek Anderson
// 02.08.2024
//
// Options for the SEnergyCorrelator module.
// ----------------------------------------------------------------------------

#ifndef ENERGYCORRELATOROPTIONS_H
#define ENERGYCORRELATOROPTIONS_H 

// c++ utilities
#include <string>
#include <utility>
#include <optional>
// analysis utilities
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelator.h"
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelatorConfig.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/SCorrelatorUtilities.h"

// make common namespacs implicit
using namespace std;
using namespace SColdQcdCorrelatorAnalysis;
using namespace SColdQcdCorrelatorAnalysis::SCorrelatorUtilities;



namespace EnergyCorrelatorOptions {

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
  const int       verbosity = 0;
  const bool      isEmbed   = false;
  const bool      doCstCuts = true;
  const bool      doDebug   = false;
  const SubEvtOpt subEvtOpt = SubEvtOpt::Everything;

  // jet acceptance
  const pair<float, float> etaJetRange = {-0.7, 0.7};

  // cst acceptance
  const pair<float, float> drCstRange  = {0.,  100.};
  const pair<float, float> eneCstRange = {0.1, 100.};



  // bundle acceptances into pairs --------------------------------------------

  pair<JetInfo, JetInfo> GetJetAccept() {

    pair<JetInfo, JetInfo> jetAccept;
    jetAccept.first.eta  = etaJetRange.first;
    jetAccept.second.eta = etaJetRange.second;
    return jetAccept;

  }  // end 'GetJetAccept()'



  pair<CstInfo, CstInfo> GetCstAccept() {

    pair<CstInfo, CstInfo> cstAccept;
    cstAccept.first.dr   = drCstRange.first;
    cstAccept.first.ene  = eneCstRange.first;
    cstAccept.second.dr  = drCstRange.second;
    cstAccept.second.ene = eneCstRange.second;
    return cstAccept;

  }  // end 'GetCstAccept()'



  // set up configurations ----------------------------------------------------

  SEnergyCorrelatorConfig GetRecoConfig(
    const vector<string> cfg_inFiles = inFiles,
    const string cfg_outFile = outFiles.first,
    const int cfg_verbosity = verbosity
  ) {

    SEnergyCorrelatorConfig cfg_reco {
      .verbosity     = cfg_verbosity,
      .isDebugOn     = doDebug,
      .isEmbed       = isEmbed,
      .isInTreeTruth = false,
      .moduleName    = "SRecoEnergyCorrelator",
      .inTreeName    = "RecoJetTree",
      .outFileName   = cfg_outFile,
      .applyCstCuts  = doCstCuts,
      .subEvtOpt     = subEvtOpt,
      .inFileNames   = cfg_inFiles,
      .nPoints       = nPoints,
      .nBinsDr       = nBinsDr,
      .drBinRange    = binRangeDr,
      .ptJetBins     = ptJetBins,
      .jetAccept     = GetJetAccept(),
      .cstAccept     = GetCstAccept()
    };
    return cfg_reco;

  }  // end 'GetRecoConfig(vector<string>, string, int)'



  SEnergyCorrelatorConfig GetTruthConfig(
    const vector<string> cfg_inFiles = inFiles,
    const string cfg_outFile = outFiles.second,
    const int cfg_verbosity = verbosity
  ) {

    SEnergyCorrelatorConfig cfg_true {
      .verbosity     = cfg_verbosity,
      .isDebugOn     = doDebug,
      .isEmbed       = isEmbed,
      .isInTreeTruth = true,
      .moduleName    = "TrueEnergyCorrelator",
      .inTreeName    = "TruthJetTree",
      .outFileName   = cfg_outFile,
      .applyCstCuts  = doCstCuts,
      .subEvtOpt     = subEvtOpt,
      .inFileNames   = cfg_inFiles,
      .nPoints       = nPoints,
      .nBinsDr       = nBinsDr,
      .drBinRange    = binRangeDr,
      .ptJetBins     = ptJetBins,
      .jetAccept     = GetJetAccept(),
      .cstAccept     = GetCstAccept()
    };
    return cfg_true;

  }  // end 'GetConfig(vector<string>, string, int)'

}  // end EnergyCorrelatorOptions namespace

#endif

// end -----------------------------------------------------------------
