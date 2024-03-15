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
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Types.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Constants.h"
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelatorConfig.h"

// make common namespacs implicit
using namespace std;
using namespace SColdQcdCorrelatorAnalysis;



namespace EnergyCorrelatorOptions {

  // options ------------------------------------------------------------------

  // input/output files
  const vector<string> inFiles = {
    "../SCorrelatorJetTree/output/condor/final_merge/correlatorJetTree.pp200py8jet10run8_trksWithRoughCuts.d26m9y2023.root"
  };
  const pair<string, string> outFiles = {
    "twoPoint.pp200py8jet10run8.refactor_afterUtilityRearchitect_reco.d15m3y2024.root",
    "twoPoint.pp200py8jet10run8.refactor_afterUtilityRearchitect_true.d15m3y2024.root"
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
  const int  verbosity = 0;
  const int  subEvtOpt = Const::SubEvtOpt::Everything;
  const bool isEmbed   = false;
  const bool doCstCuts = true;
  const bool doDebug   = false;

  // jet acceptance
  const pair<float, float> etaJetRange = {-0.7, 0.7};

  // cst acceptance
  const pair<float, float> drCstRange  = {0.,  100.};
  const pair<float, float> eneCstRange = {0.1, 100.};



  // bundle acceptances into pairs --------------------------------------------

  pair<Types::JetInfo, Types::JetInfo> GetJetAccept() {

    // create maximal range
    pair<Types::JetInfo, Types::JetInfo> jetAccept = {
      Types::JetInfo(Const::Init::Minimize),
      Types::JetInfo(Const::Init::Maximize)
    };

    // set specific bounds
    jetAccept.first.SetPT( ptJetBins.front().first );
    jetAccept.first.SetEta( etaJetRange.first );
    jetAccept.second.SetPT( ptJetBins.back().second );
    jetAccept.second.SetEta( etaJetRange.second );
    return jetAccept;

  }  // end 'GetJetAccept()'



  pair<Types::CstInfo, Types::CstInfo> GetCstAccept() {

   // create maxmimal range
    pair<Types::CstInfo, Types::CstInfo> cstAccept = {
      Types::CstInfo(Const::Init::Minimize),
      Types::CstInfo(Const::Init::Maximize)
    };

    // set specific bounds
    cstAccept.first.SetDR( drCstRange.first );
    cstAccept.first.SetEne( eneCstRange.first );
    cstAccept.second.SetDR( drCstRange.second );
    cstAccept.second.SetEne( eneCstRange.second );
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
