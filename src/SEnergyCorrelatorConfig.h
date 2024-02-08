// ----------------------------------------------------------------------------
// 'SEnergyCorrelatorConfig.h'
// Derek Anderson
// 01.18.2024
//
// Configuration struct for 'SEnergyCorrelator' module.
// ----------------------------------------------------------------------------

#ifndef SENERGYCORRELATORCONFIG_H
#define SENERGYCORRELATORCONFIG_H

// c++ utilities
#include <string>
#include <vector>
#include <utility>
// analysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/CstTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Constants.h"

// make common namespaces implicit
using namespace std;
using namespace SColdQcdCorrelatorAnalysis::SCorrelatorUtilities;



namespace SColdQcdCorrelatorAnalysis {

  // SEnergyCorrelatorConfig definition ---------------------------------------

  struct SEnergyCorrelatorConfig {

    // system options
    int            verbosity     = 0;
    bool           isDebugOn     = false;
    bool           isBatchOn     = false;
    bool           isEmbed       = false;
    bool           isStandalone  = true;
    bool           isInTreeTruth = false;
    string         moduleName    = "SEnergyCorrelator";
    string         inNodeName    = "";
    string         inTreeName    = "";
    string         outFileName   = "";
    vector<string> inFileNames;

    // calculation options
    bool                         applyCstCuts  = false;
    bool                         selectSubEvts = false;
    uint64_t                     nBinsDr       = 75;
    SubEvtOpt                    subEvtOpt     = SubEvtOpt::Everything;
    vector<int>                  subEvtsToUse  = {};
    vector<uint32_t>             nPoints       = {2};
    vector<pair<double, double>> ptJetBins     = {{0., 100.}};
    pair<double, double>         drBinRange    = {1e-5, 1.};

    // acceptance parameters
    pair<JetInfo, JetInfo> jetAccept;
    pair<CstInfo, CstInfo> cstAccept;

  };

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
