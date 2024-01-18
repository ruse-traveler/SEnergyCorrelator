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
#include "/sphenix/user/danderson/eec/SCorrelatorUtilities/CstInfo.h"
#include "/sphenix/user/danderson/eec/SCorrelatorUtilities/JetTools.h"

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
    uint16_t                     subEvtOpt     = 0;
    uint64_t                     nBinsDr;
    vector<int>                  subEvtsToUse;
    vector<size_t>               nBinsJetPt;
    vector<uint32_t>             nPoints;
    vector<pair<double, double>> ptJetBins;
    pair<double, double>         drBinRange;

    // acceptance parameters
    pair<JetInfo, JetInfo> jetAccept;
    pair<CstInfo, CstInfo> cstAccept;

  };

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
