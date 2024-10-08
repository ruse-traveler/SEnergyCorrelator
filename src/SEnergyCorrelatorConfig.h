/// ---------------------------------------------------------------------------
/*! \file    SEnergyCorrelatorConfig.h
 *  \authors Derek Anderson, Alex Clarke
 *  \date    01.18.2024
 *
 *  Configuration struct for 'SEnergyCorrelator' module.
 */
/// ---------------------------------------------------------------------------

#ifndef SENERGYCORRELATORCONFIG_H
#define SENERGYCORRELATORCONFIG_H

// make common namespaces implicit
using namespace std;



namespace SColdQcdCorrelatorAnalysis {

  // --------------------------------------------------------------------------
  //! User options for module
  // --------------------------------------------------------------------------
  struct SEnergyCorrelatorConfig {

    // general correlator parameters
    vector<uint32_t>     nPoints    { {2} };
    pair<double, double> drBinRange {1e-5, 1.};
    uint64_t             nBinsDr    {75};

    // basic options
    int  verbosity    {0};
    int  subEvtOpt    {Const::SubEvtOpt::Everything};
    bool isEmbed      {false};
    bool applyCstCuts {false};
    bool isDebugOn    {false};
    bool isBatchOn    {false};
    bool doPackageCalc {true};

    // manual calculation options
    int  momOption     {0};
    int  normOption    {0};
    bool doManualCalc  {false};
    bool doThreePoint  {false};

    // i/o options
    bool           isInTreeTruth {false};
    string         moduleName    {"SEnergyCorrelator"};
    string         inTreeName    {""};
    string         outFileName   {""};
    vector<string> inFileNames   { {} };

    // smearing options
    //   - FIXME theta should be thetaSmear
    bool  modCsts    {false};
    bool  modJets    {false};
    float effVal     {1};
    float theta      {0};
    float ptCstSmear {0};
    float ptJetSmear {0};

    // subevent options
    bool        selectSubEvts {false};
    vector<int> subEvtsToUse  { {} };

    // histogram bins
    vector<pair<double, double>> rlBins    { {{0.0, 1.0}} };
    vector<pair<double, double>> ptJetBins { {{0., 100.}} };

    // acceptance parameters
    pair<Types::JetInfo, Types::JetInfo> jetAccept;
    pair<Types::CstInfo, Types::CstInfo> cstAccept;

  };

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
