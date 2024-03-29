// ----------------------------------------------------------------------------
// 'Fun4All_RunCorrelatorJetTree.C'
// Derek Anderson
// 12.11.2022
//
// Use this to run the SCorrelatorJetTree
// class.
//
// Derived from code by Cameron Dean and
// Antonio Silva (thanks!!)
//
// NOTE: jetType sets whether or not jets
// are full (charge + neutral) or charged
//   jetType = 0: charged jets
//   jetType = 1: full jets
// ----------------------------------------------------------------------------

/****************************/
/*     MDC2 Reco for MDC2   */
/* Cameron Dean, LANL, 2021 */
/*      cdean@bnl.gov       */
/****************************/

// standard c includes
#include <vector>
#include <string>
#include <cstdlib>
#include <utility>
// f4a/sphenix includes
#include <QA.C>
#include <FROG.h>
#include <G4_Magnet.C>
#include <fun4all/Fun4AllDstInputManager.h>
#include <g4main/Fun4AllDstPileupInputManager.h>
// tracking includes
#include <Trkr_QA.C>
#include <Trkr_Reco.C>
#include <Trkr_Eval.C>
#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_Diagnostics.C>
#include <G4_TrkrSimulation.C>
#include <g4eval/SvtxEvaluator.h>
#include <g4eval/SvtxTruthRecoTableEval.h>
// calo/pf includes
#include <caloreco/RawClusterBuilderTopo.h>
#include <particleflowreco/ParticleFlowReco.h>
// user includes
#include "/sphenix/user/danderson/install/include/scorrelatorjettree/SCorrelatorJetTree.h"
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelator.h"

// load libraries
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libparticleflow.so)
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libscorrelatorjettree.so)
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libsenergycorrelator.so)

using namespace std;
using namespace SColdQcdCorrelatorAnalysis;

// global constants
static const int            NEvtDefault = 10;
static const int            VerbDefault = 0;
static const size_t         NTopoClusts = 2;
static const size_t         NTopoPar    = 3;
static const string         SOutDefault = "testingPAuInput.root";
static const vector<string> SInDefault  = {
  "DST_GLOBAL_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_TRKR_G4HIT_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_TRACKSEEDS_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_TRKR_CLUSTER_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_TRACKS_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_CALO_G4HIT_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_CALO_CLUSTER_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_TRUTH_G4HIT_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_TRUTH_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root",
  "DST_VERTEX_pythia8_Jet20_sHijing_pAu_0_10fm_500kHz_bkg_0_10fm-0000000006-00009.root"
};



void Fun4All_DoStandaloneCalculation(const vector<string>& sInput = SInDefault, const string sOutput = SOutDefault, const int nEvents = NEvtDefault, const int verbosity = VerbDefault) {

  // track & particle flow parameters
  const bool   runTracking(false);
  const bool   doTruthTableReco(false);
  const double nSigma(1.5);

  // topo cluster parameters
  const double showerR(0.025);
  const double noiseLevels[NTopoPar]   = {0.0025, 0.006, 0.03};
  const double significance[NTopoPar]  = {4.0,    2.0,   0.0};
  const double localMinE[NTopoPar]     = {1.0,    2.0,   0.5};
  const bool   enableHCal[NTopoClusts] = {false, true};
  const bool   enableECal[NTopoClusts] = {true, false};
  const bool   doSplit(true);
  const bool   allowCorners(true);

  // jet tree general parameters
  const bool isMC(true);
  const bool doDebug(true);
  const bool saveDst(true);
  const bool doQuality(true);
  const bool addTracks(true);
  const bool addECal(false);
  const bool addHCal(false);
  const bool addParticleFlow(false);

  // jet tree jet parameters
  const double       jetRes  = 0.4;
  const unsigned int jetType = 0;
  const auto         jetAlgo = SCorrelatorJetTree::ALGO::ANTIKT;
  const auto         jetReco = SCorrelatorJetTree::RECOMB::PT_SCHEME;

  // constituent acceptance
  const pair<double, double> ptParRange    = {0.,   9999.};
  const pair<double, double> etaParRange   = {-1.1, 1.1};
  const pair<double, double> ptTrackRange  = {0.2,  9999.};
  const pair<double, double> etaTrackRange = {-1.1, 1.1};
  const pair<double, double> ptFlowRange   = {0.2,  9999.};
  const pair<double, double> etaFlowRange  = {-1.1, 1.1};
  const pair<double, double> ptECalRange   = {0.3,  9999.};
  const pair<double, double> etaECalRange  = {-1.1, 1.1};
  const pair<double, double> ptHCalRange   = {0.3,  9999.};
  const pair<double, double> etaHCalRange  = {-1.1, 1.1};

  // correlator io parameters
  const vector<string> inTree = {
    "RecoJetTree",
    "TruthJetTree"
  };
  const vector<string> outFile = {
    "pa200hijing50khzrun6jet20.reco.d22m8y2023.root",
    "pa200hijing50khzrun6jet20.true.d22m8y2023.root"
  };

  // correlator module parameters
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

  // correlator jet/cst parameters
  const pair<double, double> etaJetRange = {-1., 1.};
  const pair<double, double> momCstRange = {0.,  100.};
  const pair<double, double> drCstRange  = {0.,  5.};

  // correlator jet pT bins
  const vector<pair<double, double>> ptJetBins = {{5., 10.}, {10., 15.}, {15., 20.}, {20., 30.}, {30., 50.}};

  // correlator misc parameters
  const bool isComplex = false;
  const bool inBatch   = false;

  // load libraries and create f4a server
  gSystem -> Load("libg4dst.so");
  gSystem -> Load("libFROG.so");

  FROG*          frog      = new FROG();
  Fun4AllServer* ffaServer = Fun4AllServer::instance();
  ffaServer -> Verbosity(verbosity);

  // add input files 
  for (size_t iInput = 0; iInput < sInput.size(); iInput++) {
    Fun4AllDstInputManager* inManager = new Fun4AllDstInputManager("InputDstManager" + to_string(iInput));
    inManager -> AddFile(sInput.at(iInput));
    ffaServer -> registerInputManager(inManager);
  }

  // run the tracking if not already done
  if (runTracking) {

    // enable mms
    Enable::MICROMEGAS = true;

    // initialize magnetic field
    G4MAGNET::magfield_rescale = 1.;
    MagnetInit();
    MagnetFieldInit();

    // initialize tracker cells
    Mvtx_Cells();
    Intt_Cells();
    TPC_Cells();
    Micromegas_Cells();

    // initialize tracking 
    TrackingInit();

    // do tracker clustering & reconstruction
    Mvtx_Clustering();
    Intt_Clustering();
    TPC_Clustering();
    Micromegas_Clustering();
    Tracking_Reco();
  }

  // construct track/truth table
  if (doTruthTableReco) {
    SvtxTruthRecoTableEval *tables = new SvtxTruthRecoTableEval();
    tables -> Verbosity(verbosity);
    if (runTracking) {
      ffaServer -> registerSubsystem(tables);
    }
  }

  // if using particle flow, run pf reconstruction
  if (addParticleFlow) {

    // build topo clusters
    RawClusterBuilderTopo* ecalClusterBuilder = new RawClusterBuilderTopo("EcalRawClusterBuilderTopo");
    ecalClusterBuilder -> Verbosity(verbosity);
    ecalClusterBuilder -> set_nodename("TOPOCLUSTER_EMCAL");
    ecalClusterBuilder -> set_enable_HCal(enableHCal[0]);
    ecalClusterBuilder -> set_enable_EMCal(enableECal[0]);
    ecalClusterBuilder -> set_noise(noiseLevels[0], noiseLevels[1], noiseLevels[2]);
    ecalClusterBuilder -> set_significance(significance[0], significance[1], significance[2]);
    ecalClusterBuilder -> allow_corner_neighbor(allowCorners);
    ecalClusterBuilder -> set_do_split(doSplit);
    ecalClusterBuilder -> set_minE_local_max(localMinE[0], localMinE[1], localMinE[2]);
    ecalClusterBuilder -> set_R_shower(showerR);
    ffaServer          -> registerSubsystem(ecalClusterBuilder);

    RawClusterBuilderTopo* hcalClusterBuilder = new RawClusterBuilderTopo("HcalRawClusterBuilderTopo");
    hcalClusterBuilder -> Verbosity(verbosity);
    hcalClusterBuilder -> set_nodename("TOPOCLUSTER_HCAL");
    hcalClusterBuilder -> set_enable_HCal(enableHCal[1]);
    hcalClusterBuilder -> set_enable_EMCal(enableECal[1]);
    hcalClusterBuilder -> set_noise(noiseLevels[0], noiseLevels[1], noiseLevels[2]);
    hcalClusterBuilder -> set_significance(significance[0], significance[1], significance[1]);
    hcalClusterBuilder -> allow_corner_neighbor(allowCorners);
    hcalClusterBuilder -> set_do_split(doSplit);
    hcalClusterBuilder -> set_minE_local_max(localMinE[0], localMinE[1], localMinE[2]);
    hcalClusterBuilder -> set_R_shower(showerR);
    ffaServer          -> registerSubsystem(hcalClusterBuilder);

    // do particle flow
    ParticleFlowReco *parFlowReco = new ParticleFlowReco();
    parFlowReco -> set_energy_match_Nsigma(nSigma);
    parFlowReco -> Verbosity(verbosity);
    ffaServer   -> registerSubsystem(parFlowReco);
  }

  // create correlator jet tree
  SCorrelatorJetTree *correlatorJetTree = new SCorrelatorJetTree("SCorrelatorJetTree", sOutput, isMC, doDebug);
  correlatorJetTree -> Verbosity(verbosity);
  correlatorJetTree -> SetDoQualityPlots(doQuality);
  correlatorJetTree -> SetAddTracks(addTracks);
  correlatorJetTree -> SetAddFlow(addParticleFlow);
  correlatorJetTree -> SetAddECal(addECal);
  correlatorJetTree -> SetAddHCal(addHCal);
  correlatorJetTree -> SetParPtRange(ptParRange);
  correlatorJetTree -> SetParEtaRange(etaParRange);
  correlatorJetTree -> SetTrackPtRange(ptTrackRange);
  correlatorJetTree -> SetTrackEtaRange(etaTrackRange);
  correlatorJetTree -> SetFlowPtRange(ptFlowRange);
  correlatorJetTree -> SetFlowEtaRange(etaFlowRange);
  correlatorJetTree -> SetECalPtRange(ptECalRange);
  correlatorJetTree -> SetECalEtaRange(etaECalRange);
  correlatorJetTree -> SetHCalPtRange(ptHCalRange);
  correlatorJetTree -> SetHCalEtaRange(etaHCalRange);
  correlatorJetTree -> SetJetParameters(jetRes, jetType, jetAlgo, jetReco);
  correlatorJetTree -> SetSaveDST(saveDst);
  ffaServer         -> registerSubsystem(correlatorJetTree);

  // run reconstruction & close f4a
  ffaServer -> run(nEvents);
  ffaServer -> End();
  delete ffaServer;

  // do correlator calculation on reco jets
  SEnergyCorrelator *recoCorrelator = new SEnergyCorrelator(moduleName[0], isComplex, doDebug, inBatch);
  recoCorrelator -> SetVerbosity(verbosity);
  recoCorrelator -> SetInputFile(sOutput);
  recoCorrelator -> SetInputTree(inTree[0], isTruth[0]);
  recoCorrelator -> SetOutputFile(outFile[0]);
  recoCorrelator -> SetJetParameters(ptJetBins, etaJetRange.first, etaJetRange.second);
  recoCorrelator -> SetConstituentParameters(momCstRange.first, momCstRange.second, drCstRange.first, drCstRange.second);
  recoCorrelator -> SetCorrelatorParameters(nPointCorr, nBinsDr, binRangeDr.first, binRangeDr.second);
  recoCorrelator -> Init();
  recoCorrelator -> Analyze();
  recoCorrelator -> End();

  // do correlator calculation on truth jets
  SEnergyCorrelator *trueCorrelator = new SEnergyCorrelator(moduleName[1], isComplex, doDebug, inBatch);
  trueCorrelator -> SetVerbosity(verbosity);
  trueCorrelator -> SetInputFile(sOutput);
  trueCorrelator -> SetInputTree(inTree[1], isTruth[1]);
  trueCorrelator -> SetOutputFile(outFile[1]);
  trueCorrelator -> SetJetParameters(ptJetBins, etaJetRange.first, etaJetRange.second);
  trueCorrelator -> SetConstituentParameters(momCstRange.first, momCstRange.second, drCstRange.first, drCstRange.second);
  trueCorrelator -> SetCorrelatorParameters(nPointCorr, nBinsDr, binRangeDr.first, binRangeDr.second);
  trueCorrelator -> Init();
  trueCorrelator -> Analyze();
  trueCorrelator -> End();

  // announce end & exit
  gSystem -> Exit(0);
  return;

}

// end ------------------------------------------------------------------------
