// ----------------------------------------------------------------------------
// 'SEnergyCorrelator.h'
// Derek Anderson
// 01.20.2023
//
// A module to implement Peter Komiske's EEC library
// in the sPHENIX software stack for the Cold QCD
// Energy-Energy Correlator analysis.
// ----------------------------------------------------------------------------

#ifndef SENERGYCORRELATOR_H
#define SENERGYCORRELATOR_H

// c++ utilities
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
#include <cstdlib>
#include <utility>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TString.h>
#include <TDirectory.h>
#include <Math/Vector4D.h>
// f4a utilities
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
// phool utilities
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
// fastjet libraries
#include <fastjet/PseudoJet.hh>
// eec library
#include "/sphenix/user/danderson/eec/EnergyEnergyCorrelators/eec/include/EECLongestSide.hh"
// analysis libraries
#include "SEnergyCorrealtorConfig.h"

using namespace std;
using namespace fastjet;



namespace SColdQcdCorrelatorAnalysis {

  // SEnergyCorrelator definition ---------------------------------------------

  class SEnergyCorrelator : public SubsysReco {

    public:

      // ctor/dtor
      SEnergyCorrelator(const string &name = "SEnergyCorrelator", const bool isComplex = false, const bool doDebug = false, const bool inBatch = false);
      ~SEnergyCorrelator() override;

      // F4A methods
      int Init(PHCompositeNode*)          override;
      int process_event(PHCompositeNode*) override;
      int End(PHCompositeNode*)           override;

      // standalone-only methods
      void Init();
      void Analyze();
      void End();

      // setters (inline)
      void SetConfig(const SEnergyCorrelatorConfig& config) {m_config = config;}
      void SetDoSecondCstLoop(const bool loop)              {m_doSecondCstLoop = loop;}  // FIXME remove when ready

      // system getters
      bool GetIsInputTreeTruth() {return m_isInputTreeTruth;}
      bool GetIsInputTreeEmbed() {return m_isInputTreeEmbed;}
      bool GetDoSecondCstLoop()  {return m_doSecondCstLoop;}

      // correlator getters
      SEnergyCorrelatorConfig GetConfig() {return m_config;}

    private:

      // constants
      enum CONST {
        NPARTONS = 2
      };

      // io methods (*.io.h)
      void GrabInputNode();
      void OpenInputFiles();
      void OpenOutputFile();
      void SaveOutput();
      void CloseOutputFile();

      // system methods (*.sys.h)
      void    InitializeMembers();
      void    InitializeHists();
      void    InitializeCorrs();
      void    InitializeTree();
      void    PrintMessage(const uint32_t code, const uint64_t nEvts = 0, const uint64_t event = 0);
      void    PrintDebug(const uint32_t code);
      void    PrintError(const uint32_t code, const size_t nDrBinEdges = 0, const size_t iDrBin = 0, const string sInFileName = "");
      bool    CheckCriticalParameters();
      int64_t LoadTree(const uint64_t entry);
      int64_t GetEntry(const uint64_t entry);

      // analysis methods (*.ana.h)
      void     DoCorrelatorCalculation();
      void     ExtractHistsFromCorr();
      bool     ApplyJetCuts(const double ptJet, const double etaJet);
      bool     ApplyCstCuts(const double momCst, const double drCst);
      bool     CheckIfSubEvtGood(const int embedID);
      uint32_t GetJetPtBin(const double ptJet);

      // configuration
      SEnergyCorrelatorConfig m_config;

      // io members
      TFile*  m_outFile = NULL;
      TFile*  m_inFile  = NULL;
      TChain* m_inChain = NULL;

      // output histograms
      vector<TH1D*> m_outHistVarDrAxis;
      vector<TH1D*> m_outHistErrDrAxis;
      vector<TH1D*> m_outHistVarLnDrAxis;
      vector<TH1D*> m_outHistErrLnDrAxis;

      // for weird cst check
      TH2D* hCstPtOneVsDr;
      TH2D* hCstPtTwoVsDr;
      TH2D* hCstPtFracVsDr;
      TH2D* hCstPhiOneVsDr;
      TH2D* hCstPhiTwoVsDr;
      TH2D* hCstEtaOneVsDr;
      TH2D* hCstEtaTwoVsDr;
      TH2D* hDeltaPhiOneVsDr;
      TH2D* hDeltaPhiTwoVsDr;
      TH2D* hDeltaEtaOneVsDr;
      TH2D* hDeltaEtaTwoVsDr;
      TH2D* hJetPtFracOneVsDr;
      TH2D* hJetPtFracTwoVsDr;
      TH2D* hCstPairWeightVsDr;

      // system members
      int  m_fCurrent         = 0;
      bool m_isInputTreeTruth = false;
      bool m_isInputTreeEmbed = false;
      bool m_doSecondCstLoop  = false;

      // vector of input files (standalone mode only)
      vector<string> m_inFileNames;

      // jet, cst, correlator parameters
      uint32_t                     m_nPointCorr  = 0;
      uint64_t                     m_nBinsDr     = 0;
      size_t                       m_nBinsJetPt  = 0;
      pair<double, double>         m_drBinRange  = {-999., -999.};
      pair<double, double>         m_ptJetRange  = {-999., -999.};
      pair<double, double>         m_etaJetRange = {-999., -999.};
      pair<double, double>         m_momCstRange = {-999., -999.};
      pair<double, double>         m_drCstRange  = {-999., -999.};
      vector<pair<double, double>> m_ptJetBins;
      vector<PseudoJet>            m_jetCstVector;
      vector<int>                  m_subEvtsToUse;

      // correlators
      vector<contrib::eec::EECLongestSide<contrib::eec::hist::axis::log>*> m_eecLongSide;

      // input truth tree address members
      int                  m_evtNumChrgPars              = -999;
      int                  m_partonID[CONST::NPARTONS]   = {-999,  -999};
      double               m_partonMomX[CONST::NPARTONS] = {-999., -999.};
      double               m_partonMomY[CONST::NPARTONS] = {-999., -999.};
      double               m_partonMomZ[CONST::NPARTONS] = {-999., -999.};
      double               m_evtSumPar                   = -999.;
      vector<vector<int>>* m_cstID                       = NULL;
      vector<vector<int>>* m_cstEmbedID                  = NULL;
      // input reco. tree address members
      int                  m_evtNumTrks = -999;
      double               m_evtSumECal = -999.;
      double               m_evtSumHCal = -999.;
      vector<vector<int>>* m_cstMatchID = NULL;

      // generic input tree address members
      int                     m_evtNumJets = -999;
      double                  m_evtVtxX    = -999.;
      double                  m_evtVtxY    = -999.;
      double                  m_evtVtxZ    = -999.;
      vector<unsigned long>*  m_jetNumCst  = NULL;
      vector<unsigned int>*   m_jetID      = NULL;
      vector<unsigned int>*   m_jetTruthID = NULL;
      vector<double>*         m_jetEnergy  = NULL;
      vector<double>*         m_jetPt      = NULL;
      vector<double>*         m_jetEta     = NULL;
      vector<double>*         m_jetPhi     = NULL;
      vector<double>*         m_jetArea    = NULL;
      vector<vector<double>>* m_cstZ       = NULL;
      vector<vector<double>>* m_cstDr      = NULL;
      vector<vector<double>>* m_cstEnergy  = NULL;
      vector<vector<double>>* m_cstPt      = NULL;
      vector<vector<double>>* m_cstEta     = NULL;
      vector<vector<double>>* m_cstPhi     = NULL;

      // input truth tree branch members
      TBranch* m_brPartonID[CONST::NPARTONS]   = {NULL, NULL};
      TBranch* m_brPartonMomX[CONST::NPARTONS] = {NULL, NULL};
      TBranch* m_brPartonMomY[CONST::NPARTONS] = {NULL, NULL};
      TBranch* m_brPartonMomZ[CONST::NPARTONS] = {NULL, NULL};
      TBranch* m_brEvtSumPar                   = NULL;
      TBranch* m_brCstID                       = NULL;
      TBranch* m_brCstEmbedID                  = NULL;
      // input reco. tree branch members
      TBranch* m_brEvtNumTrks = NULL;
      TBranch* m_brEvtSumECal = NULL;
      TBranch* m_brEvtSumHCal = NULL;
      TBranch* m_brCstMatchID = NULL;

      // generic input tree branch members
      TBranch* m_brEvtNumJets = NULL;
      TBranch* m_brEvtVtxX    = NULL;
      TBranch* m_brEvtVtxY    = NULL;
      TBranch* m_brEvtVtxZ    = NULL;
      TBranch* m_brJetNumCst  = NULL;
      TBranch* m_brJetID      = NULL;
      TBranch* m_brJetEnergy  = NULL;
      TBranch* m_brJetPt      = NULL;
      TBranch* m_brJetEta     = NULL;
      TBranch* m_brJetPhi     = NULL;
      TBranch* m_brJetArea    = NULL;
      TBranch* m_brCstZ       = NULL;
      TBranch* m_brCstDr      = NULL;
      TBranch* m_brCstPt      = NULL;
      TBranch* m_brCstEnergy  = NULL;
      TBranch* m_brCstEta     = NULL;
      TBranch* m_brCstPhi     = NULL;

  };

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
