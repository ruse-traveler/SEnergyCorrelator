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
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <iostream>
// root libraries
#include <TH1.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TDirectory.h>
#include <Math/Vector4D.h>
// phool utilities
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
// f4a utilities
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
// fastjet libraries
#include <fastjet/PseudoJet.hh>
// eec library
#include "/sphenix/user/danderson/eec/EnergyEnergyCorrelators/eec/include/EECLongestSide.hh"
// analysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Tools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Types.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Constants.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/Interfaces.h"
// analysis definitions
#include "SEnergyCorrelatorConfig.h"
#include "SEnergyCorrelatorLegacyInput.h"


// make common namespaces implicit
using namespace std;

// forward declarations
class PHCompositeNode;



namespace SColdQcdCorrelatorAnalysis {

  // SEnergyCorrelator definition ---------------------------------------------

  class SEnergyCorrelator : public SubsysReco {

    public:

      // ctor/dtor
      SEnergyCorrelator(SEnergyCorrelatorConfig& config);
      ~SEnergyCorrelator() override;

      // F4A methods
      int Init(PHCompositeNode*)          override;
      int process_event(PHCompositeNode*) override;
      int End(PHCompositeNode*)           override;

      // standalone-only methods
      void Init();
      void Analyze();
      void End();

      // setters
      void SetConfig(const SEnergyCorrelatorConfig& config) {m_config = config;}

      // getters
      SEnergyCorrelatorConfig GetConfig() {return m_config;}

    private:

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

      // analysis methods (*.ana.h)
      void     DoCorrelatorCalculation();
      void     ExtractHistsFromCorr();
      bool     ApplyJetCuts(const double ptJet, const double etaJet);
      bool     ApplyCstCuts(const double momCst, const double drCst);
      uint32_t GetJetPtBin(const double ptJet);

      // configuration
      SEnergyCorrelatorConfig m_config;

      // io members
      TFile*  m_outFile = NULL;
      TFile*  m_inFile  = NULL;
      TChain* m_inChain = NULL;

      // system members
      int m_fCurrent = 0;

      // output histograms
      vector<TH1D*> m_outHistVarDrAxis;
      vector<TH1D*> m_outHistErrDrAxis;
      vector<TH1D*> m_outHistVarLnDrAxis;
      vector<TH1D*> m_outHistErrLnDrAxis;

      // correlators
      vector<fastjet::contrib::eec::EECLongestSide<fastjet::contrib::eec::hist::axis::log>*> m_eecLongSide;

      // inputs
      SEnergyCorrelatorLegacyInput m_legacy;

      // input truth tree branch members
      //   - FIXME swap out for utlity types when ready
      TBranch*                 m_brEvtSumPar  = NULL;
      TBranch*                 m_brCstID      = NULL;
      TBranch*                 m_brCstEmbedID = NULL;
      pair<TBranch*, TBranch*> m_brPartonID   = {NULL, NULL};
      pair<TBranch*, TBranch*> m_brPartonMomX = {NULL, NULL};
      pair<TBranch*, TBranch*> m_brPartonMomY = {NULL, NULL};
      pair<TBranch*, TBranch*> m_brPartonMomZ = {NULL, NULL};
      // input reco. tree branch members
      TBranch* m_brEvtNumTrks = NULL;
      TBranch* m_brEvtSumECal = NULL;
      TBranch* m_brEvtSumHCal = NULL;
      TBranch* m_brCstMatchID = NULL;

      // generic input tree branch members
      //   - FIXME swap out for utlity types when ready
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

  };  // end SEnergyCorrelator

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
