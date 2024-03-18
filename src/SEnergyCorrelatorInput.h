// ----------------------------------------------------------------------------
// 'SEnergyCorrelatorInput.h'
// Derek Anderson
// 03.15.2023
//
// This struct collects the old branches from the
// input jet trees and methods for converting
// to-from the old-to-new input.
// ----------------------------------------------------------------------------

#ifndef SENERGYCORRELATORINPUT_H
#define SENERGYCORRELATORINPUT_H

// make common namespaces implicit
using namespace std;



namespace SColdQcdCorrelatorAnalysis {

  // input definition ---------------------------------------------------------

  struct SEnergyCorrelatorInput {

    // event level info
    Types::RecoInfo m_reco;
    Types::GenInfo  m_gen;

    // jet and constituent info
    vector<Types::JetInfo>         m_jets;
    vector<vector<Types::CstInfo>> m_csts; 

    void Reset() {
      m_reco.Reset();
      m_gen.Reset();
      m_jets.clear();
      m_csts.clear();
      return;
    }  // end 'Reset()'

  };  // end SEnergyCorrelatorInput 



  // legacy input definition --------------------------------------------------

  struct SEnergyCorrelatorLegacyInput {

    // input truth tree addresses
    int                  evtNumChrgPars = numeric_limits<int>::min();
    double               evtSumPar      = numeric_limits<double>::min();
    pair<int, int>       partonID       = {numeric_limits<int>::min(),    numeric_limits<int>::min()};
    pair<double, double> partonMomX     = {numeric_limits<double>::min(), numeric_limits<double>::min()};
    pair<double, double> partonMomY     = {numeric_limits<double>::min(), numeric_limits<double>::min()};
    pair<double, double> partonMomZ     = {numeric_limits<double>::min(), numeric_limits<double>::min()};
    vector<vector<int>>* cstID          = NULL;
    vector<vector<int>>* cstEmbedID     = NULL;

    // input reco. tree addresses
    int                  evtNumTrks = numeric_limits<int>::min();
    double               evtSumECal = numeric_limits<double>::min();
    double               evtSumHCal = numeric_limits<double>::min();
    vector<vector<int>>* cstMatchID = NULL;

    // generic input tree address members
    int                     evtNumJets = numeric_limits<int>::min();
    double                  evtVtxX    = numeric_limits<double>::min();
    double                  evtVtxY    = numeric_limits<double>::min();
    double                  evtVtxZ    = numeric_limits<double>::min();
    vector<unsigned long>*  jetNumCst  = NULL;
    vector<unsigned int>*   jetID      = NULL;
    vector<unsigned int>*   jetTruthID = NULL;
    vector<double>*         jetEnergy  = NULL;
    vector<double>*         jetPt      = NULL;
    vector<double>*         jetEta     = NULL;
    vector<double>*         jetPhi     = NULL;
    vector<double>*         jetArea    = NULL;
    vector<vector<double>>* cstZ       = NULL;
    vector<vector<double>>* cstDr      = NULL;
    vector<vector<double>>* cstEnergy  = NULL;
    vector<vector<double>>* cstPt      = NULL;
    vector<vector<double>>* cstEta     = NULL;
    vector<vector<double>>* cstPhi     = NULL;

    void Reset() {
      evtNumChrgPars = numeric_limits<int>::min();
      evtSumPar      = numeric_limits<double>::min();
      partonID       = make_pair(numeric_limits<int>::min(),    numeric_limits<int>::min());
      partonMomX     = make_pair(numeric_limits<double>::min(), numeric_limits<double>::min());
      partonMomY     = make_pair(numeric_limits<double>::min(), numeric_limits<double>::min());
      partonMomZ     = make_pair(numeric_limits<double>::min(), numeric_limits<double>::min());
      cstID          = NULL;
      cstEmbedID     = NULL;
      evtNumTrks     = numeric_limits<int>::min();
      evtSumECal     = numeric_limits<double>::min();
      evtSumHCal     = numeric_limits<double>::min();
      cstMatchID     = NULL;
      evtNumJets     = numeric_limits<int>::min();
      evtVtxX        = numeric_limits<double>::min();
      evtVtxY        = numeric_limits<double>::min();
      evtVtxZ        = numeric_limits<double>::min();
      jetNumCst      = NULL;
      jetID          = NULL;
      jetTruthID     = NULL;
      jetEnergy      = NULL;
      jetPt          = NULL;
      jetEta         = NULL;
      jetPhi         = NULL;
      jetArea        = NULL;
      cstZ           = NULL;
      cstDr          = NULL;
      cstEnergy      = NULL;
      cstPt          = NULL;
      cstEta         = NULL;
      cstPhi         = NULL;
      return;
    }  // end 'Reset()'

    void SetChainAddresses(TChain* chain, const bool isInTreeTruth = true) {
      if (isInTreeTruth) {
        chain -> SetBranchAddress("Parton3_ID",   &partonID.first);
        chain -> SetBranchAddress("Parton4_ID",   &partonID.second);
        chain -> SetBranchAddress("Parton3_MomX", &partonMomX.first);
        chain -> SetBranchAddress("Parton3_MomY", &partonMomY.first);
        chain -> SetBranchAddress("Parton3_MomZ", &partonMomZ.first);
        chain -> SetBranchAddress("Parton4_MomX", &partonMomX.second);
        chain -> SetBranchAddress("Parton4_MomY", &partonMomY.second);
        chain -> SetBranchAddress("Parton4_MomZ", &partonMomZ.second);
        chain -> SetBranchAddress("EvtSumParEne", &evtSumPar);
        chain -> SetBranchAddress("CstID",        &cstID);
        chain -> SetBranchAddress("CstEmbedID",   &cstEmbedID);
      } else {
        chain -> SetBranchAddress("EvtNumTrks",    &evtNumTrks);
        chain -> SetBranchAddress("EvtSumECalEne", &evtSumECal);
        chain -> SetBranchAddress("EvtSumHCalEne", &evtSumHCal);
        chain -> SetBranchAddress("CstMatchID",    &cstMatchID);
      }
      chain -> SetBranchAddress("EvtVtxX",    &evtVtxX);
      chain -> SetBranchAddress("EvtVtxY",    &evtVtxY);
      chain -> SetBranchAddress("EvtVtxZ",    &evtVtxZ);
      chain -> SetBranchAddress("EvtNumJets", &evtNumJets);
      chain -> SetBranchAddress("JetNumCst",  &jetNumCst);
      chain -> SetBranchAddress("JetID",      &jetID);
      chain -> SetBranchAddress("JetEnergy",  &jetEnergy);
      chain -> SetBranchAddress("JetPt",      &jetPt);
      chain -> SetBranchAddress("JetEta",     &jetEta);
      chain -> SetBranchAddress("JetPhi",     &jetPhi);
      chain -> SetBranchAddress("JetArea",    &jetArea);
      chain -> SetBranchAddress("CstZ",       &cstZ);
      chain -> SetBranchAddress("CstDr",      &cstDr);
      chain -> SetBranchAddress("CstEnergy",  &cstEnergy);
      chain -> SetBranchAddress("CstJt",      &cstPt);
      chain -> SetBranchAddress("CstEta",     &cstEta);
      chain -> SetBranchAddress("CstPhi",     &cstPhi);
      return;
    }  // end 'SetChainAddresses(TChain*)'

  };  // end SEnergyCorrelatorLegacyInput

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------

