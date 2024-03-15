// ----------------------------------------------------------------------------
// 'SEnergyCorrelatorLegacyInput.h'
// Derek Anderson
// 03.15.2023
//
// This struct collects the old branches from the
// input jet trees and methods for converting
// to-from the old-to-new input.
// ----------------------------------------------------------------------------

#ifndef SENERGYCORRELATORLEGACYINPUT_H
#define SENERGYCORRELATORLEGACYINPUT_H

// make common namespaces implicit


namespace SColdQcdCorrelatorAnalysis {

  struct SEnergyCorrelatorLegacyInput {

    // input truth tree addresses
    int                  evtNumChrgPars = -999;
    double               evtSumPar      = -999.;
    pair<int, int>       partonID       = {-999,  -999};
    pair<double, double> partonMomX     = {-999., -999.};
    pair<double, double> partonMomY     = {-999., -999.};
    pair<double, double> partonMomZ     = {-999., -999.};
    vector<vector<int>>* cstID          = NULL;
    vector<vector<int>>* cstEmbedID     = NULL;

    // input reco. tree addresses
    int                  evtNumTrks = -999;
    double               evtSumECal = -999.;
    double               evtSumHCal = -999.;
    vector<vector<int>>* cstMatchID = NULL;

    // generic input tree address members
    int                     evtNumJets = -999;
    double                  evtVtxX    = -999.;
    double                  evtVtxY    = -999.;
    double                  evtVtxZ    = -999.;
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
      evtNumChrgPars = -999;
      evtSumPar      = -999.;
      partonID       = {-999,  -999};
      partonMomX     = {-999., -999.};
      partonMomY     = {-999., -999.};
      partonMomZ     = {-999., -999.};
      cstID          = NULL;
      cstEmbedID     = NULL;
      evtNumTrks     = -999;
      evtSumECal     = -999.;
      evtSumHCal     = -999.;
      cstMatchID     = NULL;
      evtNumJets     = -999;
      evtVtxX        = -999.;
      evtVtxY        = -999.;
      evtVtxZ        = -999.;
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

  };  // end SEnergyCorrelatorLegacyInput

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------

