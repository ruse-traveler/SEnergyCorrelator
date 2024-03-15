// ----------------------------------------------------------------------------
// 'DoStandaloneCorrelatorCalculation.cxx'
// Derek Anderson
// 02.02.2023
//
// Use this to run the SEnergyCorrelator class in standalone mode.
// ----------------------------------------------------------------------------

#ifndef DOSTANDALONECORRELATORCALCULATION_CXX
#define DOSTANDALONECORRELATORCALCULATION_CXX

// c++ utilities
#include <string>
#include <vector>
#include <cstdlib>
#include <utility>
// analysis utilities
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelator.h"
#include "/sphenix/user/danderson/install/include/senergycorrelator/SEnergyCorrelatorConfig.h"
#include "EnergyCorrelatorOptions.h"

using namespace std;
using namespace SColdQcdCorrelatorAnalysis;

// load libraries
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libsenergycorrelator.so)
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libscorrelatorutilities.so)



// macro body -----------------------------------------------------------------

void DoStandaloneCorrelatorCalculation() {

  // get module configurations
  SEnergyCorrelatorConfig cfg_reco = EnergyCorrelatorOptions::GetRecoConfig();
  SEnergyCorrelatorConfig cfg_true = EnergyCorrelatorOptions::GetTruthConfig();

  // do correlator calculation on reco jets
  SEnergyCorrelator* recoCorrelator = new SEnergyCorrelator(cfg_reco);
  recoCorrelator -> Init();
  recoCorrelator -> Analyze();
  recoCorrelator -> End();

  // do correlator calculation on truth jets
  SEnergyCorrelator* trueCorrelator = new SEnergyCorrelator(cfg_true);
  trueCorrelator -> Init();
  trueCorrelator -> Analyze();
  trueCorrelator -> End();
  return;

}

#endif

// end ------------------------------------------------------------------------
