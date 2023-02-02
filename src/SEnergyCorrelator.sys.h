// 'SEnergyCorrelator.sys.h'
// Derek Anderson
// 01.27.2023
//
// A module to implement Peter Komiske's
// EEC library in the sPHENIX software
// stack.

#pragma once

using namespace std;
using namespace findNode;



void SEnergyCorrelator::InitializeMembers() {

  // print debug statement
  if (m_inDebugMode) {
    cout << "SEnergyCorrelator::InitializeMembers() initializing internal variables..." << endl;
  }

  m_inFile           = 0x0;
  m_inTree           = 0x0;
  m_outFile          = 0x0;
  m_verbosity        = 0;
  m_inDebugMode      = false;
  m_inComplexMode    = false;
  m_inStandaloneMode = false;
  m_inFileName       = "";
  m_inNodeName       = "";
  m_inTreeName       = "";
  m_outFileName      = "";
  return;

}  // end 'InitializeMembers()'



void SEnergyCorrelator::InitializeHists() {

  // print debug statement
  if (m_inDebugMode) {
     cout << "SEnergyCorrelator::InitializeHists() initializing histograms..." << endl;
  }

  /* output histograms wil be initialized here */
  return;

}  // end 'InitializeHists()'



void SEnergyCorrelator::InitializeCorrs() {

  /* correlators will be initialized here */
  return;

}  // end 'InitializeCorrs()'



void SEnergyCorrelator::PrintMessage() {

  // print debug statement
  if (m_inDebugMode && (m_verbosity > 5)) {
    cout << "SEnergyCorrelator::PrintMessage() printing a message..." << endl;
  }
  return;

}  // end 'PrintMessage()'



void SEnergyCorrelator::PrintError(const uint32_t code) {

  // print debug statement
  if (m_inDebugMode && (m_verbosity > 5)) {
    cout << "SEnergyCorrelator::PrintError() printing an error..." << endl;
  }

  switch (code) {
    case 0:
      cerr
      break;
  }
  return;

}  // end 'PrintError()'

// end ------------------------------------------------------------------------
