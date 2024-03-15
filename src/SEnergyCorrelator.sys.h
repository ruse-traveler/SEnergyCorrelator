// ----------------------------------------------------------------------------
// 'SEnergyCorrelator.sys.h'
// Derek Anderson
// 01.27.2023
//
// A module to implement Peter Komiske's EEC library
// in the sPHENIX software stack for the Cold QCD
// Energy-Energy Correlator analysis.
// ----------------------------------------------------------------------------

#pragma once

// make common namespaces implicit
using namespace std;
using namespace fastjet;



namespace SColdQcdCorrelatorAnalysis {

  // system methods -----------------------------------------------------------

  void SEnergyCorrelator::InitializeMembers() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(0);

    m_outHistVarDrAxis.clear();
    m_outHistErrDrAxis.clear();
    m_outHistVarLnDrAxis.clear();
    m_outHistErrLnDrAxis.clear();
    return;

  }  // end 'InitializeMembers()'



  void SEnergyCorrelator::InitializeTree() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(4);

    // check for tree
    if (!m_inChain) {
      PrintError(10);
      assert(m_inChain);
    }
    m_fCurrent = -1;
    m_inChain -> SetMakeClass(1);

    // set truth vs. reco branch addresses
    if (m_config.isInTreeTruth) {
      m_inChain -> SetBranchAddress("Parton3_ID",   &m_legacy.partonID.first,    &m_brPartonID.first);
      m_inChain -> SetBranchAddress("Parton4_ID",   &m_legacy.partonID.second,   &m_brPartonID.second);
      m_inChain -> SetBranchAddress("Parton3_MomX", &m_legacy.partonMomX.first,  &m_brPartonMomX.first);
      m_inChain -> SetBranchAddress("Parton3_MomY", &m_legacy.partonMomY.first,  &m_brPartonMomY.first);
      m_inChain -> SetBranchAddress("Parton3_MomZ", &m_legacy.partonMomZ.first,  &m_brPartonMomZ.first);
      m_inChain -> SetBranchAddress("Parton4_MomX", &m_legacy.partonMomX.second, &m_brPartonMomX.second);
      m_inChain -> SetBranchAddress("Parton4_MomY", &m_legacy.partonMomY.second, &m_brPartonMomY.second);
      m_inChain -> SetBranchAddress("Parton4_MomZ", &m_legacy.partonMomZ.second, &m_brPartonMomZ.second);
      m_inChain -> SetBranchAddress("EvtSumParEne", &m_legacy.evtSumPar,         &m_brEvtSumPar);
      m_inChain -> SetBranchAddress("CstID",        &m_legacy.cstID,             &m_brCstID);
      m_inChain -> SetBranchAddress("CstEmbedID",   &m_legacy.cstEmbedID,        &m_brCstEmbedID);
    } else {
      m_inChain -> SetBranchAddress("EvtNumTrks",    &m_legacy.evtNumTrks, &m_brEvtNumTrks);
      m_inChain -> SetBranchAddress("EvtSumECalEne", &m_legacy.evtSumECal, &m_brEvtSumECal);
      m_inChain -> SetBranchAddress("EvtSumHCalEne", &m_legacy.evtSumHCal, &m_brEvtSumHCal);
      m_inChain -> SetBranchAddress("CstMatchID",    &m_legacy.cstMatchID, &m_brCstMatchID);
    }

    // set generic branch addresses
    m_inChain -> SetBranchAddress("EvtVtxX",    &m_legacy.evtVtxX,    &m_brEvtVtxX);
    m_inChain -> SetBranchAddress("EvtVtxY",    &m_legacy.evtVtxY,    &m_brEvtVtxY);
    m_inChain -> SetBranchAddress("EvtVtxZ",    &m_legacy.evtVtxZ,    &m_brEvtVtxZ);
    m_inChain -> SetBranchAddress("EvtNumJets", &m_legacy.evtNumJets, &m_brEvtNumJets);
    m_inChain -> SetBranchAddress("JetNumCst",  &m_legacy.jetNumCst,  &m_brJetNumCst);
    m_inChain -> SetBranchAddress("JetID",      &m_legacy.jetID,      &m_brJetID);
    m_inChain -> SetBranchAddress("JetEnergy",  &m_legacy.jetEnergy,  &m_brJetEnergy);
    m_inChain -> SetBranchAddress("JetPt",      &m_legacy.jetPt,      &m_brJetPt);
    m_inChain -> SetBranchAddress("JetEta",     &m_legacy.jetEta,     &m_brJetEta);
    m_inChain -> SetBranchAddress("JetPhi",     &m_legacy.jetPhi,     &m_brJetPhi);
    m_inChain -> SetBranchAddress("JetArea",    &m_legacy.jetArea,    &m_brJetArea);
    m_inChain -> SetBranchAddress("CstZ",       &m_legacy.cstZ,       &m_brCstZ);
    m_inChain -> SetBranchAddress("CstDr",      &m_legacy.cstDr,      &m_brCstDr);
    m_inChain -> SetBranchAddress("CstEnergy",  &m_legacy.cstEnergy,  &m_brCstEnergy);
    m_inChain -> SetBranchAddress("CstJt",      &m_legacy.cstPt,      &m_brCstPt);
    m_inChain -> SetBranchAddress("CstEta",     &m_legacy.cstEta,     &m_brCstEta);
    m_inChain -> SetBranchAddress("CstPhi",     &m_legacy.cstPhi,     &m_brCstPhi);

    // announce tree setting
    if (m_config.isStandalone) PrintMessage(2);
    return;

  }  // end 'InitializeTree()'



  void SEnergyCorrelator::InitializeHists() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(5);

    for (size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {
      TH1D* hInitialVarDrAxis   = NULL;
      TH1D* hInitialErrDrAxis   = NULL;
      TH1D* hInitialVarLnDrAxis = NULL;
      TH1D* hInitialErrLnDrAxis = NULL;
      m_outHistVarDrAxis.push_back(hInitialVarDrAxis);
      m_outHistVarLnDrAxis.push_back(hInitialVarLnDrAxis);
      m_outHistErrDrAxis.push_back(hInitialErrDrAxis);
      m_outHistErrLnDrAxis.push_back(hInitialErrLnDrAxis);
    }

    // announce histogram initialization
    if (m_config.isStandalone) PrintMessage(3);
    return;

  }  // end 'InitializeHists()'



  void SEnergyCorrelator::InitializeCorrs() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(6);

    // initialize correlator for each jet pt bin
    for (size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {
      m_eecLongSide.push_back(
        new contrib::eec::EECLongestSide<contrib::eec::hist::axis::log>(
          m_config.nPoints.at(0),  // TODO enable multiple n per calculation
          m_config.nBinsDr,
          {m_config.drBinRange.first, m_config.drBinRange.second}
        )
      );
    }

    // announce correlator initialization
    if (m_config.isStandalone) PrintMessage(4);
    return;

  }  // end 'InitializeCorrs()'



  void SEnergyCorrelator::PrintMessage(const uint32_t code, const uint64_t nEvts, const uint64_t event) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) PrintDebug(22);

    switch (code) {
      case 0:
        cout << "\n  Running standalone correlator calculation...\n"
             << "    Set name & modes:\n"
             << "      module name      = " << m_config.moduleName.data() << "\n"
             << "      complex mode?    = " << !m_config.isStandalone     << "\n"
             << "      standalone mode? = " << m_config.isStandalone      << "\n"
             << "      debug mode?      = " << m_config.isDebugOn         << "\n"
             << "      batch mode?      = " << m_config.isBatchOn
             << endl;
        break;
      case 1:
        cout << "    Opened files:\n"
             << "      output = " << m_config.outFileName.data() << "\n"
             << "      inputs = {"
             << endl;
        for (const string& inFileName : m_config.inFileNames) {
          cout << "        " << inFileName.data() << endl;
        }
        cout << "      }" << endl;
        break;
      case 2:
        cout << "    Initialized input chain:\n"
             << "      tree name = " << m_config.inTreeName.data()
             << endl;
        break;
      case 3:
        cout << "    Initialized output histograms." << endl;
        break;
      case 4:
        cout << "    Initialized correlators." << endl;
        break;
      case 5:
        cout << "    Set correlator parameters:\n"
             << "      n-point = "       << m_config.nPoints.at(0)     << ", number of dR bins = " << m_config.nBinsDr           << "\n"
             << "      dR bin range = (" << m_config.drBinRange.first  << ", "                     << m_config.drBinRange.second << ")"
             << endl;
        break;
      case 6:
        cout << "    Set jet parameters:\n"
             << "      eta range = (" << m_config.jetAccept.first.GetEta() << ", " << m_config.jetAccept.second.GetEta() << ")\n"
             << "      pt range  = (" << m_config.jetAccept.first.GetPT()  << ", " << m_config.jetAccept.second.GetPT()  << ")\n"
             << "    Set pTjet bins:"
             << endl;
        for (uint32_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {
          cout << "      bin[" << iPtBin << "] = (" << m_config.ptJetBins.at(iPtBin).first << ", " << m_config.ptJetBins.at(iPtBin).second << ")" << endl;
        }
        break;
      case 7:
        cout << "    Beginning event loop: " << nEvts << " events to process..." << endl;
        break;
      case 8:
        if (m_config.isBatchOn) {
          cout << "      processing event " << (event + 1) << "/" << nEvts << "..." << endl;
        } else {
          cout << "      processing event " << (event + 1) << "/" << nEvts << "...\r" << flush;
          if ((event + 1) == nEvts) cout << endl;
        }
        break; 
      case 9:
        cout << "    Analysis finished!" << endl;
        break;
      case 10:
        cout << "    Saved output histograms." << endl;
        break;
      case 11:
        cout << "  Finished correlator calculation!\n" << endl;
        break;
      case 12:
        cout << "    Set constituent parameters:\n"
             << "      apply constituent cuts? = " << m_config.applyCstCuts   << "\n"
             << "      momentum range = ("         << m_config.cstAccept.first.GetPT() << ", " << m_config.cstAccept.second.GetPT() << ")\n"
             << "      dr range       = ("         << m_config.cstAccept.first.GetDR() << ", " << m_config.cstAccept.second.GetDR() << ")"
             << endl;
        break;
      case 13:
        cout << "    Finished event loop!" << endl;
        break;
      case 14:
        cout << "    Extracted output histograms from correlators." << endl;
        break;
      case 15:
        cout << "    Set which sub-events to use:" << endl;
        switch (m_config.subEvtOpt) {
          case 1:
            cout << "      Option " << m_config.subEvtOpt << ": use only signal event" << endl;
            break;
          case 2:
            cout << "      Option " << m_config.subEvtOpt << ": use only background events" << endl;
            break;
          case 3:
            cout << "      Option " << m_config.subEvtOpt << ": use only primary background event" << endl;
            break;
          case 4:
            cout << "      Option " << m_config.subEvtOpt << ": use only pileup events" << endl;
            break;
          case 5:
            cout << "      Option " << m_config.subEvtOpt << ": use events only with these embedding IDs: ";
            for (size_t iEvtToUse = 0; iEvtToUse < m_config.subEvtsToUse.size(); iEvtToUse++) {
              cout << m_config.subEvtsToUse[iEvtToUse];
              if ((iEvtToUse + 1) < m_config.subEvtsToUse.size()) {
                cout << ", ";
              } else {
                cout << endl;
              }
            }  // end sub-event id loop
            break;
          default:
            cout << "     Option " << m_config.subEvtOpt << ": use everything (check what you entered)" << endl;
            break;
        }
        break;
    }
    return;

  }  // end 'PrintMessage(uint32_t)'


  void SEnergyCorrelator::PrintDebug(const uint32_t code) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) {
      cout << "SEnergyCorrelator::PrintDebug(uint32_t) printing a debugging statement..." << endl;
    }

    switch (code) {
      case 0:
        cout << "SEnergyCorrelator::InitializeMembers() initializing internal variables..." << endl;
        break;
      case 1:
        cout << "SEnergyCorrelator::SEnergyCorrelator(string, bool, bool) calling ctor..." << endl;
        break;
      case 2:
        cout << "SEnergyCorrelator::Init(PHCompositeNode*) initializing..." << endl;
        break;
      case 3:
        cout << "SEnergyCorrelator::GrabInputNode() grabbing input node..." << endl;
        break;
      case 4:
        cout << "SEnergyCorrelator::InitializeTree() initializing input tree..." << endl;
        break;
      case 5:
        cout << "SEnergyCorrelator::InitializeHists() initializing histograms..." << endl;
        break;
      case 6:
        cout << "SEnergyCorrelator::InitializeCorrs() initializing correlators" << endl;
        break;
      case 7:
        cout << "SEnergyCorrelator::process_event(PHCompositeNode*) processing event..." << endl;
        break;
      case 8:
        cout << "SEnergyCorrelator::End(PHCompositeNode*) this is the end..." << endl;
        break;
      case 9:
        cout << "SEnergyCorrelator::SaveOutput() saving output..." << endl;
        break;
      case 10:
        cout << "SEnergyCorrelator::Init() initializing..." << endl;
        break;
      case 11:
        cout << "SenergyCorrelator::OpenInputFile() opening input file..." << endl;
        break;
      case 12:
        cout << "SEnergyCorrelator::Analyze() analyzing input..." << endl;
        break;
      case 13:
        cout << "SEnergyCorrelator::End() this is the end..." << endl;
        break;
      case 14:
        cout << "SEnergyCorrelator::~SEnergyCorrelator() calling dtor..." << endl;
        break;
      case 15:
        cout << "SEnergyCorrelator::OpenOutputFile() opening output file..." << endl;
        break;
      case 16:
        cout << "SEnergyCorrelator::GetEntry(uint64_t) getting tree entry..." << endl;
        break;
      case 17:
        cout << "SEnergyCorrelator::LoadTree(uint64_t) loading tree..." << endl;
        break;
      case 18:
        cout << "SEnergyCorrelator::SetInputTree(string, bool) setting input tree name..." << endl;
        break;
      case 19:
        cout << "SEnergyCorrelator::SetCorrelatorParameters(uint32_t, uint64_t, pair<double, double>) setting correlator parameters..." << endl;
        break;
      case 20:
        cout << "SEnergyCorrelator::SetJetParameters(vector<pair<double, double>>, pair<double, double>) setting jet parameters..." << endl;
        break;
      case 21:
        cout << "SEnergyCorrelators:CheckCriticalParameters() checking critical parameters..." << endl;
        break;
      case 22:
        cout << "SEnergyCorrelator::PrintMessage(uint32_t, uint64_t, uint64_t) printing a message..." << endl;
        break;
      case 23:
        cout << "SEnergyCorrelator::PrintError(uint32_t) printing an error..." << endl;
        break;
      case 24:
        cout << "SEnergyCorrelator::SetConstituentParameters(pair<double, double>, pair<double, double>) setting constituent parameters..." << endl;
        break;
      case 25:
        cout << "SEnergyCorrelator::ExtractHistsFromCorr() extracting output histograms..." << endl;
        break;
      case 26:
        cout << "SEnergyCorrelator::ApplyJetCuts(double, double) applying jet cuts..." << endl;
        break;
      case 27:
        cout << "SEnergyCorrelator::ApplyCstCuts(double, double) applying constituent cuts..." << endl;
        break;
      case 28:
        cout << "SEnergyCorrelator::GetJetPtBin(double) getting jet pT bin..." << endl;
        break;
      case 29:
        cout << "SEnergyCorrelator::CloseInputFile() closing input file..." << endl;
        break;
      case 30:
        cout << "SEnergyCorrelator::CloseOutputFile() closing output file..." << endl;
        break;
      case 31:
        cout << "SEnergyCorrelator::DoCorrelatorCalculation() looping over events and calculating correlators..." << endl;
        break;
      case 32:
        cout << "SEnergyCorrelator::SetSubEventsToUse(uint16_t, vector<int>) setting sub-events to use..." << endl;
        break;
      case 33:
        cout << "SEnergyCorrelator::CheckIfSubEvtGood(int) checking if sub-event is good..." << endl;
        break;
    }
    return;

  }  // end 'PrintDebug(uint32_t)'



  void SEnergyCorrelator::PrintError(const uint32_t code, const size_t nDrBinEdges, const size_t iDrBin, const string sInFileName) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) PrintDebug(23);

    switch (code) {
      case 0:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::Init(PHCompositeNode*) PANIC: calling complex method in standalone mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling complex method in standalone mode! Aborting!" << endl;
        }
        break;
      case 1:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::GrabInputNode() PANIC: couldn't grab node \"" << m_config.inNodeName << "\"! Aborting!" << endl;
        } else {
          cerr << "PANIC: couldn't grab node \"" << m_config.inNodeName << "\"! Aborting!" << endl;
        }
        break;
      case 2:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::GrabInputNode() PANIC: couldn't grab tree \"" << m_config.inTreeName << "\" from node \"" << m_config.inNodeName << "\"! Aborting!" << endl;
        } else {
          cerr << "PANIC: couldn't grab tree \"" << m_config.inTreeName << "\" from node \"" << m_config.inNodeName << "\"! Aborting!" << endl;
        }
        break;
      case 3:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::process_event(PHCompositeNode*) PANIC: calling complex method in standalone mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling complex method in standalone mode! Aborting!" << endl;
        }
        break;
      case 4:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::End(PHCompositeNode*) PANIC: calling complex method in standalone mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling complex method in standalone mode! Aborting!" << endl;
        }
        break;
      case 5:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::Init() PANIC: calling standalone method in complex mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling standalone method in complex mode! Aborting!" << endl;
        }
        break;
      case 6:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::OpenInputFiles() PANIC: couldn't create input TChain! Aborting" << endl;
        } else {
          cerr << "PANIC: couldn't create input TChain! Aborting!" << endl;
        }
        break;
      case 7:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::OpenInputFiles() PANIC: couldn't grab tree \"" << m_config.inTreeName << "\" from file \"" << sInFileName << "\"! Aborting!" << endl;
        } else {
          cerr << "PANIC: couldn't grab tree \"" << m_config.inTreeName << "\" from file \"" << sInFileName << "\"! Aborting!" << endl;
        }
        break;
      case 8:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::Analyze() PANIC: calling standalone method in complex mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling standalone method in complex mode! Aborting!" << endl;
        }
        break;
      case 9:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::End() PANIC: calling standalone method in complex mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling standalone method in complex mode! Aborting!" << endl;
        }
        break;
      case 10:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::InitializeTree() PANIC: no TTree! Aborting!" << endl;
        } else {
          cerr << "PANIC: no TTree! Aborting!" << endl;
        }
        break;
      case 11:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::OpenOutputFile() PANIC: couldn't open output file! Aborting!" << endl;
        } else {
          cerr << "PANIC: couldn't open output file! Aborting!" << endl;
        }
        break;
      case 12:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelator::ExtraHistsFromCorr() PANIC: number of dR bin edges is no good! Aborting!" << endl;
        } else {
          cerr << "PANIC: number of dR bin edges is no good! Aborting!\n"
               << "       nDrBinEdges = " << nDrBinEdges << ", nDrBins = " << m_config.nBinsDr
               << endl;
        }
        break;
      case 13:
        if (m_config.isStandalone) {
          cerr << "WARNING: dR bin #" << iDrBin << " with variance has a NAN as content or error..." << endl;
        }
        break;
      case 14:
        if (m_config.isStandalone) {
          cerr << "WARNING: dR bin #" << iDrBin << " with statistical error has a NAN as content or error..." << endl;
        }
        break;
      case 15:
        if (!m_config.isStandalone) {
          cerr << "SEnergyCorrelatorFile::End() PANIC: calling standalone method in complex mode! Aborting!" << endl;
        } else {
          cerr << "PANIC: calling standalone method in complex mode! Aborting!" << endl;
        }
        break;
    }
    return;

  }  // end 'PrintError(unint32_t)'



  bool SEnergyCorrelator::CheckCriticalParameters() {

    // print debugging statement
    if (m_config.isDebugOn) PrintDebug(21); 

    /* TODO checking goes here */
    return true;

  }  // end 'CheckCriticalParameters()'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
