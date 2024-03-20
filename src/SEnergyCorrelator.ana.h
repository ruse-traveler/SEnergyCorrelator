// ----------------------------------------------------------------------------
// 'SEnergyCorrelator.ana.h'
// Derek Anderson
// 02.14.2023
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

  // analysis methods ---------------------------------------------------------

  void SEnergyCorrelator::DoLocalCalculation() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(31);

    // jet loop
    for (uint64_t iJet = 0; iJet < m_input.jets.size(); iJet++) {

      // clear vector for correlator
      m_jetCstVector.clear();

      // select jet pt bin & apply jet cuts
      const bool isGoodJet = IsGoodJet( m_input.jets[iJet] );
      if (!isGoodJet) continue;

      // constituent loop
      for (uint64_t iCst = 0; iCst < m_input.csts[iJet].size(); iCst++) {

        // create cst 4-vector
        ROOT::Math::PtEtaPhiMVector pVecCst(
          m_input.csts[iJet][iCst].GetPT(),
          m_input.csts[iJet][iCst].GetEta(),
          m_input.csts[iJet][iCst].GetPhi(),
          Const::MassPion()
        );

        // if needed, apply constituent cuts
        const bool isGoodCst = IsGoodCst( m_input.csts[iJet][iCst] );
        if (!isGoodCst) continue;

        // create pseudojet
        PseudoJet constituent(
          pVecCst.Px(),
          pVecCst.Py(),
          pVecCst.Pz(),
          pVecCst.E()
        );
        constituent.set_user_index(iCst);

        // add to list
        m_jetCstVector.push_back(constituent);

      }  // end cst loop

      // run eec computation(s)
      if (m_config.doPackageCalc) {
        DoLocalCalcWithPackage( m_input.jets[iJet].GetPT()  );
      }
    }  // end jet loop
    return;

  }  // end 'DoLocalCalculation()'



  void SEnergyCorrelator::DoLocalCalcWithPackage(const double ptJet) {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(31);

    const uint32_t iPtJetBin = GetJetPtBin( ptJet );
    const bool     foundBin  = (iPtJetBin >= 0);
    const bool     hasCsts   = (m_jetCstVector.size() > 0); 
    if (foundBin && hasCsts) {
      m_eecLongSide[iPtJetBin] -> compute(m_jetCstVector);
    }
    return;

  }  // end 'DoLocalCalcWithPackage(double)'



  void SEnergyCorrelator::ExtractHistsFromCorr() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(25);

    vector<double>                       drBinEdges;
    vector<double>                       lnDrBinEdges;
    pair<vector<double>, vector<double>> histContentAndVariance;
    pair<vector<double>, vector<double>> histContentAndError;
    for (size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {

      // create names
      TString sPtJetBin("_ptJet");
      sPtJetBin += floor(m_config.ptJetBins[iPtBin].first);

      TString sVarDrAxisName("hCorrelatorVarianceDrAxis");
      TString sErrDrAxisName("hCorrelatorErrorDrAxis");
      TString sVarLnDrAxisName("hCorrelatorVarianceLnDrAxis");
      TString sErrLnDrAxisName("hCorrelatorErrorLnDrAxis");
      sVarDrAxisName.Append(sPtJetBin.Data());
      sErrDrAxisName.Append(sPtJetBin.Data());
      sVarLnDrAxisName.Append(sPtJetBin.Data());
      sErrLnDrAxisName.Append(sPtJetBin.Data());

      // clear vectors
      drBinEdges.clear();
      lnDrBinEdges.clear();
      histContentAndVariance.first.clear();
      histContentAndVariance.second.clear();
      histContentAndError.first.clear();
      histContentAndError.second.clear();

      // grab bin edges, content, and error
      drBinEdges             = m_eecLongSide[iPtBin] -> bin_edges();
      histContentAndVariance = m_eecLongSide[iPtBin] -> get_hist_vars();
      histContentAndError    = m_eecLongSide[iPtBin] -> get_hist_errs();

      // create ln(dr) bin edges and arrays
      const size_t nDrBinEdges = drBinEdges.size();
      for (size_t iDrEdge = 0; iDrEdge < nDrBinEdges; iDrEdge++) {
        const double drEdge   = drBinEdges.at(iDrEdge);
        const double lnDrEdge = log(drEdge);
        lnDrBinEdges.push_back(lnDrEdge);
      }

      // bin edge arrays for histograms
      double drBinEdgeArray[nDrBinEdges];
      double lnDrBinEdgeArray[nDrBinEdges];
      for (size_t iDrEdge = 0; iDrEdge < nDrBinEdges; iDrEdge++) {
        drBinEdgeArray[iDrEdge]   = drBinEdges.at(iDrEdge);
        lnDrBinEdgeArray[iDrEdge] = lnDrBinEdges.at(iDrEdge);
      }

      // make sure number of bin edges is reasonable
      const bool isNumBinEdgesGood = ((nDrBinEdges - 1) == m_config.nBinsDr);
      if (!isNumBinEdgesGood) {
        PrintError(12, nDrBinEdges);
        assert(isNumBinEdgesGood);
      }

      // create output histograms
      m_outHistVarDrAxis[iPtBin]   = new TH1D(sVarDrAxisName.Data(),   "", m_config.nBinsDr, drBinEdgeArray);
      m_outHistErrDrAxis[iPtBin]   = new TH1D(sErrDrAxisName.Data(),   "", m_config.nBinsDr, drBinEdgeArray);
      m_outHistVarLnDrAxis[iPtBin] = new TH1D(sVarLnDrAxisName.Data(), "", m_config.nBinsDr, lnDrBinEdgeArray);
      m_outHistErrLnDrAxis[iPtBin] = new TH1D(sErrLnDrAxisName.Data(), "", m_config.nBinsDr, lnDrBinEdgeArray);
      m_outHistVarDrAxis[iPtBin]   -> Sumw2();
      m_outHistErrDrAxis[iPtBin]   -> Sumw2();
      m_outHistVarLnDrAxis[iPtBin] -> Sumw2();
      m_outHistErrLnDrAxis[iPtBin] -> Sumw2();

      // set bin content
      for (size_t iDrEdge = 0; iDrEdge < m_config.nBinsDr; iDrEdge++) {
        const size_t iDrBin        = iDrEdge + 1;
        const double binVarContent = histContentAndVariance.first.at(iDrEdge);
        const double binVarError   = histContentAndVariance.second.at(iDrEdge);
        const double binErrContent = histContentAndError.first.at(iDrEdge);
        const double binErrError   = histContentAndError.second.at(iDrEdge);

        // check if bin with variances is good & set content/error
        const bool areVarBinValuesNans = (isnan(binVarContent) || isnan(binVarError));
        if (areVarBinValuesNans) {
          PrintError(13, 0, iDrBin);
        } else {
          m_outHistVarDrAxis[iPtBin]   -> SetBinContent(iDrBin, binVarContent);
          m_outHistVarLnDrAxis[iPtBin] -> SetBinContent(iDrBin, binVarContent);
          m_outHistVarDrAxis[iPtBin]   -> SetBinError(iDrBin, binVarError);
          m_outHistVarLnDrAxis[iPtBin] -> SetBinError(iDrBin, binVarError);
        }

        // check if bin with errors is good & set content/error
        const bool areErrBinValuesNans = (isnan(binErrContent) || isnan(binErrError));
        if (areErrBinValuesNans) {
          PrintError(14, 0, iDrBin);
        } else {
          m_outHistErrDrAxis[iPtBin]   -> SetBinContent(iDrBin, binErrContent);
          m_outHistErrLnDrAxis[iPtBin] -> SetBinContent(iDrBin, binErrContent);
          m_outHistErrDrAxis[iPtBin]   -> SetBinError(iDrBin, binErrError);
          m_outHistErrLnDrAxis[iPtBin] -> SetBinError(iDrBin, binErrError);
        }
      }  // end dr bin edge loop
    }  // end jet pt bin loop

    if (m_config.isStandalone) PrintMessage(14);
    return;

  }  // end 'ExtractHistsFromCorr()'



  bool SEnergyCorrelator::IsGoodJet(const Types::JetInfo& jet) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(26);

    const bool isGoodJet = jet.IsInAcceptance(m_config.jetAccept);
    return isGoodJet;

  }  // end 'IsGoodJet(double, double)'



  bool SEnergyCorrelator::IsGoodCst(const Types::CstInfo& cst) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(27);

    // if analyzing truth tree and needed, check embedding ID
    bool isSubEvtGood = true;
    if (m_config.isInTreeTruth && m_config.selectSubEvts) {
      isSubEvtGood = Tools::IsSubEvtGood( cst.GetEmbedID(), m_config.subEvtOpt, m_config.isEmbed );
    }

    // if needed, apply cst. cuts
    bool isInCstAccept = true;
    if (m_config.applyCstCuts) {
      isInCstAccept = cst.IsInAcceptance(m_config.cstAccept);
    }

    const bool isGoodCst = (isSubEvtGood && isInCstAccept);
    return isGoodCst;

  }  // end 'IsGoodCst(double, double)'



  int32_t SEnergyCorrelator::GetJetPtBin(const double ptJet) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(28);
 
    bool    isInPtBin(false);
    int32_t iJetPtBin(-1);
    for (size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {
      isInPtBin = ((ptJet >= m_config.ptJetBins[iPtBin].first) && (ptJet < m_config.ptJetBins[iPtBin].second));
      if (isInPtBin) {
        iJetPtBin = iPtBin;
        break; 
      }
    }
    return iJetPtBin;

  }  // end 'GetJetPtBin(double)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
