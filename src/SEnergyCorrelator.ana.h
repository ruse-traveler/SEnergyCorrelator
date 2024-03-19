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

  void SEnergyCorrelator::DoCorrelatorCalculation() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(31);

    // for eec calculation
    vector<PseudoJet> jetCstVector;

    // jet loop
    for (uint64_t iJet = 0; iJet < m_input.jets.size(); iJet++) {

      // clear vector for correlator
      jetCstVector.clear();

      // select jet pt bin & apply jet cuts
      const uint32_t iPtJetBin = GetJetPtBin( m_input.jets[iJet].GetPT() );
      const bool     isGoodJet = ApplyJetCuts( m_input.jets[iJet].GetPT(), m_input.jets[iJet].GetEta() );
      if (!isGoodJet) continue;

      // constituent loop
      for (uint64_t iCst = 0; iCst < m_input.csts[iJet].size(); iCst++) {

        // create cst 4-vector
        ROOT::Math::PtEtaPhiMVector rVecCst(
          m_input.csts[iJet][iCst].GetPT(),
          m_input.csts[iJet][iCst].GetEta(),
          m_input.csts[iJet][iCst].GetPhi(),
          Const::MassPion()
        );

        // if truth tree and needed, check embedding ID
        if (m_config.isInTreeTruth && m_config.selectSubEvts) {
          const bool isSubEvtGood = Tools::IsSubEvtGood( m_input.csts[iJet][iCst].GetEmbedID(), m_config.subEvtOpt, m_config.isEmbed );
          if (!isSubEvtGood) continue;
        }

        // if needed, apply constituent cuts
        const bool isGoodCst = ApplyCstCuts( m_input.csts[iJet][iCst].GetPT(), m_input.csts[iJet][iCst].GetEta() );
        if (m_config.applyCstCuts && !isGoodCst) continue;

        // create pseudojet
        PseudoJet constituent(
          rVecCst.Px(),
          rVecCst.Py(),
          rVecCst.Pz(),
          rVecCst.E()
        );
        constituent.set_user_index(iCst);

        // add to list
        jetCstVector.push_back(constituent);

      }  // end cst loop

      // run eec computation
      for (size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {
        const bool isInPtBin = ((m_input.jets[iJet].GetPT() >= m_config.ptJetBins[iPtBin].first) && (m_input.jets[iJet].GetPT() < m_config.ptJetBins[iPtBin].second));
        if (isInPtBin) {
          if (jetCstVector.size() > 0) {
            m_eecLongSide[iPtJetBin] -> compute(jetCstVector);
          }
        }
      }  // end pTjet bin loop
    }  // end jet loop
    return;

  }  // end 'DoCorrelatorCalculation()'



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



  bool SEnergyCorrelator::ApplyJetCuts(const double ptJet, const double etaJet) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(26);

    // grab info
    Types::JetInfo jet;
    jet.SetPT(ptJet);
    jet.SetEta(etaJet);

    const bool isInJetAccept = jet.IsInAcceptance(m_config.jetAccept);
    return isInJetAccept;

  }  // end 'ApplyJetCuts(double, double)'



  bool SEnergyCorrelator::ApplyCstCuts(const double momCst, const double drCst) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(27);

    // grab info
    Types::CstInfo cst;
    cst.SetDR(drCst);
    cst.SetEne(hypot(momCst, Const::MassPion()));

    const bool isInCstAccept = cst.IsInAcceptance(m_config.cstAccept);
    return isInCstAccept;

  }  // end 'ApplyCstCuts(double, double)'



  uint32_t SEnergyCorrelator::GetJetPtBin(const double ptJet) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(28);
 
    bool     isInPtBin(false);
    uint32_t iJetPtBin(0);
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
