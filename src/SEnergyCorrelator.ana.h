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

// c++ utilities
#include <cmath>
#include <cassert>
// root libraries
#include <TString.h>
#include <Math/Vector4D.h>
// fastjet libraries
#include <fastjet/PseudoJet.hh>
// anlysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/SCorrelatorUtilities.h"

// make common namespaces implicit
using namespace std;
using namespace fastjet;
using namespace SColdQcdCorrelatorAnalysis::SCorrelatorUtilities;



namespace SColdQcdCorrelatorAnalysis {

  // analysis methods ---------------------------------------------------------

  void SEnergyCorrelator::DoCorrelatorCalculation() {

    // print debug statement
    if (m_config.isDebugOn) PrintDebug(31);

    // announce start of event loop
    const uint64_t nEvts = m_inChain -> GetEntriesFast();
    PrintMessage(7, nEvts);

    // for eec calculation
    vector<PseudoJet> jetCstVector;

    // event loop
    uint64_t nBytes = 0;
    for (uint64_t iEvt = 0; iEvt < nEvts; iEvt++) {

      const uint64_t entry = LoadTree(m_inChain, iEvt, m_fCurrent);
      if (entry < 0) break;

      const uint64_t bytes = GetEntry(m_inChain, iEvt);
      if (bytes < 0) {
        break;
      } else {
        nBytes += bytes;
        PrintMessage(8, nEvts, iEvt);
      }

      // jet loop
      uint64_t nJets = (int) m_evtNumJets;
      for (uint64_t iJet = 0; iJet < nJets; iJet++) {

        // clear vector for correlator
        jetCstVector.clear();

        // get jet info
        const uint64_t nCsts  = m_jetNumCst -> at(iJet);
        const double   ptJet  = m_jetPt     -> at(iJet);
        const double   etaJet = m_jetEta    -> at(iJet);

        // select jet pt bin & apply jet cuts
        const uint32_t iPtJetBin = GetJetPtBin(ptJet);
        const bool     isGoodJet = ApplyJetCuts(ptJet, etaJet);
        if (!isGoodJet) continue;

        // constituent loop
        for (uint64_t iCst = 0; iCst < nCsts; iCst++) {

          // get cst info
          const double drCst  = (m_cstDr  -> at(iJet)).at(iCst);
          const double etaCst = (m_cstEta -> at(iJet)).at(iCst);
          const double phiCst = (m_cstPhi -> at(iJet)).at(iCst);
          const double ptCst  = (m_cstPt  -> at(iJet)).at(iCst);

          // create cst 4-vector
          ROOT::Math::PtEtaPhiMVector rVecCst(ptCst, etaCst, phiCst, MassPion);

          // if truth tree and needed, check embedding ID
          if (m_config.isInTreeTruth && m_config.selectSubEvts) {
            const int  embedCst     = (m_cstEmbedID -> at(iJet)).at(iCst);
            const bool isSubEvtGood = IsSubEvtGood(embedCst, m_config.subEvtOpt, m_config.isEmbed);
            if (!isSubEvtGood) continue;
          }

          // if needed, apply constituent cuts
          const bool isGoodCst = ApplyCstCuts(ptCst, drCst);
          if (m_config.applyCstCuts && !isGoodCst) continue;

          // create pseudojet & add to list
          PseudoJet constituent(rVecCst.Px(), rVecCst.Py(), rVecCst.Pz(), rVecCst.E());
          constituent.set_user_index(iCst);
          jetCstVector.push_back(constituent);
        }  // end cst loop

        // run eec computation
        for (size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++) {
          const bool isInPtBin = ((ptJet >= m_config.ptJetBins[iPtBin].first) && (ptJet < m_config.ptJetBins[iPtBin].second));
          if (isInPtBin) {
            if (jetCstVector.size() > 0) {
              m_eecLongSide[iPtJetBin] -> compute(jetCstVector);
            }
          }
        }  // end pTjet bin loop
      }  // end jet loop
    }  // end event loop

    // announce end of event loop
    PrintMessage(13);
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
    JetInfo jet;
    jet.pt  = ptJet;
    jet.eta = etaJet;

    const bool isInJetAccept = IsInAcceptance(jet, m_config.jetAccept.first, m_config.jetAccept.second);
    return isInJetAccept;

  }  // end 'ApplyJetCuts(double, double)'



  bool SEnergyCorrelator::ApplyCstCuts(const double momCst, const double drCst) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(27);

    // grab info
    CstInfo cst;
    cst.dr  = drCst;
    cst.ene = hypot(momCst, MassPion);

    const bool isInCstAccept = IsInAcceptance(cst, m_config.cstAccept.first, m_config.cstAccept.second);
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
