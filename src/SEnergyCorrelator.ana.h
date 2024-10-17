/// ---------------------------------------------------------------------------
/*! \file    SEnergyCorrelator.ana.h
 *  \authors Derek Anderson, Alex Clarke
 *  \date    01.20.2023
 *
 *  A module to run ENC calculations in the sPHENIX
 *  software stack for the Cold QCD EEC analysis.
 */
/// ---------------------------------------------------------------------------

#pragma once

// make common namespaces implicit
using namespace std;
using namespace fastjet;
using namespace ROOT::Math;



namespace SColdQcdCorrelatorAnalysis {

  // analysis methods =========================================================

  // --------------------------------------------------------------------------
  //! Run local (in-jet) ENC calculations
  // --------------------------------------------------------------------------
  /*! FIXME move smearing to a seperate function */
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

      // create jet 4-vector, smear if need be
      PtEtaPhiEVector pJetVector(
        m_input.jets[iJet].GetPT(),
        m_input.jets[iJet].GetEta(),
        m_input.jets[iJet].GetPhi(),
        m_input.jets[iJet].GetEne()
      );
      if (m_config.modJets) SmearJetMomentum(pJetVector);

      // constituent loop
      for (uint64_t iCst = 0; iCst < m_input.csts[iJet].size(); iCst++) {

        // if needed, apply constituent cuts
        const bool isGoodCst = IsGoodCst( m_input.csts[iJet][iCst] );
        if (!isGoodCst) continue;

        // create cst 4-vector, smear if need be
        PtEtaPhiMVector pVecCst(
          m_input.csts[iJet][iCst].GetPT(),
          m_input.csts[iJet][iCst].GetEta(),
          m_input.csts[iJet][iCst].GetPhi(),
          Const::MassPion()
        );
        if (m_config.modCsts) SmearCstMomentum(pVecCst);

        // apply efficiency if need be
        if (m_config.doCstEff) {
          bool surives = SurvivesEfficiency(pVecCst.Pt());
          if (!survives) continue;
        }

	//Components to use for pseudoJet
	float Px = pVecCst.Px();
	float Py = pVecCst.Py();
	float Pz = pVecCst.Pz();
	float cstE = pVecCst.E();
	float skipCst = false;

	//Create pseudoJet for original cst info
	PseudoJet CstPseudo(Px, Py, Pz, cstE);

	//Apply smearing if necessary
	if(m_config.modCsts){
	}

        // create pseudojet
        PseudoJet constituent(Px, Py, Pz, cstE);
        constituent.set_user_index(iCst);

        // add to list
        m_jetCstVector.push_back(constituent);

      }  // end cst loop

      // run eec computation(s)
      if (m_config.doPackageCalc) {
        DoLocalCalcWithPackage( jet_pT  );
      }
      if(m_config.doManualCalc) {
	DoLocalCalcManual(m_jetCstVector, pJetVector);
      }
    }  // end jet loop
    return;

  }  // end 'DoLocalCalculation()'



  // --------------------------------------------------------------------------
  //! Run local (in-jet) calculations using P. T. Komiske's EEC package
  // --------------------------------------------------------------------------
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



  // --------------------------------------------------------------------------
  //! Run local (in-jet) calculations manually
  // --------------------------------------------------------------------------
  void SEnergyCorrelator::DoLocalCalcManual(
    const vector<fastjet::PseudoJet> momentum,
    PtEtaPhiEVector normalization
  ){

    //Get norm
    double norm = GetWeight(normalization, m_config.normOption);
    //Loop over csts
    for(uint64_t iCst = 0; iCst < momentum.size(); iCst++){
      //Get weightA
      PtEtaPhiEVector cstVec(
	  momentum[iCst].pt(),
          momentum[iCst].eta(),
	  momentum[iCst].phi(),
          pow(pow(momentum[iCst].pt(), 2) + pow(Const::MassPion(), 2), 0.5)
        );
      double weightA = GetWeight(cstVec, m_config.momOption, normalization);
      //Start second cst loop
      for(uint64_t jCst = 0; jCst < momentum.size(); jCst++){
        //Get weightB
	PtEtaPhiEVector cstVecB(
	  momentum[jCst].pt(),
          momentum[jCst].eta(),
	  momentum[jCst].phi(),
          pow(pow(momentum[jCst].pt(), 2) + pow(Const::MassPion(), 2), 0.5)
        );
	double weightB = GetWeight(cstVecB, m_config.momOption, normalization);
	const double dhCstAB = (momentum[iCst].rap()-momentum[jCst].rap());
	double dfCstAB = std::fabs(momentum[iCst].phi()-momentum[jCst].phi());
	if(dfCstAB > TMath::Pi()) dfCstAB = 2*TMath::Pi() - dfCstAB;
	const double drCstAB  = sqrt((dhCstAB * dhCstAB) + (dfCstAB * dfCstAB));
	const double eecWeight = (weightA*weightB)/(norm*norm);

	//Fill manual eecs
	for(size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++){
	  bool isInPtBin = ((normalization.Pt() >= m_config.ptJetBins[iPtBin].first) && (normalization.Pt() < m_config.ptJetBins[iPtBin].second));
	  if(isInPtBin){
	    m_outManualHistErrDrAxis[iPtBin]->Fill(drCstAB, eecWeight);
	  }
	}//end of pT bin loop
	//Start of third cst Loop
	for(uint64_t kCst = 0; kCst < momentum.size() && m_config.doThreePoint; kCst++){
	  PtEtaPhiEVector cstVecC(
	    momentum[kCst].pt(),
            momentum[kCst].eta(),
	    momentum[kCst].phi(),
            pow(pow(momentum[kCst].pt(), 2) + pow(Const::MassPion(), 2), 0.5)
          );
	  double weightC = GetWeight(cstVecC, m_config.momOption, normalization);
	  const double dhCstAC = (momentum[iCst].rap()-momentum[kCst].rap());
	  double dfCstAC = std::fabs(momentum[iCst].phi()-momentum[kCst].phi());
	  if(dfCstAC > TMath::Pi()) dfCstAC = 2*TMath::Pi() - dfCstAC;
	  const double drCstAC  = sqrt((dhCstAC * dhCstAC) + (dfCstAC * dfCstAC));
	  const double dhCstBC = (momentum[jCst].rap()-momentum[kCst].rap());
	  double dfCstBC = std::fabs(momentum[jCst].phi()-momentum[kCst].phi());
	  if(dfCstBC > TMath::Pi()) dfCstBC = 2*TMath::Pi() - dfCstBC;
	  const double drCstBC  = sqrt((dhCstBC * dhCstBC) + (dfCstBC * dfCstBC));
	  const double e3cWeight = (weightA*weightB*weightC)/(norm*norm*norm);

	  //Determine RL and RS
	  double RL = std::max(std::max(drCstAB, drCstAC), drCstBC);
	  double RS = std::min(std::min(drCstAB, drCstAC), drCstBC);
	  
	  //Fill Projected E3C
	  for(size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++){
	    bool isInPtBin = ((normalization.Pt() >= m_config.ptJetBins[iPtBin].first) && (normalization.Pt() < m_config.ptJetBins[iPtBin].second));
	    if(isInPtBin){
	      m_outProjE3C[iPtBin]->Fill(RL, e3cWeight);
	    }
	  }//end of ptBin loop

	  //Get RM
	  double RM = drCstAB; //set RM default value
	  if((drCstAB >= drCstAC && drCstAB <= drCstBC) || (drCstAB <= drCstAC && drCstAB >= drCstBC)) RM = drCstAB;
	  if((drCstAC >= drCstAB && drCstAC <= drCstBC) || (drCstAC <= drCstAB && drCstAC >= drCstBC)) RM = drCstAC;
	  if((drCstBC >= drCstAB && drCstBC <= drCstAC) || (drCstBC <= drCstAB && drCstBC >= drCstAC)) RM = drCstBC;

	  //skip in case RS or RM are 0
	  if(RS == 0 || RM == 0) continue;

	  //Get Parameterization
	  const double xi = RS/RM;
	  const double phi = std::asin(sqrt(1 - pow(RL-RM, 2)/(RS*RS)));

	  //Fill E3Cs
	  for(size_t iPtBin = 0; iPtBin < m_config.ptJetBins.size(); iPtBin++){
	    bool isInPtBin = ((normalization.Pt() >= m_config.ptJetBins[iPtBin].first) && (normalization.Pt() < m_config.ptJetBins[iPtBin].second));
	    if(isInPtBin){
	      for(size_t jRLBin = 0; jRLBin < m_config.rlBins.size(); jRLBin++){
		if(RL >= m_config.rlBins[jRLBin].first && RL < m_config.rlBins[jRLBin].second){
		  m_outE3C[iPtBin][jRLBin]->Fill(xi, phi, e3cWeight);
		}
	      }
	    }
	  }//end of ptBin loop
	}//end of third cst loop
      }//end of second cst loop
    }//end of first cst loop

  }  // end 'DoLocalCalcManual(vector<fastjet::PseudoJet>, PtEtaPhiEVector)'



  // --------------------------------------------------------------------------
  //! Extract histograms from ENC package
  // --------------------------------------------------------------------------
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

      TString sVarDrAxisName("hPackageCorrelatorVarianceDrAxis");
      TString sErrDrAxisName("hPackageCorrelatorErrorDrAxis");
      TString sVarLnDrAxisName("hPackageCorrelatorVarianceLnDrAxis");
      TString sErrLnDrAxisName("hPackageCorrelatorErrorLnDrAxis");

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
      m_outPackageHistVarDrAxis[iPtBin]   = new TH1D(sVarDrAxisName.Data(),   "", m_config.nBinsDr, drBinEdgeArray);
      m_outPackageHistErrDrAxis[iPtBin]   = new TH1D(sErrDrAxisName.Data(),   "", m_config.nBinsDr, drBinEdgeArray);
      m_outPackageHistVarLnDrAxis[iPtBin] = new TH1D(sVarLnDrAxisName.Data(), "", m_config.nBinsDr, lnDrBinEdgeArray);
      m_outPackageHistErrLnDrAxis[iPtBin] = new TH1D(sErrLnDrAxisName.Data(), "", m_config.nBinsDr, lnDrBinEdgeArray);
      m_outPackageHistVarDrAxis[iPtBin]   -> Sumw2();
      m_outPackageHistErrDrAxis[iPtBin]   -> Sumw2();
      m_outPackageHistVarLnDrAxis[iPtBin] -> Sumw2();
      m_outPackageHistErrLnDrAxis[iPtBin] -> Sumw2();

      // set bin content
      for (size_t iDrEdge = 0; iDrEdge < m_config.nBinsDr; iDrEdge++) {
        const size_t iDrBin        = iDrEdge;
        const double binVarContent = histContentAndVariance.first.at(iDrEdge);
        const double binVarError   = histContentAndVariance.second.at(iDrEdge);
        const double binErrContent = histContentAndError.first.at(iDrEdge);
        const double binErrError   = histContentAndError.second.at(iDrEdge);

        // check if bin with variances is good & set content/error
        const bool areVarBinValuesNans = (isnan(binVarContent) || isnan(binVarError));
        if (areVarBinValuesNans) {
          PrintError(13, 0, iDrBin);
        } else {
          m_outPackageHistVarDrAxis[iPtBin]   -> SetBinContent(iDrBin, binVarContent);
          m_outPackageHistVarLnDrAxis[iPtBin] -> SetBinContent(iDrBin, binVarContent);
          m_outPackageHistVarDrAxis[iPtBin]   -> SetBinError(iDrBin, binVarError);
          m_outPackageHistVarLnDrAxis[iPtBin] -> SetBinError(iDrBin, binVarError);
        }

        // check if bin with errors is good & set content/error
        const bool areErrBinValuesNans = (isnan(binErrContent) || isnan(binErrError));
        if (areErrBinValuesNans) {
          PrintError(14, 0, iDrBin);
        } else {
          m_outPackageHistErrDrAxis[iPtBin]   -> SetBinContent(iDrBin, binErrContent);
          m_outPackageHistErrLnDrAxis[iPtBin] -> SetBinContent(iDrBin, binErrContent);
          m_outPackageHistErrDrAxis[iPtBin]   -> SetBinError(iDrBin, binErrError);
          m_outPackageHistErrLnDrAxis[iPtBin] -> SetBinError(iDrBin, binErrError);
        }
      }  // end dr bin edge loop
    }  // end jet pt bin loop

    PrintMessage(14);
    return;

  }  // end 'ExtractHistsFromCorr()'



  // --------------------------------------------------------------------------
  //! Check if a jet satisfies cuts
  // --------------------------------------------------------------------------
  bool SEnergyCorrelator::IsGoodJet(const Types::JetInfo& jet) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(26);

    const bool isGoodJet = jet.IsInAcceptance(m_config.jetAccept);
    return isGoodJet;

  }  // end 'IsGoodJet(Types::JetInfo&)'



  // --------------------------------------------------------------------------
  //! Check if a constituent satisfies cuts
  // --------------------------------------------------------------------------
  bool SEnergyCorrelator::IsGoodCst(const Types::CstInfo& cst) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(27);

    // if analyzing truth tree and needed, check embedding ID
    bool isSubEvtGood = true;
    if (m_config.isInTreeTruth && m_config.selectSubEvts) {
      isSubEvtGood = Tools::IsSubEvtGood(
        cst.GetEmbedID(),
        m_config.subEvtOpt,
        m_config.isEmbed
      );
    }

    // if needed, apply cst. cuts
    bool isInCstAccept = true;
    if (m_config.applyCstCuts) {
      isInCstAccept = cst.IsInAcceptance(m_config.cstAccept);
    }

    const bool isGoodCst = (isSubEvtGood && isInCstAccept);
    return isGoodCst;

  }  // end 'IsGoodCst(Types::CstInfo&)'



  // --------------------------------------------------------------------------
  //! Check if a value (e.g. cst pt) passes efficiency
  // --------------------------------------------------------------------------
  bool SEnergyCorrelator::SurvivesEfficiency(const double value) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(37);

    // apply efficiency if need be
    const float rand = m_rando          -> Uniform(0., 1.);
    const float eff  = m_config.funcEff -> Eval(value);

    // return if value passed
    return (rand <= eff);

  }  // end 'SurvivesEfficiency(double)'



  // --------------------------------------------------------------------------
  //! Get weight of a point (e.g. a constituent) wrt. a reference (e.g. a jet)
  // --------------------------------------------------------------------------
  double SEnergyCorrelator::GetWeight(
    PtEtaPhiEVector momentum,
    int option,
    optional<PtEtaPhiEVector> reference
  ){

    double weight = 1;
    double Et = 1;
    if(reference.has_value()){
      TVector3 pRef(reference.value().X(), reference.value().Y(), reference.value().Z());
      TVector3 pMom(momentum.X(), momentum.Y(), momentum.Z());
      TVector3 pEt = pMom - (pMom*pRef/(pRef*pRef))*pRef;
      Et = pow(pEt.Mag2()+pow(Const::MassPion(),2),0.5);
    }
    switch(option){
      case Norm::Et:
	weight = Et;
	break;
      case Norm::E:
	weight = momentum.E();
	break;
      case Norm::Pt:
	[[fallthrough]];
      default:
	weight = momentum.Pt();
	break;
    }
    return weight;

  }  // end 'GetWeight(PtEtaPhiEVector, int, optional<PtEtaPhiEVector>)'



  // --------------------------------------------------------------------------
  //! Get bin no. for a given jet pt
  // --------------------------------------------------------------------------
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



  // --------------------------------------------------------------------------
  //! Smear jet momentum
  // --------------------------------------------------------------------------
  void SEnergyCorrelator::SmearJetMomentum(PtEtaPhiEVector& pJet) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(35);

    // grab unsmeared jet pt, E & calculate mass
    const double ptOrig = pJet.Pt();
    const double eOrig  = pJet.E();
    const double mJet   = pJet.M();

    // apply smearing
    const double ptSmear = ptOrig + m_rando -> Gaus(0., m_config.ptJetSmear * ptOrig);
    const double eSmear  = sqrt( hypot(ptSmear, mJet) ); 

    // update 4-vector and exit
    pJet.SetPt(ptSmear);
    pJet.SetE(eSmear);
    return;

  }  // end 'SmearJetMomentum(PtEtaPhiEVector&)'




  // --------------------------------------------------------------------------
  //! Smear constituent momentum
  // --------------------------------------------------------------------------
  void SEnergyCorrelator::SmearCstMomentum(PtEtaPhiEVector& pCst) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 7)) PrintDebug(36);


    // grab unsmeared cst pt
    const double ptOrig = pCst.Pt();

    // apply smearings
    //   - FIXME should I add flags for each type of modification?
    //     e.g. pt smearing and angular smearing?
    //   - Or could I do something w/ lambdas...
    const double ptSmear ptOrig + m_rando -> Gaus(0., m_config.ptCstSmear * ptOrig);

    /* REFRENCE
	  //Apply pT smearing if needed
	  float newpT = m_input.csts[iJet][iCst].GetPT();
	  if(m_config.ptCstSmear != 0) newpT+=shift->Gaus(0, m_config.ptCstSmear*m_input.csts[iJet][iCst].GetPT());
	  PtEtaPhiMVector rVecCstCopy(newpT, m_input.csts[iJet][iCst].GetEta(), m_input.csts[iJet][iCst].GetPhi(), Const::MassPion());
	  TVector3 p(rVecCstCopy.X(), rVecCstCopy.Y(), rVecCstCopy.Z());
	  TVector3 rotation_axis(rVecCstCopy.X(), rVecCstCopy.Y(), rVecCstCopy.Z());
	  //Apply angular smearing if needed
	  if(m_config.theta != 0){
	    float Theta0 = p.Theta();
	    float deltaTheta = shift->Gaus(0, m_config.theta);
	    float deltaPhi = shift->Uniform(-TMath::Pi(), TMath::Pi());
	    p.SetTheta(Theta0+deltaTheta);
	    p.Rotate(deltaPhi, rotation_axis);
	  }
	  //Change values to use in pseudoJet
	  Px = p.X();
	  Py = p.Y();
	  Pz = p.Z();
	  cstE = rVecCstCopy.E(); 
    */

    // update 4-vector and exit
    pCst.SetPt(ptSmear);
    return;

  }  // end 'SmearCstMomentum(PtEtaPhiEVector&)'

}  // end SColdQcdCorrelatorAnalysis namespace

// end ------------------------------------------------------------------------
