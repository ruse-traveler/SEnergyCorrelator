// 'SEnergyCorrelator.h'
// Derek Anderson
// 01.20.2023
//
// A module to implement Peter Komiske's
// EEC library in the sPHENIX software
// stack.

#ifndef SENERGYCORRELATOR_H
#define SENERGYCORRELATOR_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wpragmas"

// standard c includes
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
// f4a includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllHistoManager.h>
// phool includes
#include <phool/phool.h>
#include <phool/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
// root includes
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TDirectory.h>
// fastjet includes
#include <fastjet/PseudoJet.hh>
// eec include
#include "/sphenix/user/danderson/eec/EnergyEnergyCorrelators/eec/include/EEC.hh"

#pragma GCC diagnostic pop

using namespace std;

// forward declarations
class TH1;
class TFile;
class TTree;
class PHCompositeNode;



// SEnergyCorrelator definition -----------------------------------------------

class SEnergyCorrelator : public SubsysReco {

  public:

    // ctor/dtor
    SEnergyCorrelator(const string &name = "SEnergyCorrelator", const bool isComplex = false, const bool doDebug = false);
    ~SEnergyCorrelator() override;

    // F4A methods
    int Init(PHCompositeNode *topNode)          override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode)           override;

    // standalone-only methods
    void Init();
    void Analyze();
    void End();

    // setters
    void SetVerbosity(const int verb)           {m_verbosity   = verb;}
    void SetInputFile(const string &iFileName)  {m_inFileName  = iFileName;}
    void SetInputNode(const string &iNodeName)  {m_inNodeName  = iNodeName;}
    void SetInputTree(const string &iTreeName)  {m_inTreeName  = iTreeName;}
    void SetOutputFile(const string &oFileName) {m_outFileName = oFileName;}

    // system getters
    int    GetVerbosity()        {return m_verbosity;}
    bool   GetInDebugMode()      {return m_inDebugMode;}
    bool   GetInComplexMode()    {return m_inComplexMode;}
    bool   GetInStandaloneMode() {return m_inStandaloneMode;}
    string GetInputFileName()    {return m_inFileName;}
    string GetInputNodeName()    {return m_inNodeName;}
    string GetInputTreeName()    {return m_inTreeName;}
    string GetOutputFileName()   {return m_outFileName;}

  private:

    // io methods
    void GrabInputNode();
    void OpenInputFile();
    void OpenOutputFile();
    void SaveOutput();

    // system methods
    void InitializeMembers();
    void InitializeHists();
    void InitializeCorrs();
    void InitializeTree();
    void PrintMessage(const uint32_t code);
    void PrintError(const uint32_t code);

    // io members
    TFile *m_outFile;
    TFile *m_inFile;
    TTree *m_inTree;

    // system members
    int    m_fCurrent;
    int    m_verbosity;
    bool   m_inDebugMode;
    bool   m_inComplexMode;
    bool   m_inStandaloneMode;
    bool   m_isInputTreeTruth;
    string m_inFileName;
    string m_inNodeName;
    string m_inTreeName;
    string m_outFileName;

    // histogram members
    size_t         m_nJetPtBins;
    size_t         m_nOutHistBins;
    vector<double> m_histBinEdges;
    vector<double> m_outputHists;

    // input truth tree address members
    int    m_truParton3_ID;
    int    m_truParton4_ID;
    double m_truParton3_MomX;
    double m_truParton3_MomY;
    double m_truParton3_MomZ;
    double m_truParton4_MomX;
    double m_truParton4_MomY;
    double m_truParton4_MomZ;
    // input reco. tree address members
    int    m_recParton3_ID;
    int    m_recParton4_ID;
    double m_recParton3_MomX;
    double m_recParton3_MomY;
    double m_recParton3_MomZ;
    double m_recParton4_MomX;
    double m_recParton4_MomY;
    double m_recParton4_MomZ;

    // generic input tree address members
    int                     m_evtNumJets;
    vector<unsigned long>  *m_jetNumCst;
    vector<unsigned int>   *m_jetID;
    vector<unsigned int>   *m_jetTruthID;
    vector<double>         *m_jetEnergy;
    vector<double>         *m_jetPt;
    vector<double>         *m_jetEta;
    vector<double>         *m_jetPhi;
    vector<double>         *m_jetArea;
    vector<vector<double>> *m_cstZ;
    vector<vector<double>> *m_cstDr;
    vector<vector<double>> *m_cstEnergy;
    vector<vector<double>> *m_cstJt;
    vector<vector<double>> *m_cstEta;
    vector<vector<double>> *m_cstPhi;

    // input truth tree branch members
    TBranch *m_brTruParton3_ID;
    TBranch *m_brTruParton4_ID;
    TBranch *m_brTruParton3_MomX;
    TBranch *m_brTruParton3_MomY;
    TBranch *m_brTruParton3_MomZ;
    TBranch *m_brTruParton4_MomX;
    TBranch *m_brTruParton4_MomY;
    TBranch *m_brTruParton4_MomZ;
    // input reco. tree branch members
    TBranch *m_brRecParton3_ID;
    TBranch *m_brRecParton4_ID;
    TBranch *m_brRecParton3_MomX;
    TBranch *m_brRecParton3_MomY;
    TBranch *m_brRecParton3_MomZ;
    TBranch *m_brRecParton4_MomX;
    TBranch *m_brRecParton4_MomY;
    TBranch *m_brRecParton4_MomZ;

    // generic input tree branch members
    TBranch *m_brEvtNumJets;
    TBranch *m_brJetNumCst;
    TBranch *m_brJetID;
    TBranch *m_brJetTruthID;
    TBranch *m_brJetEnergy;
    TBranch *m_brJetPt;
    TBranch *m_brJetEta;
    TBranch *m_brJetPhi;
    TBranch *m_brJetArea;
    TBranch *m_brCstZ;
    TBranch *m_brCstDr;
    TBranch *m_brCstEnergy;
    TBranch *m_brCstJt;
    TBranch *m_brCstEta;
    TBranch *m_brCstPhi;

};

#endif

// end ------------------------------------------------------------------------
