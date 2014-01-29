#define eFakeSel_cxx

#include "eFakeSel.h"
#include <TH2.h>
#include <TStyle.h>
#include "TSystem.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
using namespace std;

string  selection       = "nunuGamma_GGF";//"datadriven_doubleMu";
string  period          = "2012";
string  suffix          = "SUFFIX";

bool    verbose         = false;
bool    mc              = false;
bool    qcd             = false;

reweight::LumiReWeighting LumiWeightsD_;
reweight::LumiReWeighting LumiWeightsD_sys_;

TF1* fjetreso = new TF1("jetReso","TMath::Sqrt(((sign([0])*(([0]/x)^2))+(([1]^2)*(x^([2]-1.))))+([3]^2))",0,15000);
TF1* fmex     = new TF1("Mex Sigma", "[0]*TMath::Sqrt(x) + [1]",0,5000);
TF1* fmey     = new TF1("Mey Sigma", "[0]*TMath::Sqrt(x) + [1]",0,5000);
TF1* fphoton  = new TF1("Photon_Sigma","[0]*pow(x,[1])",0,10000);
TF1* fmuon    = new TF1("ptsigma_barrel", "pol3", 0, 10000);

int n_vertex_count = 0;
int n_photon_count = 0;
int n_hltmatched_count = 0;
int pass_filter = 0;

void eFakeSel::Begin(TTree *tree)
{

   TString option = GetOption();

   //PU Reweighting:

   Float_t DataDist_2012D[60] = {1768.77, 3651.68, 7356.92, 14793.3, 51246.8, 489173, 2.63921e+06, 6.82361e+06, 1.88317e+07, 5.19794e+07, 1.08899e+08, 1.88257e+08, 2.57316e+08, 3.01258e+08, 3.27492e+08, 3.44354e+08, 3.59374e+08, 3.74823e+08, 3.90058e+08, 4.0217e+08, 4.08643e+08, 4.11422e+08, 4.11151e+08, 4.05573e+08, 3.92856e+08, 3.72226e+08, 3.44284e+08, 3.10504e+08, 2.7226e+08, 2.30759e+08, 1.8777e+08, 1.45857e+08, 1.07763e+08, 7.56393e+07, 5.0551e+07, 3.23928e+07, 2.0178e+07, 1.25011e+07, 7.95461e+06, 5.38133e+06, 3.95572e+06, 3.15182e+06, 2.66553e+06, 2.33493e+06, 2.08e+06, 1.86352e+06, 1.669e+06, 1.48952e+06, 1.32229e+06, 1.16642e+06, 1.02175e+06, 888282, 766092, 655189, 555473, 466702, 388487, 320302, 261506, 211367};

   Float_t DataDist_2012D_sys[60] = {1632.3, 3240.95, 6296.46, 12068, 30517.4, 227305, 1.49497e+06, 4.52239e+06, 1.0476e+07, 2.96274e+07, 6.82538e+07, 1.28218e+08, 2.00114e+08, 2.55292e+08, 2.90327e+08, 3.11739e+08, 3.26279e+08, 3.39668e+08, 3.53444e+08, 3.67108e+08, 3.78444e+08, 3.85013e+08, 3.88036e+08, 3.88742e+08, 3.85622e+08, 3.76942e+08, 3.61739e+08, 3.39978e+08, 3.12601e+08, 2.80828e+08, 2.45645e+08, 2.08081e+08, 1.69728e+08, 1.32723e+08, 9.92184e+07, 7.08564e+07, 4.84379e+07, 3.18819e+07, 2.04316e+07, 1.29836e+07, 8.39661e+06, 5.69316e+06, 4.14131e+06, 3.24848e+06, 2.71199e+06, 2.35979e+06, 2.10068e+06, 1.88918e+06, 1.70364e+06, 1.53423e+06, 1.37666e+06, 1.22917e+06, 1.0912e+06, 962650, 843537, 733913, 633795, 543116, 461705, 389280};
  
   Float_t MCDist_Summer2012_S10[60] = {2.560E-06,5.239E-06,1.420E-05,5.005E-05,1.001E-04,2.705E-04,1.999E-03,6.097E-03,1.046E-02,1.383E-02,1.685E-02,2.055E-02,2.572E-02,3.262E-02,4.121E-02,4.977E-02,5.539E-02,5.725E-02,5.607E-02,5.312E-02,5.008E-02,4.763E-02,4.558E-02,4.363E-02,4.159E-02,3.933E-02,3.681E-02,3.406E-02,3.116E-02,2.818E-02,2.519E-02,2.226E-02,1.946E-02,1.682E-02,1.437E-02,1.215E-02,1.016E-02,8.400E-03,6.873E-03,5.564E-03,4.457E-03,3.533E-03,2.772E-03,2.154E-03,1.656E-03,1.261E-03,9.513E-04,7.107E-04,5.259E-04,3.856E-04,2.801E-04,2.017E-04,1.439E-04,1.017E-04,7.126E-05,4.948E-05,3.405E-05,2.322E-05,1.570E-05,5.005E-06};

   std::vector<float> DataDistD;
   std::vector<float> DataDistD_sys;
   std::vector<float> MCDist;
  
   for( int i=0; i<60; ++i) {
     DataDistD.push_back(DataDist_2012D[i]);
     DataDistD_sys.push_back(DataDist_2012D_sys[i]);
     MCDist.push_back(MCDist_Summer2012_S10[i]);
   }

   LumiWeightsD_ = reweight::LumiReWeighting(MCDist, DataDistD);
   LumiWeightsD_sys_ = reweight::LumiReWeighting(MCDist, DataDistD_sys);

   // Get trigger names from jobTree                                                                                                                                                
   vector<string>* triggerNames = 0;
   TFile   *inFile         = tree->GetCurrentFile();
   TTree   *jobTree        = (TTree*)inFile->Get("ntupleProducer/jobTree");

   jobTree->SetBranchAddress("triggerNames", &triggerNames);
   jobTree->GetEntry();
   triggerSelector = new TriggerSelector(selection, period, *triggerNames);

   histoFile = new TFile("higgsHistograms_SUFFIX.root","RECREATE");
   histoFile->cd();
  
   //Setting up new tree
   histoFile->mkdir("Analyzer", "Analyzer");
   AnDir = (TDirectory*)histoFile->GetDirectory("Analyzer");


   if(verbose) cout << "CREATING TREE..." << endl;

   TZee = new TTree("TZee", "Tree with Z Candidates");
   TZep = new TTree("TZep", "Tree with Z Candidates");
   TZek = new TTree("TZek", "Tree with Z Candidates");
   TZpp = new TTree("TZpp", "Tree with Z Candidates");
   TZpk = new TTree("TZpk", "Tree with Z Candidates");
   TZkk = new TTree("TZkk", "Tree with Z Candidates");

   Zee = new TLorentzVector(-100, -100, 0, 0);
   Zep = new TLorentzVector(-100, -100, 0, 0);
   Zek = new TLorentzVector(-100, -100, 0, 0);
   Zpp = new TLorentzVector(-100, -100, 0, 0);
   Zpk = new TLorentzVector(-100, -100, 0, 0);
   Zkk = new TLorentzVector(-100, -100, 0, 0);

   Zee_E1 = new TLorentzVector(-100, -100, 0, 0);
   Zee_E2 = new TLorentzVector(-100, -100, 0, 0);
   Zep_E = new TLorentzVector(-100, -100, 0, 0);
   Zep_P = new TLorentzVector(-100, -100, 0, 0);
   Zek_E = new TLorentzVector(-100, -100, 0, 0);
   Zek_K = new TLorentzVector(-100, -100, 0, 0);
   Zpp_P1 = new TLorentzVector(-100, -100, 0, 0);
   Zpp_P2 = new TLorentzVector(-100, -100, 0, 0);
   Zpk_P = new TLorentzVector(-100, -100, 0, 0);
   Zpk_K = new TLorentzVector(-100, -100, 0, 0);
   Zkk_K1 = new TLorentzVector(-100, -100, 0, 0);
   Zkk_K2 = new TLorentzVector(-100, -100, 0, 0);

   pfMET = new TVector2(-10,-10);

   b_Zee          = TZee->Branch("Zee", "TLorentzVector", &Zee);
   b_Zep          = TZep->Branch("Zep", "TLorentzVector", &Zep);
   b_Zek          = TZek->Branch("Zek", "TLorentzVector", &Zek);
   b_Zpp          = TZpp->Branch("Zpp", "TLorentzVector", &Zpp);
   b_Zpk          = TZpk->Branch("Zpk", "TLorentzVector", &Zpk);
   b_Zkk          = TZkk->Branch("Zkk", "TLorentzVector", &Zkk);
   b_Zee_E1       = TZee->Branch("Zee_O1", "TLorentzVector", &Zee_E1);
   b_Zee_E2       = TZee->Branch("Zee_O2", "TLorentzVector", &Zee_E2);
   b_Zep_E        = TZep->Branch("Zep_O1", "TLorentzVector", &Zep_E );
   b_Zep_P        = TZep->Branch("Zep_O2", "TLorentzVector", &Zep_P );
   b_Zek_E        = TZek->Branch("Zek_O1", "TLorentzVector", &Zek_E );
   b_Zek_K        = TZek->Branch("Zek_O2", "TLorentzVector", &Zek_K );
   b_Zpp_P1       = TZpp->Branch("Zpp_O1", "TLorentzVector", &Zpp_P1);
   b_Zpp_P2       = TZpp->Branch("Zpp_O2", "TLorentzVector", &Zpp_P2);
   b_Zpk_P        = TZpk->Branch("Zpk_O1", "TLorentzVector", &Zpk_P );
   b_Zpk_K        = TZpk->Branch("Zpk_O2", "TLorentzVector", &Zpk_K );
   b_Zkk_K1       = TZkk->Branch("Zkk_O1", "TLorentzVector", &Zkk_K1);
   b_Zkk_K2       = TZkk->Branch("Zkk_O2", "TLorentzVector", &Zkk_K2);

   b_Zee_ntrks       = TZee->Branch("Zee_ntrks", &ntrks, "ntrks/F");
   b_Zep_ntrks       = TZep->Branch("Zep_ntrks", &ntrks, "ntrks/F");
   b_Zek_ntrks       = TZek->Branch("Zek_ntrks", &ntrks, "ntrks/F");
   b_Zpp_ntrks       = TZpp->Branch("Zpp_ntrks", &ntrks, "ntrks/F");
   b_Zpk_ntrks       = TZpk->Branch("Zpk_ntrks", &ntrks, "ntrks/F");
   b_Zkk_ntrks       = TZkk->Branch("Zkk_ntrks", &ntrks, "ntrks/F");
   b_Zee_MET         = TZee->Branch("Zee_MET", "TVector2", &pfMET);
   b_Zep_MET         = TZep->Branch("Zep_MET", "TVector2", &pfMET);
   b_Zek_MET         = TZek->Branch("Zek_MET", "TVector2", &pfMET);
   b_Zpp_MET         = TZpp->Branch("Zpp_MET", "TVector2", &pfMET);
   b_Zpk_MET         = TZpk->Branch("Zpk_MET", "TVector2", &pfMET);
   b_Zkk_MET         = TZkk->Branch("Zkk_MET", "TVector2", &pfMET);


if(verbose) cout << "INITIALIZING TREE VARIABLES..." << endl;

 // ch      nh       ph
  float EAPhoTemp[7][3] = {
    {0.012,  0.030,   0.148}, //         eta < 1.0  
    {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479   
    {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0  
    {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2 
    {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3  
    {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4 
    {0.012,  0.072,   0.266}  // 2.4   < eta       
  };

  for (unsigned int i =0; i<7; i++){
    for (unsigned int j =0; j<3; j++){
      EAPho[i][j] = EAPhoTemp[i][j];
    }
  }


}

void eFakeSel::SlaveBegin(TTree * /*tree*/)
{
   TString option = GetOption();
}

Bool_t eFakeSel::Process(Long64_t entry)
{
  pfMET = new TVector2(-10,-10);
  GetEntry(entry);

  if( entry % 1000 == 0 ) cout << "Processing event number: " << entry << endl;
  if(verbose) cout << "Processing event number: " << entry << endl;

  TLorentzVector Photon;

  if(mc)   fmex->SetParameters(0.53,3.08);
  else fmex->SetParameters(0.61,0.37);
  
  if(mc)   fmey->SetParameters(0.53,3.4); 
  else fmey->SetParameters(0.62,0.17);
 
  fphoton->SetParameters(0.1138,-0.449);
  fmuon->SetParameters(0.006821,0.0001047,-1.213e-08,1.874e-12);
 
  vector<TLorentzVector> Good_Photons;
  vector<TLorentzVector> Good_Electrons;
  vector<TLorentzVector> Good_eFakes;

  pfMET = new TVector2(-10,-10);
  TVector2 pfMet;
 
  
  if(!mc){
    //Trigger Status
    
    if (NoiseFilters_isScraping) { if(verbose) cout << "NoiseFilters_isScraping: " << NoiseFilters_isScraping << endl; return kTRUE;}
    if (NoiseFilters_isNoiseHcalHBHE) {if(verbose) cout << "NoiseFilters_isNoiseHcalHBHE: " << NoiseFilters_isNoiseHcalHBHE <<endl; return kTRUE;}
    if (NoiseFilters_isNoiseHcalLaser) {if(verbose) cout << "NoiseFilters_isNoiseHcalLaser: " <<NoiseFilters_isNoiseHcalLaser <<endl; return kTRUE;}
    if (NoiseFilters_isNoiseEcalTP) {if(verbose) cout << "NoiseFilters_isNoiseEcalTP: " << NoiseFilters_isNoiseEcalTP<< endl; return kTRUE;}
    if (NoiseFilters_isNoiseEcalBE) {if(verbose) cout<< "NoiseFilters_isNoiseEcalBE: " << NoiseFilters_isNoiseEcalBE << endl; return kTRUE;}
    if (NoiseFilters_isCSCTightHalo) {if(verbose) cout<< "NoiseFilters_isCSCTightHalo: " << NoiseFilters_isCSCTightHalo <<endl; return kTRUE;}
    if (NoiseFilters_isNoiseEEBadSc) {if(verbose) cout << "NoiseFilters_isNoiseEEBadSc: " << NoiseFilters_isNoiseEEBadSc <<endl; return kTRUE;}
    if (!NoiseFilters_isNoisetrkPOG1){if(verbose) cout << "!NoiseFilters_isNoisetrkPOG1: " << !NoiseFilters_isNoisetrkPOG1 << endl; return kTRUE;}
    if (!NoiseFilters_isNoisetrkPOG2){if(verbose) cout << "!NoiseFilters_isNoisetrkPOG2: " << !NoiseFilters_isNoisetrkPOG2 << endl; return kTRUE;}
    if (!NoiseFilters_isNoisetrkPOG3) {if(verbose) cout << "!NoiseFilters_isNoisetrkPOG3: " << !NoiseFilters_isNoisetrkPOG3 << endl;return kTRUE;}
    pass_filter++;
    
  }

  
  ////////////////////////////
  //Check the event vertices//  
  //////////////////////////// 
  if(verbose) cout << "Checking for good vertices..." << endl; 

  ntrks = -10;
  int ngVer = 0;
  vector<TVector3> goodVertices;
  vector<TCPrimaryVtx> pvVect;
  for (int i = 0; i < primaryVtx->GetSize(); ++i) {
    TCPrimaryVtx* pVtx = (TCPrimaryVtx*) primaryVtx->At(i);
    if (!pVtx->IsFake()
	&& pVtx->NDof() > 4.
	&& fabs(pVtx->z()) <= 24.
	&& fabs(pVtx->Perp()) <= 2.
	){
      		goodVertices.push_back(*pVtx);
		if( ngVer == 0 ) ntrks = pVtx->Ntracks();
		ngVer++;
	 }
  }
  if (goodVertices.size() < 1) return kTRUE; 
  
  pvPosition = new TVector3();
  *pvPosition = goodVertices[0];
 
  n_vertex_count++;
 
  pfMet = TVector2(corrMET->Px(), corrMET->Py());
  pfMET = &pfMet;

  if(verbose) cout << "Checking for good photons..." << endl; 
 
//Filling Photons and eFakes

  for (Int_t i = 0; i < recoPhotons->GetSize(); ++i) {
    TCPhoton* thisPhoton = (TCPhoton*) recoPhotons->At(i);
    if( PassPhotonID(thisPhoton, mediumPhID) && PassPhotonIso(thisPhoton,mediumPhIso) && thisPhoton->Pt() > 5.){
	Good_Photons.push_back(*thisPhoton);
    }
    if( PassPhotonID_eFake(thisPhoton, mediumPhID) && PassPhotonIso(thisPhoton,mediumPhIso) && thisPhoton->Pt() > 5.){
	Good_eFakes.push_back(*thisPhoton);
    }
  }

//Filling Electrons

  for (Int_t i = 0; i < recoElectrons->GetSize(); i++) {
     TCElectron* thisElec = (TCElectron*) recoElectrons->At(i);
     if( PassElectronID(thisElec, looseElID) && PassElectronIso(thisElec, looseElIso) && thisElec->Pt() > 5.){
	Good_Electrons.push_back(*thisElec);
     } 
  }

//Filling Z_ee
if(Good_Electrons.size() > 1){
//  cout << "NUMBER OF ELECTRONS > 1!!" << endl;
  for (Int_t i = 0; i < Good_Electrons.size(); i++){
     for (Int_t j = i+1; j < Good_Electrons.size(); j++){
	TLorentzVector ZTemp = Good_Electrons[i] + Good_Electrons[j];
	TLorentzVector* ztemp = new TLorentzVector(ZTemp.Px(), ZTemp.Py(), ZTemp.Pz(), ZTemp.E());
	Zee = ztemp;
	Zee_E1 = &Good_Electrons[i];
	Zee_E2 = &Good_Electrons[j];
	TZee->Fill();
     }
  }
}
//Filling Z_pp
if(Good_Photons.size() > 1){
  for (Int_t i = 0; i < Good_Photons.size(); i++){
     for (Int_t j = i+1; j < Good_Photons.size(); j++){
	TLorentzVector ZTemp = Good_Photons[i] + Good_Photons[j];
	TLorentzVector* ztemp = new TLorentzVector(ZTemp.Px(), ZTemp.Py(), ZTemp.Pz(), ZTemp.E());
	Zpp = ztemp;
	Zpp_P1 = &Good_Photons[i];
	Zpp_P2 = &Good_Photons[j];
	TZpp->Fill();
     }
  }
}
//Filling Z_pepe
if(Good_eFakes.size() > 1){
  for (Int_t i = 0; i < Good_eFakes.size(); i++){
     for (Int_t j = i+1; j < Good_eFakes.size(); j++){
	TLorentzVector ZTemp = Good_eFakes[i] + Good_eFakes[j];
	TLorentzVector* ztemp = new TLorentzVector(ZTemp.Px(), ZTemp.Py(), ZTemp.Pz(), ZTemp.E());
	Zkk = ztemp;
	Zkk_K1 = &Good_eFakes[i];
	Zkk_K2 = &Good_eFakes[j];
	TZkk->Fill();
     }
  }
}
//Filling Zp_e
if(Good_Photons.size() > 0 && Good_Electrons.size() > 0) {
  for (Int_t i = 0; i < Good_Photons.size(); i++){
     for (Int_t j = 0; j < Good_Electrons.size(); j++){
	TLorentzVector ZTemp = Good_Photons[i] + Good_Electrons[j];
	TLorentzVector* ztemp = new TLorentzVector(ZTemp.Px(), ZTemp.Py(), ZTemp.Pz(), ZTemp.E());
	Zep = ztemp;
	Zep_E = &Good_Electrons[j];
	Zep_P = &Good_Photons[i];
	TZep->Fill();
     }
  }
}
//Filling Zp_pe
if(Good_Photons.size() > 0 && Good_eFakes.size() > 0) {
  for (Int_t i = 0; i < Good_Photons.size(); i++){
     for (Int_t j = 0; j < Good_eFakes.size(); j++){
	TLorentzVector ZTemp = Good_Photons[i] + Good_eFakes[j];
	TLorentzVector* ztemp = new TLorentzVector(ZTemp.Px(), ZTemp.Py(), ZTemp.Pz(), ZTemp.E());
	Zpk = ztemp;
	Zpk_K = &Good_eFakes[j];
	Zpk_P = &Good_Photons[i];
	TZpk->Fill();
     }
  }
}
//Filling Zpe_e
if(Good_eFakes.size() > 0 && Good_Electrons.size() > 0) {
  for (Int_t i = 0; i < Good_eFakes.size(); i++){
     for (Int_t j = 0; j < Good_Electrons.size(); j++){
	TLorentzVector ZTemp = Good_eFakes[i] + Good_Electrons[j];
	TLorentzVector* ztemp = new TLorentzVector(ZTemp.Px(), ZTemp.Py(), ZTemp.Pz(), ZTemp.E());
	Zek = ztemp;
	Zek_K = &Good_eFakes[i];
	Zek_E = &Good_Electrons[j];
	TZek->Fill();
     }
  }
}
  return kTRUE;
}

void eFakeSel::SlaveTerminate()
{
}

void eFakeSel::Terminate()
{
  histoFile->cd();
  AnDir->cd();
  TZee->Write();
  TZep->Write();
  TZek->Write();
  TZpp->Write();
  TZpk->Write();
  TZkk->Write();
}


bool eFakeSel::PassQCDID(TCPhoton *ph){

  //Iso:
  float chEA,nhEA,phEA,chIsoCor,nhIsoCor,phIsoCor,tmpEta;
  
  if (fabs(tmpEta) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }

  chIsoCor = ph->IsoMap("chIso03")-rhoFactor*chEA;
  nhIsoCor = ph->IsoMap("nhIso03")-rhoFactor*nhEA;
  phIsoCor = ph->IsoMap("phIso03")-rhoFactor*phEA;

  //Id

  bool denoID=false;
              
  bool upperBound=false;
  bool lowerBound =false;

  tmpEta = ph->SCEta();

  double  maxPFCharged= TMath::Min(5.0*(2.6) , 0.20*ph->Pt());
  double  maxPFPhoton = TMath::Min(5.0*(1.3+0.005*ph->Pt()) , 0.20*ph->Pt());
  double  maxPFNeutral= TMath::Min(5.0*(3.5+0.04*ph->Pt()) , 0.20*ph->Pt());


  upperBound= (max((double)chIsoCor,0.) < maxPFCharged
	       && max((double)nhIsoCor,0.) < maxPFNeutral
	       && max((double)phIsoCor,0.) < maxPFPhoton
	       && ph->HadOverEm() < 0.05
	       && ph->SigmaIEtaIEta() < 0.015
	       && ph->ConversionVeto() == 1);
  lowerBound= (max((double)chIsoCor,0.) > 2.6
	       || max((double)nhIsoCor,0.) > 3.5+0.04 *ph->Pt()
	       || max((double)phIsoCor,0.) > 1.3+0.005*ph->Pt());
  
  if( (upperBound) && (lowerBound) )denoID = true;

  return denoID; 

}


bool eFakeSel::PassPhotonID(TCPhoton *ph, phIDCuts cutLevel){
  float tmpEta;
  bool phoPass = false;
  tmpEta = ph->SCEta();
  if (fabs(tmpEta) > 2.5) return phoPass;
  if (fabs(tmpEta) > 1.4442 && fabs(tmpEta) < 1.566) return phoPass;
  if(
     (fabs(tmpEta)  < 1.4442
      //&& ph->ConversionVeto()       == 0
      && ph->ConversionVeto()       == 1 //cutLevel.PassedEleSafeVeto[0]
      && ph->TrackVeto()	    == 0	
      && ph->HadOverEm()               < cutLevel.HadOverEm[0]
      && ph->SigmaIEtaIEta()           < cutLevel.sigmaIetaIeta[0]
      ) ||
     (fabs(tmpEta)  > 1.566
      //&& ph->ConversionVeto()       == 0
      && ph->ConversionVeto()       == 1 //cutLevel.PassedEleSafeVeto[1]
      && ph->TrackVeto()	    == 0
      && ph->HadOverEm()               < cutLevel.HadOverEm[1]
      && ph->SigmaIEtaIEta()           < cutLevel.sigmaIetaIeta[1]
      )
     ) phoPass = true;
  return phoPass;
}

bool eFakeSel::PassPhotonID_eFake(TCPhoton *ph, phIDCuts cutLevel){
  float tmpEta;
  bool phoPass = false;
  bool ConvNotPass = false;
  tmpEta = ph->SCEta();
  if (fabs(tmpEta) > 2.5) return phoPass;
  if (fabs(tmpEta) > 1.4442 && fabs(tmpEta) < 1.566) return phoPass;
  if( ph->ConversionVeto() == 0 || ph->TrackVeto() == 1 ) ConvNotPass = true;
  if(
     (fabs(tmpEta)  < 1.4442
      //&& ph->ConversionVeto()       == 0
      && ConvNotPass		       == true
      && ph->HadOverEm()               < cutLevel.HadOverEm[0]
      && ph->SigmaIEtaIEta()           < cutLevel.sigmaIetaIeta[0]
      ) ||
     (fabs(tmpEta)  > 1.566
      //&& ph->ConversionVeto()       == 0
      && ConvNotPass		       == true
      && ph->HadOverEm()               < cutLevel.HadOverEm[1]
      && ph->SigmaIEtaIEta()           < cutLevel.sigmaIetaIeta[1]
      )
     ) phoPass = true;
  return phoPass;
}

bool eFakeSel::PassPhotonIso(TCPhoton *ph, phIsoCuts cutLevel){
  float chEA,nhEA,phEA,chIsoCor,nhIsoCor,phIsoCor,tmpEta;
  bool isoPass = false;
  tmpEta = ph->SCEta();

  if(fabs(tmpEta) > 2.5) return isoPass;

  if (fabs(tmpEta) < 1.0){
    chEA = EAPho[0][0];
    nhEA = EAPho[0][1];
    phEA = EAPho[0][2];
  }else if (fabs(tmpEta) < 1.479){
    chEA = EAPho[1][0];
    nhEA = EAPho[1][1];
    phEA = EAPho[1][2];
  }else if (fabs(tmpEta) < 2.0){
    chEA = EAPho[2][0];
    nhEA = EAPho[2][1];
    phEA = EAPho[2][2];
  }else if (fabs(tmpEta) < 2.2){
    chEA = EAPho[3][0];
    nhEA = EAPho[3][1];
    phEA = EAPho[3][2];
  }else if (fabs(tmpEta) < 2.3){
    chEA = EAPho[4][0];
    nhEA = EAPho[4][1];
    phEA = EAPho[4][2];
  }else if (fabs(tmpEta) < 2.4){
    chEA = EAPho[5][0];
    nhEA = EAPho[5][1];
    phEA = EAPho[5][2];
  }else{
    chEA = EAPho[6][0];
    nhEA = EAPho[6][1];
    phEA = EAPho[6][2];
  }

  chIsoCor = ph->IsoMap("chIso03")-rhoFactor*chEA;
  nhIsoCor = ph->IsoMap("nhIso03")-rhoFactor*nhEA;
  phIsoCor = ph->IsoMap("phIso03")-rhoFactor*phEA;

  if (cutLevel.cutName == "loosePhIso"){
    if (
        (fabs(tmpEta) < 1.4442
	 //(fabs(ph->Eta())  < 1.566                                                                                                                                                
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[0]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[0] + 0.04*ph->Pt()
         && max((double)phIsoCor,0.)          < cutLevel.phIso03[0] + 0.005*ph->Pt()
	 ) ||
	(fabs(tmpEta) > 1.566
	 //(fabs(ph->Eta())  > 1.566                                                                                                                                                
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[1]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[1] + 0.04*ph->Pt()
         //&& phoCut["phIso03"]/ph->Pt() < nuthin                                                                                                                                  
	 )
	) isoPass = true;
  } else {
    if (
	//(fabs(ph->Eta())  < 1.566                                                                                                                                                
        (fabs(tmpEta) < 1.4442
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[0]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[0] + 0.04*ph->Pt()
         && max((double)phIsoCor,0.)          < cutLevel.phIso03[0] + 0.005*ph->Pt()
	 ) ||
        //(fabs(ph->Eta())  > 1.566                                                                                                                                                
        (fabs(tmpEta) > 1.566
         && max((double)chIsoCor,0.)          < cutLevel.chIso03[1]
         && max((double)nhIsoCor,0.)          < cutLevel.nhIso03[1] + 0.04*ph->Pt()
         && max((double)phIsoCor,0.)          < cutLevel.phIso03[1] + 0.005*ph->Pt()
	 )
	) isoPass = true;
  }
  return isoPass;
}

bool eFakeSel::PassElectronIso(TCElectron *el, elIsoCuts cutLevel){
  float thisEA = 0;
  if (fabs(el->Eta())     <  1.0) thisEA = EAEle[0];
  else if (fabs(el->Eta())     <  1.5) thisEA = EAEle[1];
  else if (fabs(el->Eta())     <  2.0) thisEA = EAEle[2];
  else if (fabs(el->Eta())     <  2.2) thisEA = EAEle[3];
  else if (fabs(el->Eta())     <  2.3) thisEA = EAEle[4];
  else if (fabs(el->Eta())     <  2.4) thisEA = EAEle[5];
  else if (fabs(el->Eta())     >  2.4) thisEA = EAEle[6];

  float combIso = (el->IsoMap("pfChIso_R04")
    + max(0.,(double)el->IsoMap("pfNeuIso_R04") + el->IsoMap("pfPhoIso_R04") - rhoFactor*thisEA));
  bool isoPass = false;
  if (combIso/el->Pt() < cutLevel.relCombIso04) isoPass = true;
  return isoPass;
}


bool eFakeSel::PassMuonID(TCMuon *mu, muIDCuts cutLevel){

  bool muPass = false;
  if (
    fabs(mu->Eta()) < 2.4
      && mu->IsPF()                          == cutLevel.IsPF
      && mu->IsGLB()                         == cutLevel.IsGLB
      && mu->NormalizedChi2()                < cutLevel.NormalizedChi2
      && mu->NumberOfValidMuonHits()         > cutLevel.NumberOfValidMuonHits
      && mu->NumberOfMatchedStations()       > cutLevel.NumberOfMatchedStations
      && mu->NumberOfValidPixelHits()        > cutLevel.NumberOfValidPixelHits
      && mu->TrackLayersWithMeasurement()    > cutLevel.TrackLayersWithMeasurement
      && fabs(mu->Dxy(pvPosition))           < cutLevel.dxy
      && fabs(mu->Dz(pvPosition))            < cutLevel.dz
    ) muPass = true;
  return muPass;

}

bool eFakeSel::PassMuonIso(TCMuon *mu, muIsoCuts cutLevel){

  float combIso;

  combIso = (mu->IsoMap("pfChargedHadronPt_R04")
    + max(0.,(double)mu->IsoMap("pfNeutralHadronEt_R04") + mu->IsoMap("pfPhotonEt_R04") - 0.5*mu->IsoMap("pfPUPt_R04")));

  bool isoPass = false;
  if (combIso/mu->Pt() < cutLevel.relCombIso04) isoPass = true;
  return isoPass;
}



bool eFakeSel::PassElectronID(TCElectron *el, elIDCuts cutLevel)
{
  bool elPass = false;
  if (fabs(el->SCEta()) > 2.5) return elPass;
  if (fabs(el->SCEta()) > 1.4442 && fabs(el->SCEta()) < 1.566) return elPass;
  if (
      (fabs(el->Eta()) < 1.566
       && fabs(el->SCDeltaEta())    < cutLevel.dEtaIn[0]
       && fabs(el->SCDeltaPhi())    < cutLevel.dPhiIn[0]
       && el->SigmaIEtaIEta()             < cutLevel.sigmaIetaIeta[0]
       && el->HadOverEm()                 < cutLevel.HadOverEm[0]
       && fabs(el->Dxy(pvPosition))       < cutLevel.dxy[0]
       && fabs(el->Dz(pvPosition))        < cutLevel.dz[0]
       && el->IdMap("fabsEPDiff")         < cutLevel.fabsEPDiff[0]
       && el->ConversionMissHits()        <= cutLevel.ConversionMissHits[0]
       && el->PassConversionVeto()        == cutLevel.PassedConversionProb[0]
       ) ||
    (fabs(el->Eta()) > 1.566  
      && fabs(el->SCDeltaEta())    < cutLevel.dEtaIn[1]
      && fabs(el->SCDeltaPhi())    < cutLevel.dPhiIn[1]
      && el->SigmaIEtaIEta()             < cutLevel.sigmaIetaIeta[1]
      && el->HadOverEm()                 < cutLevel.HadOverEm[1]
      && fabs(el->Dxy(pvPosition))       < cutLevel.dxy[1]
      && fabs(el->Dz(pvPosition))        < cutLevel.dz[1]
      && el->IdMap("fabsEPDiff")         < cutLevel.fabsEPDiff[1]
      && el->ConversionMissHits()        <= cutLevel.ConversionMissHits[1]
     && el->PassConversionVeto()         == cutLevel.PassedConversionProb[1]
      )
    ) elPass = true; 
  for (int j = 0; j < recoMuons->GetSize(); ++ j)
       {
         TCMuon* thisMuon = (TCMuon*) recoMuons->At(j);    
         if (thisMuon->DeltaR(*el) < 0.05){
           elPass = false;
           break;
         }
       }
       return elPass;
}


double  eFakeSel::func(const double* par){
/*
  Ndim = reco_pt.size();
  double px = 0, py =0, arg = 0;  
  for(int i=0; i<Ndim; i++){
    px += par[i]*cos(reco_phi[i]);
    py += par[i]*sin(reco_phi[i]);
    arg += pow((reco_pt[i]-par[i])/(sigma_ms[i]),2);
  }
  return arg + ((px*px)/(sigma_mex*sigma_mex) + (py*py)/(sigma_mey*sigma_mey));
*/
return 0;
}

double  eFakeSel::jetreso(TF1 *fjetreso_T, double pt, double eta_or){
  double reso = 0;
  //TF1 *fjetreso = new TF1 ("jetReso","TMath::Sqrt(((sign([0])*(([0]/x)^2))+(([1]^2)*(x^([2]-1.))))+([3]^2))",0,15000);
  double scale = 0;

  double eta = TMath::Abs(eta_or);

  if (eta <= 0.3)              {fjetreso_T->SetParameters(2.866,0.3118,0.4075,0.01823); scale = 1.052;}
  if (eta > 0.3 && eta <= 0.5) {fjetreso_T->SetParameters(2.91,0.2793,0.4629,0.001049); scale = 1.052;}
  if (eta > 0.5 && eta <= 0.8) {fjetreso_T->SetParameters(2.768,0.3797,0.3144,0.02803); scale = 1.057;}
  if (eta > 0.8 && eta <= 1.1) {fjetreso_T->SetParameters(2.934,0.3251,0.4401,0.0079); scale = 1.057;}
  if (eta > 1.1 && eta <= 1.4) {fjetreso_T->SetParameters(2.617,0.736,0.0899,-0.04179); scale = 1.096;}
  if (eta > 1.4 && eta <= 1.7) {fjetreso_T->SetParameters(0.1406,1.477,-0.2062,-0.03656); scale = 1.096;}
  if (eta > 1.7 && eta <= 2.0) {fjetreso_T->SetParameters(1.959,1.099,-0.1357,-0.02382); scale = 1.134;}
  if (eta > 2.0 && eta <= 2.3) {fjetreso_T->SetParameters(4.113,0.4146,0.1918,0.02413); scale = 1.134;}
  if (eta > 2.3 && eta <= 2.8) {fjetreso_T->SetParameters(5.817,0.1547,0.5529,0.001136); scale = 1.288;}
  if (eta > 2.8 && eta <= 3.2) {fjetreso_T->SetParameters(4.894,0.3666,0.4251,-0.00215); scale = 1.288;}
  if (eta > 3.2 && eta <= 4.1) {fjetreso_T->SetParameters(3.624,0.2542,0.6046,0.02232); scale = 1.288;}
  if (eta > 4.1)               {fjetreso_T->SetParameters(2.727,1.035,-0.1662,0); scale = 1.288;}


  if(mc)   reso = pt * fjetreso_T->Eval(pt);
  else reso = pt * scale * fjetreso_T->Eval(pt);
  return reso;

}


float eFakeSel::computeLICDT(TCPhoton *ph){
  
  vector<TCPhoton::CrystalInfo> savedCrystals = ph->GetCrystalVect();

  int seed = -1;
  float seedE = -999999, LICTD = 0;
  
  for(int k=0;k<ph->GetNCrystals() && k < 100; ++k){
    if(savedCrystals[k].energy > seedE){seed = k; seedE = savedCrystals[seed].energy;}
    if(seed<0) LICTD = -99;
    if(seed==k) continue;
    if (savedCrystals[k].energy > 1.) {
      if(TMath::Abs(savedCrystals[seed].time - savedCrystals[k].time) > TMath::Abs(LICTD)) LICTD = savedCrystals[seed].time - savedCrystals[k].time;
    }
  }
  return LICTD;
  savedCrystals.clear();
}
