#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TCut.h"
#include <iostream>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
using namespace std;

void histogram_maker_plotter_categorization()
{
  bool is_gpart = true;
  //GParT
  double L_WP[] = {0.95,0.975,0.99, 0.95};
  double H_WP[] = {0.975,0.99,1.0, 1.0};
  //PNet
  //double L_WP[] = {0.955, 0.98, 0.992};
  //double H_WP[] = {0.98, 0.992, 1.0};
  TString CAT[] = {"CAT3", "CAT2","CAT1","CAT0"};

for (int jj = 0; jj < (int)(sizeof(CAT) / sizeof(CAT[0])); jj++)
  {
        std::cout << "Running on categories: " << CAT[jj] << std::endl;
        // TString input_dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zbbg/";
	TString input_dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zbbg";
	
        TString subdir;
	if(is_gpart) subdir = "GParT/";
	else subdir = "PNet/";

        // TString output_dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zbbg/Histogram/"+subdir+CAT[jj]+"/";
        TString output_dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zbbg/Histogram/"+subdir+CAT[jj]+"/";
	
        // 2022 / 2022EE
	TString era = "2022";

	double LUMI = 35.00;
	Double_t tot_weight[] = {
		(1.18502e+10/5.0),1.98309e+07,1.43534e+06,
		1.39551e+06,2.85479e+06,
		3.97481e+07, 4.85838e+06, 1.19719e+06, 667596,	
		1.0
	};

	TString proc_name[] = {
	    "GJ_PTG-200to400_TuneCP5_13p6TeV_amcatnlo-pythia8",
	    "GJ_PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8",
	    "GJ_PTG-600_TuneCP5_13p6TeV_amcatnlo-pythia8",
	    "ZGto2QG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
	    "WGto2QG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
	    "Zto2Q-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8", 
	    "Zto2Q-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8",    
	    "Zto2Q-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8",
	    "Zto2Q-4Jets_HT-800_TuneCP5_13p6TeV_madgraphMLM-pythia8",  
	    "EGamma" // "EGamma_ZbbG"
	};



const int nproc = sizeof(proc_name)/sizeof(proc_name[0])-1;
const int nwgt = sizeof(tot_weight)/sizeof(tot_weight[0])-1;
if(nproc != nwgt)  std::cout << "Mismatch" << std::endl;

Bool_t debug = false;	

for(int kk = 0; kk < nproc + 1 ; kk++)
{
        std::cout << proc_name[kk] << "  " ;
        TFile* final_file = TFile::Open(output_dir +  proc_name[kk] +"_"+era+"_v2.root","RECREATE");

        //photon variables
        TH1F* hist_photon_eta     = new TH1F("hist_photon_eta","",50, -2.5, 2.5);           hist_photon_eta->Sumw2(); 
        TH1F* hist_photon_pt      = new TH1F("hist_photon_pt","",36, 200.0, 2000.0);        hist_photon_pt->Sumw2();
         
        //jet variables 
        TH1F* hist_ak8_mass       = new TH1F("hist_ak8_mass","",180, 20.0, 200.0);          hist_ak8_mass->Sumw2();
        TH1F* hist_ak8_mSD        = new TH1F("hist_ak8_mSD","",180, 20.0, 200.0);           hist_ak8_mSD->Sumw2();
        TH1F* hist_ak8_gpart_mass = new TH1F("hist_ak8_gpart_mass","",180, 20.0, 200.0);    hist_ak8_gpart_mass->Sumw2();
        TH1F* hist_ak8_pnet       = new TH1F("hist_ak8_pnet","",100, 0.0, 1.0);             hist_ak8_pnet->Sumw2();
        TH1F* hist_ak8_gpart      = new TH1F("hist_ak8_gpart","",100, 0.0, 1.0);            hist_ak8_gpart->Sumw2();
        TH1F* hist_ak8_pT         = new TH1F("hist_ak8_pT","",35, 100.0, 2000.0);           hist_ak8_pT->Sumw2();
        TH1F* hist_ak8_eta        = new TH1F("hist_ak8_eta","",50, -2.5, 2.5);              hist_ak8_eta->Sumw2();
        TH1F* hist_ak8_phi        = new TH1F("hist_ak8_phi","",50, -3.141, 3.141);          hist_ak8_phi->Sumw2();
 
        //jet-photon dr 
        TH1F* hist_photon_jet_dr  = new TH1F("hist_photon_jet_dr", "",300,0.0,3.0 );        hist_photon_jet_dr->Sumw2();
        //HT 
        TH1F* hist_HT             = new TH1F("hist_HT", "", 36, 200.0, 2000.0);             hist_HT->Sumw2();

        
        // build robust file path (ensure trailing slash)
        TString file_path = input_dir;
        if(!file_path.EndsWith("/")) file_path += "/";
        file_path += proc_name[kk] + TString("_") + era + TString(".root");

        TFile *file = TFile::Open(file_path);
        if(!file || file->IsZombie()) {
                std::cerr << "Error: could not open file '" << file_path << "'\n";
                if(file) { file->Close(); }
                continue;
        }

        TTree *mtree = dynamic_cast<TTree*>(file->Get("tree"));
        if(!mtree) {
                // try to find any TTree in the file as a fallback
                TIter next(file->GetListOfKeys());
                TKey *key;
                while((key = (TKey*)next())) {
                        TObject *obj = key->ReadObj();
                        if(obj && obj->InheritsFrom("TTree")) {
                                mtree = dynamic_cast<TTree*>(obj);
                                std::cerr << "Warning: TTree 'tree' not found. Using '" << mtree->GetName() << "' from file '" << file->GetName() << "'\n";
                                break;
                        }
                }
                if(!mtree) {
                        std::cerr << "Error: no TTree found in file '" << file->GetName() << "'\n";
                        file->Close();
                        continue;
                }
        }

        Float_t         wgt;
        Float_t         xsec;
        Float_t         weight;
        Float_t         HT;
        Float_t         AK8_ptJet0;
        Float_t         AK8_etaJet0;
        Float_t         AK8_phiJet0;
        Float_t         AK8_sdmJet0;
        Float_t         AK8_legacy_mass0;
        Float_t         AK8_pnet_mass0;
        Float_t         AK8_gpart_mass0;
        Float_t         AK8_glopart_probHbbVsQCDTT;
        Float_t         AK8_glopart_probHbbVsQCD;
        Float_t         AK8_pnet_probHbbVsQCD;
        Float_t         photon_pt;
        Float_t         photon_eta;
        Float_t         photon_phi;
        Bool_t          is_back_to_back;
        Int_t           Zbb;
        Int_t           Zcc;
        Int_t           Zqq;
        Float_t         VpTweight_QCD;
        Float_t         VpTweight_EWK;
        Float_t         PhoIDISO_SF;
        Float_t         puweight;
        Float_t         puweight_up;
        Float_t         puweight_dn;
        Bool_t          true_z_match;
        Float_t         PS_ISR_up;
        Float_t         PS_ISR_dn;
        Float_t         PS_FSR_up;
        Float_t         PS_FSR_dn;
        Float_t         weight_lead_id;
        Float_t         weight_subl_id;
        Float_t         weight_lead_reco;

        TBranch        *b_wgt;   //!
        TBranch        *b_xsec;   //!
        TBranch        *b_weight;   //!
        TBranch        *b_HT;   //!
        TBranch        *b_AK8_ptJet0;   //!
        TBranch        *b_AK8_etaJet0;   //!
        TBranch        *b_AK8_phiJet0;   //!
        TBranch        *b_AK8_sdmJet0;   //!
        TBranch        *b_AK8_legacy_mass0;   //!
        TBranch        *b_AK8_pnet_mass0;   //!
        TBranch        *b_AK8_gpart_mass0;   //!
        TBranch        *b_AK8_glopart_probHbbVsQCDTT;   //!
        TBranch        *b_AK8_glopart_probHbbVsQCD;   //!
        TBranch        *b_AK8_pnet_probHbbVsQCD;   //!
        TBranch        *b_photon_pt;   //!
        TBranch        *b_photon_eta;   //!
        TBranch        *b_photon_phi;   //!
        TBranch        *b_is_back_to_back;   //!
        TBranch        *b_Zbb;   //!
        TBranch        *b_Zcc;   //!
        TBranch        *b_Zqq;   //!
        TBranch        *b_VpTweight_QCD;   //!
        TBranch        *b_VpTweight_EWK;   //!
        TBranch        *b_PhoIDISO_SF;   //!
        TBranch        *b_puweight;   //!
        TBranch        *b_puweight_up;   //!
        TBranch        *b_puweight_dn;   //!
        TBranch        *b_true_z_match;   //!
        TBranch        *b_PS_ISR_up;   //!
        TBranch        *b_PS_ISR_dn;   //!
        TBranch        *b_PS_FSR_up;   //!
        TBranch        *b_PS_FSR_dn;   //!
        TBranch        *b_weight_lead_id;   //!
        TBranch        *b_weight_subl_id;   //!
        TBranch        *b_weight_lead_reco;   //

        mtree->SetBranchAddress("weight", &weight, &b_weight);
        mtree->SetBranchAddress("HT", &HT, &b_HT);
        mtree->SetBranchAddress("AK8_ptJet0", &AK8_ptJet0, &b_AK8_ptJet0);
        mtree->SetBranchAddress("AK8_etaJet0", &AK8_etaJet0, &b_AK8_etaJet0);
        mtree->SetBranchAddress("AK8_phiJet0", &AK8_phiJet0, &b_AK8_phiJet0);
        mtree->SetBranchAddress("AK8_sdmJet0", &AK8_sdmJet0, &b_AK8_sdmJet0);
        mtree->SetBranchAddress("AK8_legacy_mass0", &AK8_legacy_mass0, &b_AK8_legacy_mass0);
        mtree->SetBranchAddress("AK8_pnet_mass0", &AK8_pnet_mass0, &b_AK8_pnet_mass0);
        mtree->SetBranchAddress("AK8_gpart_mass0", &AK8_gpart_mass0, &b_AK8_gpart_mass0);
        mtree->SetBranchAddress("AK8_glopart_probHbbVsQCDTT", &AK8_glopart_probHbbVsQCDTT, &b_AK8_glopart_probHbbVsQCDTT);
        mtree->SetBranchAddress("AK8_glopart_probHbbVsQCD", &AK8_glopart_probHbbVsQCD, &b_AK8_glopart_probHbbVsQCD);
        mtree->SetBranchAddress("AK8_pnet_probHbbVsQCD", &AK8_pnet_probHbbVsQCD, &b_AK8_pnet_probHbbVsQCD);
        mtree->SetBranchAddress("photon_pt", &photon_pt, &b_photon_pt);
        mtree->SetBranchAddress("photon_eta", &photon_eta, &b_photon_eta);
        mtree->SetBranchAddress("photon_phi", &photon_phi, &b_photon_phi);
        mtree->SetBranchAddress("is_back_to_back", &is_back_to_back, &b_is_back_to_back);

        if(proc_name[kk]!="EGamma_ZbbG")
                {
                        mtree->SetBranchAddress("wgt", &wgt, &b_wgt);
                        mtree->SetBranchAddress("xsec", &xsec, &b_xsec);

                        mtree->SetBranchAddress("Zbb", &Zbb, &b_Zbb);
                        mtree->SetBranchAddress("Zcc", &Zcc, &b_Zcc);
                        mtree->SetBranchAddress("Zqq", &Zqq, &b_Zqq);
                        mtree->SetBranchAddress("VpTweight_QCD", &VpTweight_QCD, &b_VpTweight_QCD);
                        mtree->SetBranchAddress("VpTweight_EWK", &VpTweight_EWK, &b_VpTweight_EWK);
                        mtree->SetBranchAddress("PhoIDISO_SF", &PhoIDISO_SF, &b_PhoIDISO_SF);
                        mtree->SetBranchAddress("puweight", &puweight, &b_puweight);
                        mtree->SetBranchAddress("puweight_up", &puweight_up, &b_puweight_up);
                        mtree->SetBranchAddress("puweight_dn", &puweight_dn, &b_puweight_dn);
                        mtree->SetBranchAddress("true_z_match", &true_z_match, &b_true_z_match);
                        mtree->SetBranchAddress("weight_lead_id", &weight_lead_id, &b_weight_lead_id);
                        mtree->SetBranchAddress("weight_subl_id", &weight_subl_id, &b_weight_subl_id);
                        mtree->SetBranchAddress("weight_lead_reco", &weight_lead_reco, &b_weight_lead_reco);
                }

                Long64_t nn; 
                Float_t nnlo_kfactor = 0.94;
                nn = mtree->GetEntries();

                for(Long64_t j =0; j < nn ; j++)
                {
                        mtree->GetEntry(j);
                        if(HT < 500) continue;
                        // if(AK8_etaJet0 > 450) continue;
                        if(AK8_ptJet0 < 300) continue;
                        if (proc_name[kk] == "EGamma_ZbbG") weight = 1.0;
                        else if (proc_name[kk].Index("Zto2Q") == 0) weight = (weight * PhoIDISO_SF * nnlo_kfactor * puweight * VpTweight_QCD * VpTweight_EWK * LUMI ) / tot_weight[kk];
                        else if (proc_name[kk].Index("ZGto2QG") == 0 || proc_name[kk].Index("WGto2QG") == 0 )weight = (weight * PhoIDISO_SF * nnlo_kfactor * puweight * VpTweight_EWK * LUMI ) / tot_weight[kk];
                        else weight = (weight * PhoIDISO_SF * puweight * LUMI ) / tot_weight[kk];

                        Float_t b_tag_disc;
                        if(!is_gpart) b_tag_disc = AK8_pnet_probHbbVsQCD;
                        else b_tag_disc = AK8_glopart_probHbbVsQCD;
                        
                        if(b_tag_disc  > L_WP[jj] && b_tag_disc < H_WP[jj] )
                        {
                                hist_photon_eta->Fill(photon_eta, weight);
                                hist_photon_pt->Fill(photon_pt, weight);
                                hist_HT->Fill(HT, weight);
                                hist_ak8_mSD->Fill(AK8_sdmJet0, weight);
                                hist_ak8_mass->Fill(AK8_legacy_mass0, weight);
                                hist_ak8_gpart_mass->Fill(AK8_gpart_mass0, weight);
                                hist_ak8_pnet->Fill(AK8_pnet_probHbbVsQCD, weight);
                                hist_ak8_gpart->Fill(AK8_glopart_probHbbVsQCD, weight);
                                hist_ak8_pT->Fill(AK8_ptJet0, weight);
                                hist_ak8_eta->Fill(AK8_etaJet0, weight);
                                hist_ak8_phi->Fill(AK8_phiJet0, weight);
                                TLorentzVector jet; jet.SetPtEtaPhiM(AK8_ptJet0,AK8_etaJet0,AK8_phiJet0,AK8_gpart_mass0);
                                TLorentzVector pho; pho.SetPtEtaPhiM(photon_pt,photon_eta,photon_phi,0.0);
                                hist_photon_jet_dr->Fill(jet.DeltaR(pho),weight);
                        }
                }
                        //event loopDYto2L-4Jets
                        final_file->Write();
                        final_file->cd();
                        std::cout << hist_ak8_gpart_mass->Integral() << std::endl;

                        //if(proc_name[kk]=="Muon") std::cout << hist_dimu_mass->Integral() << std::endl;
                        hist_photon_eta ->Write();
                        hist_photon_pt  ->Write();
                        hist_ak8_mSD    ->Write(); 
                        hist_ak8_mass   ->Write();
                        hist_ak8_gpart_mass->Write();
                        hist_ak8_pnet   ->Write();
                        hist_ak8_gpart  ->Write();
                        hist_ak8_pT     ->Write();
                        hist_ak8_eta    ->Write();
                        hist_ak8_phi    ->Write();
                        hist_photon_jet_dr->Write();
                        hist_HT         ->Write();
                final_file->Close();
        }//kk
  }
}

         
