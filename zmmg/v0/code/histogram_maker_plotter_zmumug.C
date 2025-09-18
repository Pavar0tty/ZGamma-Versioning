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
#include "TCut.h"
#include <iostream>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
using namespace std;

void histogram_maker_plotter_zmumug()
{
// TString input_dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zmmg/";
// TString output_dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zmmg/Histogram/";

TString input_dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zmmg/";
TString output_dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zmmg/Histogram/";

//2022EE
TString era = "2022EE";
double LUMI = 27.00;
Double_t tot_weight[] = {
        235227, 14627, 2667.39, 
	2.00495e+11, 8.20022e+10, 2.12091e+11, 2.82106e+08, 4.93834e+07, 3.65275e+06, 250513,
	1.0
};

TString proc_name[] = {
    "DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
    "DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
    "DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",    
    "DYto2L-4Jets_MLL-50to120_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50to120_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50to120_HT-100to400_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50to120_HT-800to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50to120_HT-1500to2500_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50to120_HT-2500_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "EGamma"
};


const int nproc = sizeof(proc_name)/sizeof(proc_name[0])-1;
const int nwgt = sizeof(tot_weight)/sizeof(tot_weight[0])-1;
if(nproc != nwgt)  std::cout << "Mismatch" << std::endl;

for(int kk = 0; kk < nproc + 1 ; kk++)
{
Bool_t debug = false;	
          std::cout << proc_name[kk] << "  " ;
	  TFile* final_file = TFile::Open(output_dir +  proc_name[kk] +"_"+era+"_v2.root","RECREATE");
	  TH1F* hist_photon_eta       = new TH1F("hist_photon_eta","",50, -2.5, 2.5);           hist_photon_eta->Sumw2(); 
	  TH1F* hist_photon_pt        = new TH1F("hist_photon_pt","",36, 200.0, 2000.0);        hist_photon_pt->Sumw2();
          
	  TH1F* hist_dimu_pt          = new TH1F("hist_dimu_pt","",36, 200.0, 2000.0);          hist_dimu_pt->Sumw2();
	  TH1F* hist_dimu_eta         = new TH1F("hist_dimu_eta","",50, -2.5, 2.5);             hist_dimu_eta->Sumw2();
	  TH1F* hist_dimu_mass        = new TH1F("hist_dimu_mass","", 100, 70.0, 110.0);        hist_dimu_mass->Sumw2();
          TH1F* hist_dimu_dr          = new TH1F("hist_dimu_dr","",200,0.0,2.0);                hist_dimu_dr->Sumw2();
          TH1F* hist_leadmu_photon_dr = new TH1F("hist_leadmu_photon_dr", "",300,0.0,3.0 );     hist_leadmu_photon_dr->Sumw2();
          TH1F* hist_sublmu_photon_dr = new TH1F("hist_sublmu_photon_dr", "",300,0.0,3.0 );     hist_sublmu_photon_dr->Sumw2();

	  TH1F* hist_lead_pt          = new TH1F("hist_lead_pt", "", 20, 0.0, 500.0);           hist_lead_pt->Sumw2();
	  TH1F* hist_subl_pt          = new TH1F("hist_subl_pt", "", 20, 0.0, 500.0);           hist_subl_pt->Sumw2();
    
	  TH1F* hist_lead_eta         = new TH1F("hist_lead_eta", "", 50, -2.5, 2.5);           hist_lead_eta->Sumw2();
	  TH1F* hist_subl_eta         = new TH1F("hist_subl_eta", "", 50, -2.5, 2.5);           hist_subl_eta->Sumw2();
          TH1F* hist_HT               = new TH1F("hist_HT", "", 36, 200.0, 2000.0);             hist_HT->Sumw2();


              TFile *file = TFile::Open(input_dir+ proc_name[kk] + "_"+era+".root");
              if(!file || file->IsZombie()) {
                        std::cerr << "❓Error: could not open file '" << (input_dir + proc_name[kk] + "_" + era + ".root") << "❓'\n";
                        // ensure file pointer is cleaned up if partially opened
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
                                        std::cerr << "⚠Warning: TTree 'tree' not found. Using '" << mtree->GetName() << "' from file '" << file->GetName() << "⚠'\n";
                                        break;
                                }
                                // if not a TTree, continue
                        }
                        if(!mtree) {
                                std::cerr << "❌Error: no TTree found in file '" << file->GetName() << "'❌\n";
                                file->Close();
                                continue;
                        }
              }


   Float_t         wgt;
   Float_t         xsec;
   Float_t         weight;
   Float_t         HT;
   Float_t         lead_muon_pt;
   Float_t         lead_muon_eta;
   Float_t         lead_muon_phi;
   Float_t         subl_muon_pt;
   Float_t         subl_muon_eta;
   Float_t         subl_muon_phi;
   Float_t         dimuon_cand_pt;
   Float_t         dimuon_cand_eta;
   Float_t         dimuon_cand_mass;
   Float_t         dimuon_cand_dr;
   Float_t         leadmu_photon_dr;
   Float_t         sublmu_photon_dr;
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
   Float_t         weight_lead_id;
   Float_t         weight_subl_id;
   Float_t         weight_lead_reco;

   TBranch        *b_wgt;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_lead_muon_pt;   //!
   TBranch        *b_lead_muon_eta;   //!
   TBranch        *b_lead_muon_phi;   //!
   TBranch        *b_subl_muon_pt;   //!
   TBranch        *b_subl_muon_eta;   //!
   TBranch        *b_subl_muon_phi;   //!
   TBranch        *b_dimuon_cand_pt;   //!
   TBranch        *b_dimuon_cand_eta;   //!
   TBranch        *b_dimuon_cand_mass;   //!
   TBranch        *b_dimuon_cand_dr;   //!
   TBranch        *b_leadmu_photon_dr;   //!
   TBranch        *b_sublmu_photon_dr;   //!
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
   TBranch        *b_weight_lead_id;   //!
   TBranch        *b_weight_subl_id;   //!
   TBranch        *b_weight_lead_reco;   //!
 
   mtree->SetBranchAddress("weight", &weight, &b_weight);
   mtree->SetBranchAddress("HT", &HT, &b_HT);
   mtree->SetBranchAddress("lead_muon_pt", &lead_muon_pt, &b_lead_muon_pt);
   mtree->SetBranchAddress("lead_muon_eta", &lead_muon_eta, &b_lead_muon_eta);
   mtree->SetBranchAddress("lead_muon_phi", &lead_muon_phi, &b_lead_muon_phi);
   mtree->SetBranchAddress("subl_muon_pt", &subl_muon_pt, &b_subl_muon_pt);
   mtree->SetBranchAddress("subl_muon_eta", &subl_muon_eta, &b_subl_muon_eta);
   mtree->SetBranchAddress("subl_muon_phi", &subl_muon_phi, &b_subl_muon_phi);
   mtree->SetBranchAddress("dimuon_cand_pt", &dimuon_cand_pt, &b_dimuon_cand_pt);
   mtree->SetBranchAddress("dimuon_cand_eta", &dimuon_cand_eta, &b_dimuon_cand_eta);
   mtree->SetBranchAddress("dimuon_cand_mass", &dimuon_cand_mass, &b_dimuon_cand_mass);
   mtree->SetBranchAddress("dimuon_cand_dr", &dimuon_cand_dr, &b_dimuon_cand_dr);
   mtree->SetBranchAddress("leadmu_photon_dr", &leadmu_photon_dr, &b_leadmu_photon_dr);
   mtree->SetBranchAddress("sublmu_photon_dr", &sublmu_photon_dr, &b_sublmu_photon_dr);
   mtree->SetBranchAddress("photon_pt", &photon_pt, &b_photon_pt);
   mtree->SetBranchAddress("photon_eta", &photon_eta, &b_photon_eta);
   mtree->SetBranchAddress("photon_phi", &photon_phi, &b_photon_phi);
   mtree->SetBranchAddress("is_back_to_back", &is_back_to_back, &b_is_back_to_back);
   if(proc_name[kk]!="EGamma")
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
	nn = mtree->GetEntries();
        for(Long64_t j =0; j < nn ; j++)
          {
                mtree->GetEntry(j);
		Float_t nnlo_kfactor = 0.94;
		//if(dimuon_cand_pt < 250 || dimuon_cand_pt > 450 ) continue;
		//if(dimuon_cand_pt < 450 || dimuon_cand_pt > 550 ) continue;
		//std::cout << weight << "  " << nnlo_kfactor << "  " << weight_subl_id << "  " << weight_lead_id << "  " << weight_lead_reco << "  " << VpTweight_QCD << "  " << VpTweight_EWK << std::endl;
		if (proc_name[kk] == "EGamma") weight = 1.0;
	        else if (proc_name[kk].Index("DYto2L") == 0) weight = (weight * PhoIDISO_SF * nnlo_kfactor * weight_subl_id * weight_lead_id * weight_lead_reco * LUMI * puweight * VpTweight_QCD * VpTweight_EWK ) / tot_weight[kk];
		else if (proc_name[kk].Index("DYGto2LG") == 0) weight = (weight * PhoIDISO_SF * nnlo_kfactor * weight_subl_id * weight_lead_id * weight_lead_reco * LUMI ) / tot_weight[kk];
		else weight = (weight * weight_subl_id * weight_lead_id * weight_lead_reco * LUMI * puweight) / tot_weight[kk];	
                
                hist_photon_eta->Fill(photon_eta, weight);
		hist_photon_pt->Fill(photon_pt, weight);
		hist_dimu_pt->Fill(dimuon_cand_pt, weight);
		hist_dimu_eta->Fill(dimuon_cand_eta, weight);
                hist_dimu_mass->Fill(dimuon_cand_mass, weight);
		hist_dimu_dr->Fill(dimuon_cand_dr, weight);
		hist_leadmu_photon_dr->Fill(leadmu_photon_dr, weight);
                hist_sublmu_photon_dr->Fill(sublmu_photon_dr, weight);
		hist_lead_pt->Fill(lead_muon_pt, weight);
                hist_subl_pt->Fill(subl_muon_pt, weight);
		hist_lead_eta->Fill(lead_muon_eta, weight);
                hist_subl_eta->Fill(subl_muon_eta, weight);
		hist_HT->Fill(HT, weight);
         }//event loopDYto2L-4Jets
          final_file->Write();
          final_file->cd();
	        std::cout << hist_dimu_mass->Integral() << std::endl;
	        //if(proc_name[kk]=="Muon") std::cout << hist_dimu_mass->Integral() << std::endl;
	        hist_photon_eta->Write();
	        hist_photon_pt->Write();
		hist_dimu_pt->Write();
		hist_dimu_eta->Write();
		hist_dimu_mass->Write();
		hist_lead_pt->Write();
		hist_subl_pt->Write();
		hist_lead_eta->Write();
		hist_subl_eta->Write();
		hist_HT->Write();
          final_file->Close();
    }//kk

}

         
