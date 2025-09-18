#include "TTree.h"
#include "TBrowser.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include <TStyle.h>
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include <vector>
#include "TLegend.h"
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include<iostream>
#include<cmath>
#include <fstream>
#include<TGraphAsymmErrors.h>
#include "CustomColors.h"
#include "CMS_lumi.C"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"

static TGraphAsymmErrors* staterror(TH1* h, TH1* /*tmpl*/){
    if(!h) return nullptr;
    int n = h->GetNbinsX();
    TGraphAsymmErrors* gr = new TGraphAsymmErrors(n);
    for(int i=1; i<=n; ++i){
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);
        double ex = h->GetBinWidth(i)/2.0;
        double ey = h->GetBinError(i);
        gr->SetPoint(i-1, x, y);
        gr->SetPointError(i-1, ex, ex, ey, ey);
    }
    return gr;
}

// Imposta a zero eventuali bin negativi (può capitare con campioni NLO con pesi negativi)
static void SanitizeForStack(TH1* h){
  if(!h) return;
  for(int i=1; i<=h->GetNbinsX(); ++i){
    if(h->GetBinContent(i) < 0){
      h->SetBinContent(i, 0.0);
    }
  }
}


void plotter_ZmmG()
{
   InitializeColors();
  TString era = "2022EE";
  Float_t lumi = 27.0; // pb^-1 scale consistent with histogram_maker

   // plotting dir
   TString plot_dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Plots/zmmg/";

   TString cat = "";
   //TString dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zbbg/Histogram/GParT/"+cat+"/";
  //TString dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zbbg/Histogram/";
  TString dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zmmg/Histogram/";

   //TString Var_Name[] = {"dimu_dr","photon_eta","photon_pt","dimu_mass","dimu_pt","dimu_eta","lead_pt","lead_eta","subl_pt","subl_eta","lead_eta","subl_eta","HT"};
   // Variabili disponibili dagli istogrammi prodotti da histogram_maker_plotter_zmumug.C
   TString Var_Name[] = {"photon_pt","photon_eta","dimu_mass","dimu_pt","dimu_eta","dimu_dr","leadmu_photon_dr","sublmu_photon_dr","lead_pt","subl_pt","lead_eta","subl_eta","HT"};
  bool normToData = false; // metti a true se vuoi scalare MC all'integrale dei dati (confronto di forma)
   
   const int nVars = sizeof(Var_Name)/sizeof(Var_Name[0]);
   for(int variable = 0; variable < nVars; variable++)
   {
    std::cout << Var_Name[variable] << std::endl;
    TString var = "hist_" + Var_Name[variable];
    // Helper per aprire e sommare istogrammi da una lista di file
    auto getSummedHisto = [&](const std::vector<TString>& files, const TString& hname, const char* outName){
        TH1F* sum = nullptr;
        for(const auto& f : files){
            TFile* tf = TFile::Open(dir + f);
            if(!tf || tf->IsZombie()){ std::cerr << "[warn] cannot open " << (dir+f) << std::endl; if(tf) tf->Close(); continue; }
            TH1F* h = (TH1F*) tf->Get(hname);
            if(!h){ std::cerr << "[warn] histogram '" << hname << "' not found in " << tf->GetName() << std::endl; tf->Close(); continue; }
            if(!sum){ sum = (TH1F*) h->Clone(outName); sum->SetDirectory(0); }
            else { sum->Add(h); }
            tf->Close();
        }
        if(!sum) { sum = new TH1F(outName, "", 1, 0, 1); sum->SetDirectory(0); sum->SetBinContent(1,0); sum->SetBinError(1,0); }
        return sum;
    };

    // Liste dei file MC (v2) e data
    std::vector<TString> files_DYG = {
      "DYGto2LG-1Jets_MLL-50_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_"+era+"_v2.root",
      "DYGto2LG-1Jets_MLL-50_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_"+era+"_v2.root",
      "DYGto2LG-1Jets_MLL-50_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_"+era+"_v2.root"
    };
    std::vector<TString> files_DY = {
      "DYto2L-4Jets_MLL-50to120_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root",
      "DYto2L-4Jets_MLL-50to120_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root",
      "DYto2L-4Jets_MLL-50to120_HT-100to400_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root",
      "DYto2L-4Jets_MLL-50to120_HT-400to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root",
      "DYto2L-4Jets_MLL-50to120_HT-800to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root",
      "DYto2L-4Jets_MLL-50to120_HT-1500to2500_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root",
      "DYto2L-4Jets_MLL-50to120_HT-2500_TuneCP5_13p6TeV_madgraphMLM-pythia8_"+era+"_v2.root"
    };
    TString data_file = "EGamma_"+era+"_v2.root";

    // Costruisci istogrammi
  TH1F* h_dyg = getSummedHisto(files_DYG, var, "h_dyg");
  TH1F* h_dy  = getSummedHisto(files_DY,  var, "h_dy");
    TFile *fdata = TFile::Open(dir+data_file);
    if(!fdata || fdata->IsZombie()) { std::cerr << "[error] cannot open data file " << (dir+data_file) << std::endl; if(fdata) fdata->Close(); continue; }
    TH1F *hdata  = (TH1F*)fdata->Get(var);
    if(!hdata){ std::cerr << "[error] histogram '" << var << "' not found in data file.\n"; fdata->Close(); continue; }
  hdata = (TH1F*)hdata->Clone("hdata"); hdata->SetDirectory(0); fdata->Close();

  // Rebin adattivo per migliorare la leggibilità
  auto doRebin = [&](TH1* h){ if(!h) return; };
  int r = 1;
  if(Var_Name[variable].Contains("_eta")) r = 5;            // 50 -> 10 bin
  else if(Var_Name[variable].Contains("_dr")) r = 10;       // 200 -> 20 bin
  else if(Var_Name[variable]=="dimu_mass") r = 2;          // 100 -> 50
  else if(Var_Name[variable]=="lead_pt" || Var_Name[variable]=="subl_pt") r = 2; // 20 -> 10
  else if(Var_Name[variable]=="photon_pt" || Var_Name[variable]=="dimu_pt" || Var_Name[variable]=="HT") r = 2; // 36 -> 18
  if(r>1){ h_dy->Rebin(r); h_dyg->Rebin(r); hdata->Rebin(r); }

  // Evita negativi in stack/log
  SanitizeForStack(h_dy);
  SanitizeForStack(h_dyg);

    /*
   double chi2_pre = 0;
    for (int bin = 1; bin <= zbb->GetNbinsX(); ++bin) {
        double content1 = zbb->GetBinContent(bin);
        double content2 = zcc->GetBinContent(bin);
        double error1 = zbb->GetBinError(bin);
        double error2 = zcc->GetBinError(bin);
        // Calculate chi-squared contribution for this bin
        double bin_chi2; 
        if ((error1 * error1 + error2 * error2) ==0 ) bin_chi2 = 0.0;
        else bin_chi2 = (content1 - content2) * (content1 - content2) / (error1 * error1 + error2 * error2);
        chi2_pre += bin_chi2/zbb->GetNbinsX();
    }

    double chi2_post = 0;
    for (int bin = 1; bin <= zbb->GetNbinsX(); ++bin) {
        double content1 = zbb->GetBinContent(bin);
        double content2 = zqq->GetBinContent(bin);
        double error1 = zbb->GetBinError(bin);
        double error2 = zqq->GetBinError(bin);
        // Calculate chi-squared contribution for this bin
        double bin_chi2;
        if ((error1 * error1 + error2 * error2) ==0 ) bin_chi2 = 0.0;
        else bin_chi2 = (content1 - content2) * (content1 - content2) / (error1 * error1 + error2 * error2);
        chi2_post += bin_chi2/zbb->GetNbinsX();
    }

   std::cout << chi2_pre << "   " << chi2_post << std::endl;
   zbb->GetXaxis()->SetTitle("#sump_{Z} (b1 + b2 + q1 + q2) [GeV]");
   */

gStyle->SetLineWidth(2);
TCanvas *canv = new TCanvas("canv", " ",370,119,800,800);
canv->Range(-25.95707,-4.739513,11.50268,-0.9363005);
canv->SetFillColor(0);
canv->SetBorderMode(0);
canv->SetBorderSize(3);
canv->SetGridy();
//canv->SetLogy();
canv->SetGridx();
canv->SetGridy();
canv->SetTickx(1);
canv->SetTicky(1);
canv->SetLeftMargin(0.1090258);
canv->SetRightMargin(0.04011461);
canv->SetFrameBorderMode(0);
canv->SetFrameLineWidth(3);
canv->SetFrameBorderMode(0);
canv->SetFrameLineWidth(3);
canv->SetFrameBorderMode(0);
canv->SetFrameLineWidth(3);
canv->SetFrameBorderMode(0);


      TPad *pad1 = new TPad("pad1", "pad1",0,0.3,1,1);
      pad1->Draw();
      pad1->Range(-148.318,2.314742,1071.865,8.438788);
      pad1->SetFillColor(0);
      pad1->SetBorderMode(0);
      pad1->SetBorderSize(2);
  // Usa scala log solo per variabili a larga dinamica
  Int_t setLog = (Var_Name[variable]=="photon_pt" || Var_Name[variable]=="dimu_pt" || Var_Name[variable]=="HT") ? 1 : 0;
  pad1->SetLogy(setLog);
      pad1->SetGridx();
      pad1->SetGridy();
      pad1->SetTickx(1);
      pad1->SetTicky(1);
      pad1->SetLeftMargin(0.1215539);
      pad1->SetRightMargin(0.05889724);
      pad1->SetBottomMargin(0.005321326);
      pad1->SetFrameBorderMode(0);
      pad1->SetFrameBorderMode(0);


      TPad *pad2 = new TPad("pad2", "newpad",0,0,1,0.3);
      pad2->Draw();
      pad2->Range(-148.318,-0.4344828,1071.865,1.009655);
      pad2->SetGridx();
      pad2->SetGridy();
      pad2->SetTickx(1);
      pad2->SetTicky(1);
      pad2->SetFillColor(0);
      pad2->SetBorderMode(0);
      pad2->SetBorderSize(2);
      pad2->SetLeftMargin(0.1215539);
      pad2->SetRightMargin(0.06015038);
      pad2->SetTopMargin(0.006685769);
      pad2->SetBottomMargin(0.3008596);
      pad2->SetFrameBorderMode(0);
      pad2->SetFrameBorderMode(0);

      pad1->cd();

  hdata->SetStats(0);
  hdata->SetTitle("");
  hdata->SetMarkerStyle(20);
  hdata->SetMarkerColor(kBlack);
  hdata->SetLineColor(kBlack);
  hdata->SetLineWidth(2);
  hdata->SetMarkerSize(1.2);
  hdata->GetYaxis()->SetTitle("Events");
  // Titoli degli assi per variabile
  auto xTitle = [&](const TString& v){
    if(v=="photon_pt") return TString("p_{T}(#gamma) [GeV]");
    if(v=="photon_eta") return TString("#eta(#gamma)");
    if(v=="dimu_mass") return TString("m(#mu#mu) [GeV]");
    if(v=="dimu_pt") return TString("p_{T}(#mu#mu) [GeV]");
    if(v=="dimu_eta") return TString("#eta(#mu#mu)");
    if(v=="dimu_dr") return TString("#DeltaR(#mu,#mu)");
    if(v=="leadmu_photon_dr") return TString("#DeltaR(#mu_{lead},#gamma)");
    if(v=="sublmu_photon_dr") return TString("#DeltaR(#mu_{sublead},#gamma)");
    if(v=="lead_pt") return TString("p_{T}(#mu_{lead}) [GeV]");
    if(v=="subl_pt") return TString("p_{T}(#mu_{sublead}) [GeV]");
    if(v=="lead_eta") return TString("#eta(#mu_{lead})");
    if(v=="subl_eta") return TString("#eta(#mu_{sublead})");
    if(v=="HT") return TString("H_{T} [GeV]");
    return Var_Name[variable];
  };
  hdata->GetXaxis()->SetTitle(xTitle(Var_Name[variable]));
  hdata->GetXaxis()->SetNdivisions(506);
  hdata->GetXaxis()->SetLabelFont(42);
  hdata->GetXaxis()->SetLabelSize(0.050);
  hdata->GetXaxis()->SetTitleSize(0.05);
  hdata->GetXaxis()->SetTickLength(0.025);
  hdata->GetXaxis()->SetTitleFont(42);
  hdata->GetXaxis()->SetTitleOffset(1.00);
  hdata->GetYaxis()->SetNdivisions(506);
  hdata->GetYaxis()->SetLabelFont(42);
  hdata->GetYaxis()->SetLabelSize(0.050);
  hdata->GetYaxis()->SetTitleSize(0.05);
  hdata->GetYaxis()->SetTitleFont(42);
  hdata->GetYaxis()->SetTitleOffset(1.20);
  if(setLog){
    hdata->GetYaxis()->SetRangeUser(std::max(1e-3, 0.1*hdata->GetMaximum()/100.0), 100.0*hdata->GetMaximum());
  } else {
    hdata->GetYaxis()->SetRangeUser(0.0, 1.5*hdata->GetMaximum());
  }
  hdata->Draw("EPX0");

  // Stile MC
  h_dy->SetLineWidth(2);
  h_dy->SetLineColor(TColor::GetColor("#ffa90e"));
  h_dy->SetFillColor(TColor::GetColor("#ffa90e"));

  h_dyg->SetLineWidth(2);
  h_dyg->SetLineColor(colorIndices[2]);
  h_dyg->SetFillColor(colorIndices[2]);


  THStack *hs = new THStack("hs","");
                // Ordine: DY (fondo) + DYG (Z+#gamma)
                hs->Add(h_dy,"hist");
                hs->Add(h_dyg,"hist");
                hs->Draw("hist,SAME");
                hdata->Draw("EPX0,SAME");

                   double xmax = hdata->GetXaxis()->GetXmax();
                   double xmin = hdata->GetXaxis()->GetXmin();

  TH1F *totalMC = new TH1F();
  totalMC->Sumw2();
  totalMC = (TH1F*)h_dy->Clone();
  totalMC->Add(h_dyg);
  canv->SetBorderSize(2);
  double scaletodata = (totalMC->Integral()>0) ? hdata->Integral()/totalMC->Integral() : 1.0;
  if(normToData){ h_dy->Scale(scaletodata); h_dyg->Scale(scaletodata); totalMC->Scale(scaletodata); }

  // Create graph for uncertainties
    int nBins = totalMC->GetNbinsX();
    TGraphAsymmErrors* gr_unc = new TGraphAsymmErrors(nBins);

    for (int i = 1; i <= nBins; i++) {
        double binCenter = totalMC->GetBinCenter(i);
        double binContent = totalMC->GetBinContent(i);
        double binError = totalMC->GetBinError(i); // Stat uncertainty

        gr_unc->SetPoint(i - 1, binCenter, binContent);
        gr_unc->SetPointError(i - 1, totalMC->GetBinWidth(i) / 2, totalMC->GetBinWidth(i) / 2, binError, binError);
    }

    // Style the uncertainty band
    gr_unc->SetFillColorAlpha(kBlack, 0.35); // Semi-transparent black
    gr_unc->SetFillStyle(3004); // Hatched pattern
    gr_unc->SetLineWidth(0); // No line
    gr_unc->Draw("E2,SAME");
    hdata->Draw("EPX0,SAME");
  gSystem->Sleep(1);
TLegend *legend1 = new TLegend(0.50, 0.6, 0.85, 0.85);
legend1->SetTextFont(42);
//legend1->SetTextFont(13);
legend1->SetLineColor(0);
legend1->SetTextSize(0.055);
legend1->SetFillStyle(0);
legend1->AddEntry(hdata, "Data", "epx0");
legend1->AddEntry(h_dy, "DY #rightarrow #mu#mu (jets/fake #gamma)", "F");
legend1->AddEntry(h_dyg, "DY+#gamma #rightarrow #mu#mu#gamma", "F");
legend1->AddEntry(gr_unc, "MC uncertainties", "F");
legend1->Draw();

canv->cd();
pad2->cd();
TH1F *h_A = (TH1F*)  hdata->Clone();
h_A->Divide(totalMC);

   h_A->SetStats(0);
   h_A->SetTitle("");
   h_A->SetMarkerColor(kBlack);
   h_A->SetLineColor(kBlack);
   h_A->SetLineWidth(2);
  h_A->GetXaxis()->SetTitle(xTitle(Var_Name[variable]));
   h_A->GetXaxis()->SetLabelFont(42);
   h_A->GetXaxis()->SetLabelSize(0.09);
   h_A->GetXaxis()->SetTitleSize(0.11);
   h_A->GetXaxis()->SetTickLength(0.08);
   h_A->GetXaxis()->SetTitleOffset(1.05);
   h_A->GetXaxis()->SetTitleFont(42);
   h_A->GetYaxis()->SetTitle("Data/MC");
   h_A->GetYaxis()->SetNdivisions(207);
   h_A->GetYaxis()->SetLabelFont(42);
   h_A->GetYaxis()->SetLabelSize(0.09);
   h_A->GetYaxis()->SetTitleSize(0.11);
   h_A->GetYaxis()->SetTickLength(0.02);
   h_A->GetYaxis()->SetTitleOffset(0.42);
   h_A->GetYaxis()->SetTitleFont(42);
   h_A->GetXaxis()->SetRangeUser(xmin,xmax);

   h_A->GetYaxis()->SetRangeUser(0.0,2.0);
   h_A->SetMarkerStyle(20);
   h_A->SetMarkerColor(kBlack);
   h_A->SetLineColor(kBlack);
   h_A->SetMarkerSize(1.2);
   h_A->Draw("EP");

  // Banda di incertezza sul rapporto: centrata a 1 con errori relativi del MC
  int nBinsRatio = totalMC->GetNbinsX();
  TGraphAsymmErrors *gr_ratio = new TGraphAsymmErrors(nBinsRatio);
  for (int i=1; i<=nBinsRatio; ++i){
      double x  = totalMC->GetBinCenter(i);
      double ex = totalMC->GetBinWidth(i)/2.0;
      double mc = totalMC->GetBinContent(i);
      double emc= totalMC->GetBinError(i);
      double rel = (mc>0) ? emc/mc : 0.0;
      gr_ratio->SetPoint(i-1, x, 1.0);
      gr_ratio->SetPointError(i-1, ex, ex, rel, rel);
  }
  gr_ratio->SetFillColor(kGray+1);
  gr_ratio->SetLineColor(kGray+1);
  gr_ratio->SetFillStyle(3144);
  gr_ratio->Draw ("e2 same");
  h_A->Draw("EP,SAME");

  canv->Update();
  writeExtraText = true;
  extraText  = "Preliminary";
  // Evita duplicazione dell'energia; CMS_lumi aggiunge già (13 TeV)/(13.6 TeV) in base a iPeriod
  custom_lumi_13TeV = "27.0 fb^{-1}";
  CMS_lumi(canv, 4, 11);

  // plotting dir
  canv->Print(plot_dir+Var_Name[variable]+"_"+era+".pdf");
  canv->Print(plot_dir+Var_Name[variable]+"_"+era+".png");
  }
}



