#include "TTree.h"
#include "TBrowser.h"
#include "TH1D.h"
#include "TSystem.h"
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


void plotter_categorization()
{
   InitializeColors();
   TString era = "2022";
   Float_t lumi;
   if(era == "2023preEE") lumi = 8.0;
   else if(era == "2023") lumi = 27.0;

   TString CAT[] = {"CAT3", "CAT2","CAT1","CAT0"};

   for (int jj = 0; jj < (int)(sizeof(CAT) / sizeof(CAT[0])); jj++) {

   // plotting dir
   TString plot_dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Plots/zbbg/GParT/"+CAT[jj]+"/";

   //TString dir = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/HH4b/Zbb_calibration/NTuples_Zmumu/zbbg/Histogram/GParT/"+cat+"/";
   TString dir = "/gwpool/users/fballo/HiggsAnalysis/code/boostedZgamma/Ntuples_zbbg/Histogram/GParT/"+CAT[jj]+"/";

   //TString Var_Name[] = {"dimu_dr","photon_eta","photon_pt","dimu_mass","dimu_pt","dimu_eta","lead_pt","lead_eta","subl_pt","subl_eta","lead_eta","subl_eta","HT"};
   //TString Var_Name[] = {"ak8_pt","dimu_mass","dimu_pt","dimu_dr"};
   TString Var_Name[] = {"photon_pt","photon_eta","ak8_pnet","ak8_gpart","ak8_mass","ak8_mSD","ak8_gpart_mass","ak8_pT","ak8_eta","photon_jet_dr","HT"};
   

   for(int variable = 0; variable <11; variable++)
   {
    std::cout << Var_Name[variable] << std::endl;
    TString var = "hist_" + Var_Name[variable];
    // ensure plot directory exists
    TString outdir = plot_dir;
    if(!outdir.EndsWith("/")) outdir += "/";
    TString dir_to_create = outdir;
    if(dir_to_create.EndsWith("/")) dir_to_create.Chop();
    gSystem->mkdir(dir_to_create.Data(), kTRUE);

    TFile *dy = TFile::Open(dir+"ZGto2QG_2022.root");
    TFile *vv = TFile::Open(dir+"Zto2Q_2022.root");
    TFile *gj = TFile::Open(dir+"GJ_2022.root");
    TFile *fdata = TFile::Open(dir+"EGamma_2022.root");

    if(!dy || dy->IsZombie() || !vv || vv->IsZombie() || !gj || gj->IsZombie() || !fdata || fdata->IsZombie()) {
      std::cerr << "Error: one or more input files missing for variable '" << Var_Name[variable] << "' in dir '" << dir << "'\n";
      if(dy) { dy->Close(); }
      if(vv) { vv->Close(); }
      if(gj) { gj->Close(); }
      if(fdata) { fdata->Close(); }
      continue;
    }


    TH1F *hdy = (TH1F*)dy->Get(var);
    TH1F *hgj = (TH1F*)gj->Get(var);
    TH1F *hvv = (TH1F*)vv->Get(var);
    TH1F *hdata  = (TH1F*)fdata->Get(var);

    if(!hdy || !hgj || !hvv || !hdata) {
      std::cerr << "Error: missing histogram '" << var << "' in one of the input files; skipping.\n";
      // close files and continue
      dy->Close(); vv->Close(); gj->Close(); fdata->Close();
      continue;
    }
    double binEdges[9] = {300,350,400,450,500,600,800,1000,1500};
    int nBins_ = sizeof(binEdges) / sizeof(binEdges[0]) - 1;
    if(Var_Name[variable] == "ak8_pnet" || Var_Name[variable] == "ak8_gpart") 
    {
      hdy->Rebin(1); 
      hvv->Rebin(1);
      hgj->Rebin(1);
      hdata->Rebin(1); 
    }
    if(Var_Name[variable] == "ak8_mass" || Var_Name[variable] == "ak8_mSD" ||  Var_Name[variable] == "ak8_gpart_mass")
    {
      hdy->Rebin(5);
      hvv->Rebin(5);
      hgj->Rebin(5);
      hdata->Rebin(5);
    }
    else if(Var_Name[variable] == "photon_jet_dr")
    {
      hdy->Rebin(10);
      hvv->Rebin(10);
      hgj->Rebin(10);
      hdata->Rebin(10);
    }
    if(Var_Name[variable] == "dimu_pt")
    {
       hdy= (TH1F*) hdy->Rebin(nBins_, "hdy", binEdges);
       hvv= (TH1F*) hvv->Rebin(nBins_, "hvv", binEdges);
       hdata= (TH1F*) hdata->Rebin(nBins_, "hdata", binEdges); 
    } 
    /*
    AddOverflowToLastBin(zbb);
    AddOverflowToLastBin(zcc);
    AddOverflowToLastBin(zqq);
    AddOverflowToLastBin(wqq);
    AddOverflowToLastBin(qcd);
    AddOverflowToLastBin(hdata);
    */
    //Float_t qcd_fraction = (hdata->Integral() - zbb->Integral() - zcc->Integral() - zqq->Integral() - wqq->Integral())/qcd->Integral(); 
    //qcd->Scale(qcd_fraction);
    std::cout << hdata->Integral() << "  " << hgj->Integral() <<  std::endl;
    double gjInt = hgj->Integral();
    if(gjInt == 0) {
      std::cerr << "Warning: GJ integral is zero for '"<<Var_Name[variable]<<"' - skipping scale.\n";
      dy->Close(); vv->Close(); gj->Close(); fdata->Close();
      continue;
    }
    Float_t qcd_fraction = (hdata->Integral() - hdy->Integral() - hvv->Integral() )/ gjInt; 
    std::cout << qcd_fraction << std::endl;
    hgj->Scale(qcd_fraction);
    //

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
      if(Var_Name[variable]=="gpart" || Var_Name[variable]=="pnet") pad1->SetLogy();
      pad1->SetLogy();
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
  hdata->GetXaxis()->SetTitle("D_{VBF} ");
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
  hdata->GetYaxis()->SetRangeUser(0.001,100.0*hdata->GetMaximum());
  //hdata->GetYaxis()->SetRangeUser(0.0,2.0*hdata->GetMaximum());
  if(Var_Name[variable]=="gpart" || Var_Name[variable]=="pnet") hdata->GetYaxis()->SetRangeUser(0.001,1000.0*hdata->GetMaximum());
  hdata->Draw("EPX0");

  hvv->SetLineWidth(2);
  hvv->SetLineColor(TColor::GetColor("#ffa90e"));
  hvv->SetFillColor(TColor::GetColor("#ffa90e"));

  hdy->SetLineWidth(2);
  hdy->SetLineColor(colorIndices[2]);
  hdy->SetFillColor(colorIndices[2]);

  hgj->SetLineWidth(2);
  hgj->SetLineColor(TColor::GetColor("#3f90da"));
  hgj->SetFillColor(TColor::GetColor("#3f90da"));


  THStack *hs = new THStack("hs",",");

                hs->Add(hgj,"hist"); 
                hs->Add(hvv,"hist");
                hs->Add(hdy,"hist");

                hs->Draw("hist");
                hdata->Draw("EPX0 same");

                double xmax = hdata->GetXaxis()->GetXmax();
                double xmin = hdata->GetXaxis()->GetXmin();

  TH1F *totalMC = (TH1F*)hdy->Clone();
  totalMC->Sumw2();
  totalMC->Add(hvv);
  totalMC->Add(hgj);
  canv->SetBorderSize(2);
  double totalInt = totalMC->Integral();
  double scaletodata = 1.0;
  if(totalInt != 0) scaletodata = hdata->Integral()/totalInt;
  else std::cerr << "Warning: totalMC integral is zero for '"<<Var_Name[variable]<<"' - scaletodata set to 1.\n";

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
    gr_unc->Draw("E2 same");
  hdata->Draw("EPX0 same");
  gSystem->Sleep(1);
TLegend *legend1 = new TLegend(0.50, 0.6, 0.85, 0.85);
legend1->SetTextFont(42);
//legend1->SetTextFont(13);
legend1->SetLineColor(0);
legend1->SetTextSize(0.055);
legend1->SetFillStyle(0);
legend1->AddEntry(hdata, "Data", "epx0");
// legend: match stack order (gamma+jets, Z->bb, Z->bb+gamma)
legend1->AddEntry(hgj, "#gamma+jets", "F");
legend1->AddEntry(hvv, "Z #rightarrow bb", "F");
legend1->AddEntry(hdy, "Z #rightarrow bb + #gamma", "F");
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
   h_A->GetXaxis()->SetTitle(Var_Name[variable]);
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

  TGraphAsymmErrors *staterree = (TGraphAsymmErrors*)staterror(totalMC,totalMC);
  staterree->SetFillColor(kGray+1);
  staterree->SetLineColor(kGray+1);
  staterree->SetFillStyle(3144);
  staterree->Draw ("e2 same");
  h_A->Draw("EP same");

  canv->Update();
  writeExtraText = true;
  extraText  = "Preliminary";
  custom_lumi_13TeV = "35.0 fb^{-1},";
  CMS_lumi(canv, 4, 11);

  // save both linear and log-scale versions with tailored Y ranges
  // compute dynamic min/max
  double maxData = hdata->GetMaximum();
  double maxMCplot = hs->GetMaximum();
  double ymax = std::max(maxData, maxMCplot);
  double minPos = 1e9;
  for(int ib=1; ib<=totalMC->GetNbinsX(); ++ib){
    double v1 = totalMC->GetBinContent(ib);
    double v2 = hdata->GetBinContent(ib);
    if(v1>0 && v1<minPos) minPos = v1;
    if(v2>0 && v2<minPos) minPos = v2;
  }
  if(!(minPos>0 && minPos<1e9)) minPos = 1e-3;

  // Linear
  pad1->SetLogy(0);
  hdata->GetYaxis()->SetRangeUser(0.0, 1.3*ymax);
  pad1->Modified();
  pad1->Update();
  canv->Modified();
  canv->Update();
  canv->Print(plot_dir+Var_Name[variable]+"_"+era+"_lin.pdf");
  canv->Print(plot_dir+Var_Name[variable]+"_"+era+"_lin.png");

  // Log
  pad1->SetLogy(1);
  hdata->GetYaxis()->SetRangeUser(0.5*minPos, 10.0*ymax);
  pad1->Modified();
  pad1->Update();
  canv->Modified();
  canv->Update();
  canv->Print(plot_dir+Var_Name[variable]+"_"+era+"_log.pdf");
  canv->Print(plot_dir+Var_Name[variable]+"_"+era+"_log.png");
  }

}

}



