import os
import sys
import glob
import argparse
import subprocess
import shutil
import time
import ROOT
import numpy as np
from array import array 

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-dir', type=str, default='', help='input directory with files')
parser.add_argument('-c', '--cat-id', type=str, default='', choices=['cat0','cat1','cat2','cat3','cat4'], help='category identifier')
parser.add_argument('-o', '--output-directory', type=str, default='', help='name of the output directory')
parser.add_argument('-t', '--tagger', type=str, default='', choices=['pnet','gpart'], help='tagger to be considered')
parser.add_argument('-s', '--sum-matched-unmatched', action='store_true', help='sum matched and umatched. Default is unmatched go in gj backgorund')
parser.add_argument('-y', '--year', type=str, default='2022', choices=['2022','2023'], help='data taking period')
parser.add_argument('--mass-obs', type=str, default='mSD', help='mass observable to fit')
parser.add_argument('--rebin-pdf', type=int, default=10, help='rebin factor for pdf')
parser.add_argument('--rebin-data', type=int, default=2, help='rebin factor for data')
parser.add_argument('--jms-uncertainty', type=float, default=0.03, help='value of the jet mass scale uncertainty to be used')
parser.add_argument('--jmr-uncertainty', type=float, default=0.10, help='value of the jet mass resolution uncertainty to be used')
parser.add_argument('--float-mass-peak', action='store_true', help='float peak position')
parser.add_argument('--float-mass-reso', action='store_true', help='float peak resolution')
parser.add_argument('--use-poly-bkg', action='store_true', help='use polynomial background for pass cat0 and cat1 instead of a gaussian')
parser.add_argument('--plot-bkg', action='store_true', help='plot the bkg pdf')
args = parser.parse_args()

lumi_dict = {
    '2022': 35.,
    '2023': 27.7
}
    
ROOT.gInterpreter.ProcessLine('#include "CMS_style.h"')
ROOT.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
os.system("mkdir -p "+args.output_directory)

def createCardTemplate(cat_score,f_ws,sys={},ws_name="w"):

    dc_name = "datacard_zbb_"+cat_score+".txt"

    f = open(dc_name,"w")
    f.write("imax * \n")
    f.write("jmax * \n")
    f.write("kmax * \n")
    f.write("------------ \n")
    f.write("shapes zjet "+f_ws+" "+ws_name+":"+"pdf_zjet\n")
    f.write("shapes wjet "+f_ws+" "+ws_name+":"+"pdf_wjet\n")
    f.write("shapes zg "+f_ws+" "+ws_name+":"+"pdf_zg\n")
    f.write("shapes gj "+f_ws+" "+ws_name+":"+"pdf_gj\n")
    f.write("shapes data_obs "+f_ws+" "+ws_name+":"+"data\n")
    f.write("------------ \n")
    f.write("bin   \n")
    f.write("observation -1 \n")        
    f.write("------------ \n")
    f.write("process   zjet  wjet  zg  gj \n")
    f.write("process   0     1     2   3 \n")
    f.write("rate   1     1     1   1 \n")
    f.write("------------ \n")
    f.write("lumi_13p6TeV  lnN  1.014  1.014  1.014  -\n")
    if "gj_scale_zjet" in sys:
        f.write("gj_scale_zjet  lnN  %.2f  -  -  -\n"%(sys["gj_scale_zjet"]))
    else:
        ## uncertainty on gj taken from NLO Z+jets sample in genXsecAnalyzer
        f.write("gj_scale_zjet  lnN  0.94/1.05  -  -  -\n")
    if "gj_scale_wjet" in sys:
        f.write("gj_scale_wjet  lnN  -  %.2f  -  -\n"%(sys["gj_scale_wjet"]))
    else:
        ## uncertainty on gj taken from NLO W+jets sample in genXsecAnalyzer
        f.write("gj_scale_wjet  lnN  -  0.95/1.05  -  -\n")
    if "gj_scale_zg" in sys:
        f.write("gj_scale_zg  lnN  -  -  %.2f  -\n"%(sys["gj_scale_zg"]))
    else:
        ## uncertainty on gj taken from NLO Zγ sample in genXsecAnalyzer
        f.write("gj_scale_zg  lnN  -  -  0.94/1.05  -\n")
    if "NLOEWK_corr_zjet" in sys:
        f.write("NLOEWK_corr_zjet  lnN  %.2f  -  -  -\n"%(sys["NLOEWK_corr_zjet"]))
    else:
        f.write("NLOEWK_corr_zjet  lnN  0.92/1.09  -  -  -\n")
    if "NLOEWK_corr_wjet" in sys:
        f.write("NLOEWK_corr_wjet  lnN  -  %.2f  -  -\n"%(sys["NLOEWK_corr_wjet"]))
    else:        
        f.write("NLOEWK_corr_wjet  lnN  -  0.92/1.09  -  -\n")
    if "NLOEWK_corr_zg" in sys:
        f.write("NLOEWK_corr_zg  lnN  -  -  %.2f  -\n"%(sys["NLOEWK_corr_zg"]))
    else:        
        f.write("NLOEWK_corr_zg  lnN  -  -  0.92/1.09  -\n")
    if "pdf_zjet" in sys:
        f.write("pdf_zjet lnN  %.2f  -  -  -\n"%(sys["pdf_zjet"]))
    else:
        f.write("pdf_zjet lnN  0.98/1.02  -  -  -\n")
    if "pdf_wjet" in sys:
        f.write("pdf_wjet lnN  -  %.2f  -  -\n"%(sys["pdf_wjet"]))
    else:
        f.write("pdf_wjet lnN  -  0.98/1.02  -  -\n")
    if "pdf_zg" in sys:
        f.write("pdf_zg lnN  -  -  %.2f  -\n"%(sys["pdf_zg"]))
    else:
        f.write("pdf_zg lnN  -  -  0.98/1.02  -\n")
    f.write("CMS_scale_j  lnN  0.97/1.03 0.96/1.04  0.97/1.03  -\n")
    f.write("CMS_res_j  lnN  0.98/1.01   0.97/1.02  0.98/1.01  -\n")
    f.write("CMS_trigger_zbb  lnN  0.98/1.02  0.97/1.02  0.98/1.02  -\n")
    if args.float_mass_peak:
        f.write("CMS_jms_unc flatParam\n")
    else:
        f.write("CMS_jms_unc param 0 1\n")
    if args.float_mass_reso:
        f.write("CMS_jmr_unc flatParam\n")
    else:
        f.write("CMS_jmr_unc param 0 1\n")
        
    f.close()
    return dc_name

# def createCardTemplate(cat_score,cat_sel,f_ws,sys={},ws_name="w"):

#     dc_name = "datacard_zbb_"+cat_score+"_"+cat_sel+".txt"

#     f = open(dc_name,"w")
#     f.write("imax * \n")
#     f.write("jmax * \n")
#     f.write("kmax * \n")
#     f.write("------------ \n")
#     f.write("shapes zjet "+cat_sel+" "+f_ws+" "+ws_name+":"+"pdf_zjet_"+cat_sel+"\n")
#     f.write("shapes wjet "+cat_sel+" "+f_ws+" "+ws_name+":"+"pdf_wjet_"+cat_sel+"\n")
#     f.write("shapes qcd "+cat_sel+" "+f_ws+" "+ws_name+":"+"pdf_qcd_"+cat_sel+"\n")
#     f.write("shapes data_obs "+cat_sel+" "+f_ws+" "+ws_name+":"+"data_"+cat_sel+"\n")
#     f.write("------------ \n")
#     f.write("bin   "+cat_sel+" \n")
#     f.write("observation -1 \n")        
#     f.write("------------ \n")
#     f.write("bin   "+cat_sel+" "+cat_sel+" "+cat_sel+" \n")
#     f.write("process   zjet  wjet  qcd \n")
#     f.write("process   0     1     2 \n")
#     f.write("rate   1     1     1 \n")
#     f.write("------------ \n")
#     f.write("lumi_13p6TeV  lnN  1.014  1.014  -\n")
#     if "QCD_scale_zjet" in sys:
#         f.write("QCD_scale_zjet  lnN  %.2f  -  -\n"%(sys["QCD_scale_zjet"]))
#     else:
#         ## uncertainty on QCD taken from NLO Z+jets sample in genXsecAnalyzer
#         f.write("QCD_scale_zjet  lnN  0.94/1.05  -  -\n")
#     if "QCD_scale_wjet" in sys:
#         f.write("QCD_scale_wjet  lnN  -  %.2f  -\n"%(sys["QCD_scale_wjet"]))
#     else:
#         ## uncertainty on QCD taken from NLO W+jets sample in genXsecAnalyzer
#         f.write("QCD_scale_wjet  lnN  -  0.95/1.05  -\n")
#     if "NLOEWK_corr_zjet" in sys:
#         f.write("NLOEWK_corr_zjet  lnN  %.2f  -  -\n"%(sys["NLOEWK_corr_zjet"]))
#     else:
#         f.write("NLOEWK_corr_zjet  lnN  0.92/1.09  -  -\n")
#     if "NLOEWK_corr_wjet" in sys:
#         f.write("NLOEWK_corr_wjet  lnN  -  %.2f  -\n"%(sys["NLOEWK_corr_zjet"]))
#     else:        
#         f.write("NLOEWK_corr_wjet  lnN  -  0.92/1.09  -\n")
#     if "pdf_zjet" in sys:
#         f.write("pdf_zjet lnN  %.2f  -  -\n"%(sys["pdf_zjet"]))
#     else:
#         f.write("pdf_zjet lnN  0.98/1.02  -  -\n")
#     if "pdf_wjet" in sys:
#         f.write("pdf_wjet lnN  -  %.2f  -\n"%(sys["pdf_wjet"]))
#     else:
#         f.write("pdf_wjet lnN  -  0.98/1.02  -\n")
#     f.write("CMS_scale_j  lnN  0.97/1.03 0.96/1.04  -\n")
#     f.write("CMS_res_j  lnN  0.98/1.01   0.97/1.02  -\n")
#     f.write("CMS_trigger_zbb  lnN  0.98/1.02  0.97/1.02  -\n")
#     if args.float_mass_peak:
#         f.write("CMS_jms_unc flatParam\n")
#     else:
#         f.write("CMS_jms_unc param 0 1\n")
#     if args.float_mass_reso:
#         f.write("CMS_jmr_unc flatParam\n")
#     else:
#         f.write("CMS_jmr_unc param 0 1\n")
        
#     f.close()
#     return dc_name

def getChi2(h,h_fit):
    chi2 = 0
    ndf  = 0
    for i in range(0,h.GetNbinsX()):
        if h.GetBinContent(i+1) == 0 : continue
        res = h.GetBinContent(i+1)-h_fit.GetBinContent(i+1)
        chi2 += (res**2)/(h.GetBinError(i+1)**2)
        ndf = ndf + 1
    return chi2,ndf

def getMaxAndHWHM(h):
    binMax = h.GetMaximumBin()
    xMax = h.GetXaxis().GetBinCenter(binMax)
    yMax = h.GetBinContent(binMax)    
    halfMax = yMax / 2.0
    binLeft = binMax
    while binLeft > 1 and h.GetBinContent(binLeft) > halfMax:
        binLeft = binLeft-1
    binRight = binMax
    while binRight > 1 and h.GetBinContent(binRight) > halfMax:
        binRight = binRight+1

    xLeft=h.GetXaxis().GetBinCenter(binLeft)
    xRight=h.GetXaxis().GetBinCenter(binRight)
    hwhm = (xRight-xLeft)/2.

    return xMax,hwhm

def sum_histograms_from_files(file_list, hist_name):
    """Open each file in file_list, extract hist_name (or nearest mass hist), and return summed TH1 (clone).
    Returns None if no histogram found.
    """
    h_sum = None
    for fn in file_list:
        try:
            f = ROOT.TFile(fn, 'READ')
        except Exception as e:
            print('Warning: cannot open', fn, e)
            continue
        if not f or not f.IsOpen():
            print('Warning: cannot open', fn)
            continue
        # try the exact name first
        h = None
        if f.GetListOfKeys().Contains(hist_name):
            h = f.Get(hist_name)
        else:
            # try common alternatives containing 'mass' or the root of hist_name
            for k in f.GetListOfKeys():
                kn = k.GetName()
                if hist_name in kn:
                    h = f.Get(kn)
                    break
                if 'mass' in kn.lower() and 'hist' in kn.lower():
                    h = f.Get(kn)
                    break
        if not h:
            # fallback: pick any TH1 in the file
            for k in f.GetListOfKeys():
                obj = f.Get(k.GetName())
                if obj.InheritsFrom('TH1'):
                    h = obj
                    break
        if not h:
            f.Close()
            continue

        if h_sum is None:
            h_sum = h.Clone(os.path.basename(fn) + '_' + h.GetName())
            h_sum.SetDirectory(0)
        else:
            try:
                h_sum.Add(h)
            except Exception as e:
                print('Warning: failed to add histogram from', fn, e)
        f.Close()
    return h_sum

cat_candidates = [
    args.input_dir,
    os.path.join(args.input_dir, args.cat_id),
    os.path.join(args.input_dir, args.cat_id.upper()),
    os.path.join(os.path.dirname(args.input_dir), args.cat_id),
    os.path.join(os.path.dirname(args.input_dir), args.cat_id.upper())
]
cat_dir = None
for p in cat_candidates:
    if p and os.path.isdir(p):
        if any(f.endswith('.root') for f in os.listdir(p)):
            cat_dir = p
            break

if cat_dir:
    print('Detected CAT directory:', cat_dir)
    files = sorted(glob.glob(os.path.join(cat_dir, '*.root')))
    year = args.year
    preferred = {
        'z': [f'Zto2Q_{year}.root', 'Zto2Q_2022.root', 'Zto2Q_2023.root'],
        'w': [f'WGto2QG_{year}.root', 'Wto2Q_'+year+'.root', 'WGto2QG_2022.root'],
        'gj': [f'GJ_{year}.root', 'GJ_2022.root', 'GJ_2023.root'],
        'zg': [f'ZGto2QG_{year}.root', 'ZGto2QG_2022.root', 'ZGto2QG_2023.root'],
        'data': [f'EGamma_{year}.root', f'JetMET_{year}.root']
    }

    z_files = []
    w_files = []
    gj_files = []
    zg_files = []
    data_files = []

    for name in preferred['z']:
        p = os.path.join(cat_dir, name)
        if os.path.exists(p):
            z_files.append(p)
    for name in preferred['w']:
        p = os.path.join(cat_dir, name)
        if os.path.exists(p):
            w_files.append(p)
    for name in preferred['gj']:
        p = os.path.join(cat_dir, name)
        if os.path.exists(p):
            gj_files.append(p)
    for name in preferred['data']:
        p = os.path.join(cat_dir, name)
        if os.path.exists(p):
            data_files.append(p)
    for name in preferred['zg']:
        p = os.path.join(cat_dir, name)
        if os.path.exists(p):
            zg_files.append(p)

    z_files += [f for f in files if os.path.basename(f).lower().startswith('zto2q-4jets')]

    if not z_files:
        z_files = [f for f in files if 'zto2q' in os.path.basename(f).lower() or 'zgto2qg' in os.path.basename(f).lower()]
    if not w_files:
        w_files = [f for f in files if 'wgto2qg' in os.path.basename(f).lower() or os.path.basename(f).lower().startswith('w')]
    if not gj_files:
        gj_files = [f for f in files if os.path.basename(f).lower().startswith('gj') or 'gjet' in os.path.basename(f).lower()]
    if not zg_files:
        zg_files = [f for f in files if os.path.basename(f).lower().startswith('zg') or 'zgjet' in os.path.basename(f).lower()]
    if not data_files:
        data_files = [f for f in files if 'egamma' in os.path.basename(f).lower() or 'jetmet' in os.path.basename(f).lower()]

    hist_name = 'hist_ak8_' + args.mass_obs

    h_zjet = sum_histograms_from_files(z_files, hist_name)
    h_wjet = sum_histograms_from_files(w_files, hist_name)
    h_gj = sum_histograms_from_files(gj_files, hist_name)
    h_zg = sum_histograms_from_files(zg_files, hist_name)
    
    h_data = sum_histograms_from_files(data_files, hist_name)

    f_zjet = None
    f_wjet = None
    f_gj = None
    f_zg = None
    f_data = None
else:
    if args.year == "2022":
        pref_z = os.path.join(args.input_dir, 'Zto2Q_2022.root')
        pref_w = os.path.join(args.input_dir, 'WGto2QG_2022.root')
        pref_gj = os.path.join(args.input_dir, 'GJ_2022.root')
        pref_zg = os.path.join(args.input_dir, 'ZGto2QG_2022.root')
        pref_data = os.path.join(args.input_dir, 'EGamma_2022.root')

        if os.path.exists(pref_z):
            f_zjet = ROOT.TFile(pref_z, 'READ')
        else:
            f_zjet = ROOT.TFile(os.path.join(args.input_dir, 'Zto2Q_2022.root'), 'READ')
        if os.path.exists(pref_w):
            f_wjet = ROOT.TFile(pref_w, 'READ')
        else:
            f_wjet = ROOT.TFile(os.path.join(args.input_dir, 'Wto2Q_2022.root'), 'READ')
        if os.path.exists(pref_gj):
            f_gj = ROOT.TFile(pref_gj, 'READ')
        else:
            f_gj = ROOT.TFile(os.path.join(args.input_dir, 'gj_2022.root'), 'READ')
        if os.path.exists(pref_zg):
            f_zg = ROOT.TFile(pref_zg, 'READ')
        else:
            f_zg = ROOT.TFile(os.path.join(args.input_dir, 'ZGto2QG_2022.root'), 'READ')
        if os.path.exists(pref_data):
            f_data = ROOT.TFile(pref_data, 'READ')
        else:
            f_data = ROOT.TFile(os.path.join(args.input_dir, 'JetMET_2022.root'), 'READ')
    elif args.year == "2023":
        pref_z = os.path.join(args.input_dir, 'Zto2Q_2023.root')
        pref_w = os.path.join(args.input_dir, 'WGto2QG_2023.root')
        pref_gj = os.path.join(args.input_dir, 'GJ_2023.root')
        pref_zg = os.path.join(args.input_dir, 'ZGto2QG_2023.root')
        pref_data = os.path.join(args.input_dir, 'EGamma_2023.root')
        if os.path.exists(pref_z):
            f_zjet = ROOT.TFile(pref_z, 'READ')
        else:
            f_zjet = ROOT.TFile(os.path.join(args.input_dir, 'Zto2Q_2023.root'), 'READ')
        if os.path.exists(pref_w):
            f_wjet = ROOT.TFile(pref_w, 'READ')
        else:
            f_wjet = ROOT.TFile(os.path.join(args.input_dir, 'Wto2Q_2023.root'), 'READ')
        if os.path.exists(pref_gj):
            f_gj = ROOT.TFile(pref_gj, 'READ')
        else:
            f_gj = ROOT.TFile(os.path.join(args.input_dir, 'gj_2023.root'), 'READ')
        if os.path.exists(pref_zg):
            f_zg = ROOT.TFile(pref_zg, 'READ')
        else:
            f_zg = ROOT.TFile(os.path.join(args.input_dir, 'ZGto2QG_2023.root'), 'READ')
        if os.path.exists(pref_data):
            f_data = ROOT.TFile(pref_data, 'READ')
        else:
            f_data = ROOT.TFile(os.path.join(args.input_dir, 'JetMET_2023.root'), 'READ')

if cat_dir is None:
    h_zjet = f_zjet.Get("hist_ak8_"+args.mass_obs)
    h_wjet = f_wjet.Get("hist_ak8_"+args.mass_obs)
    h_gj = f_gj.Get("hist_ak8_"+args.mass_obs)
    h_zg = f_zg.Get("hist_ak8_"+args.mass_obs)
    h_data = f_data.Get("hist_ak8_"+args.mass_obs)

missing = []
for name, obj in [('h_data', locals().get('h_data', None)), ('h_zjet', locals().get('h_zjet', None)), ('h_wjet', locals().get('h_wjet', None)), ('h_gj', locals().get('h_gj', None)), ('h_zg', locals().get('h_zg', None))]:
    if obj is None:
        missing.append(name)
if missing:
    print('ERROR: missing required histograms:', missing)
    print('Input directory was:', args.input_dir)
    print('If you are using GParT structure, run with -i <path>/GParT/CAT0 and ensure histograms named hist_ak8_'+args.mass_obs+' exist in the ROOTs')
    sys.exit(1)
    
obs = ROOT.RooRealVar(args.mass_obs,"",h_data.GetXaxis().GetBinCenter(h_data.GetMaximumBin()),h_data.GetXaxis().GetXmin(),h_data.GetXaxis().GetXmax())
obs.setBins(h_data.GetNbinsX())
obs_list = ROOT.RooArgList()
obs_list.add(obs)

w = ROOT.RooWorkspace("w","")

## data histogram conversion
rh_data = ROOT.RooDataHist("data","",obs_list,h_data)
w.Import(rh_data)

label = ROOT.TLatex()
label.SetTextAlign(12)
label.SetNDC()
label.SetTextSize(label.GetTextSize()*0.8)

c = ROOT.TCanvas("c","",600,600)

######################
## Z+jets modelling ##
######################

rh_zjet = ROOT.RooDataHist("zjet","",obs_list,h_zjet)

## Relativistic Breit-Wigner
zjet_mz = ROOT.RooRealVar("zjet_mz","zjet_mz",91.18)
zjet_gammaz = ROOT.RooRealVar("zjet_gammaz","zjet_gammaz",2.49)
pdf_zjet_bw = ROOT.RooGenericPdf("pdf_zjet_bw","pdf_zjet_bw","@0/(pow(@0*@0 - @1*@1,2) + @2*@2*@0*@0*@0*@0/(@1*@1))",ROOT.RooArgList(obs,zjet_mz,zjet_gammaz))
zjet_mz.setConstant(True)
zjet_gammaz.setConstant(True)

## Gaussian resolution
zjet_gmean = ROOT.RooRealVar("zjet_gmean","zjet_gmean",0.,-20,20)
zjet_gsigmaL = ROOT.RooRealVar("zjet_gsigmaL","zjet_gsigmaL",2.5,1,15)
zjet_gsigmaR = ROOT.RooRealVar("zjet_gsigmaR","zjet_gsigmaR",2.5,1,15)

if args.float_mass_peak:
    CMS_jms_unc = ROOT.RooRealVar("CMS_jms_unc","CMS_jms_unc",0.,-1.,1.)
else:
    CMS_jms_unc = ROOT.RooRealVar("CMS_jms_unc","CMS_jms_unc",0.,-5.,5.)
if args.float_mass_reso:
    CMS_jmr_unc = ROOT.RooRealVar("CMS_jmr_unc","CMS_jmr_unc",0.,-1.,1.)
else:
    CMS_jmr_unc = ROOT.RooRealVar("CMS_jmr_unc","CMS_jmr_unc",0.,-5.,5.)
CMS_jms_unc.setConstant(True)
CMS_jmr_unc.setConstant(True)

zjet_jms_unc = ROOT.RooRealVar("zjet_jms_unc","zjet_jms_unc",args.jms_uncertainty)
zjet_jmr_unc = ROOT.RooRealVar("zjet_jmr_unc","zjet_jmr_unc",args.jmr_uncertainty)
zjet_jms_unc.setConstant(True)
zjet_jmr_unc.setConstant(True)

zjet_gpeak = ROOT.RooFormulaVar("zjet_gpeak","","@0*(1+@1*(@0+@3)*@2)",ROOT.RooArgList(zjet_gmean,zjet_jms_unc,CMS_jms_unc,zjet_mz))
zjet_gresoL = ROOT.RooFormulaVar("zjet_gresoL","","@0*(1+@1*@2)",ROOT.RooArgList(zjet_gsigmaL,zjet_jmr_unc,CMS_jmr_unc))
zjet_gresoR = ROOT.RooFormulaVar("zjet_gresoR","","@0*(1+@1*@2)",ROOT.RooArgList(zjet_gsigmaR,zjet_jmr_unc,CMS_jmr_unc))

pdf_zjet_reso = ROOT.RooBifurGauss("pdf_zjet_reso","pdf_zjet_reso",obs,zjet_gpeak,zjet_gresoL,zjet_gresoR)

pdf_zjet_sig = ROOT.RooFFTConvPdf("pdf_zjet_sig","pdf_zjet_sig",obs,pdf_zjet_bw,pdf_zjet_reso)
    
## combinatorial backgorund
if args.use_poly_bkg:
    zjet_coef_1 = ROOT.RooRealVar("zjet_coef_1","zjet_coef_1",0.001,-10,10)
    zjet_coef_2 = ROOT.RooRealVar("zjet_coef_2","zjet_coef_2",0.001,-10,10)
    pdf_zjet_bkg = ROOT.RooChebychev("pdf_zjet_bkg","pdf_zjet_bkg",obs,ROOT.RooArgList(zjet_coef_1,zjet_coef_2))
else:
    zjet_coef_mean = ROOT.RooRealVar("zjet_coef_mean","zjet_coef_mean",91.,80.,110.)
    zjet_coef_sigmaL = ROOT.RooRealVar("zjet_coef_sigmaL","zjet_coef_sigmaL",20.,10.,50)
    zjet_coef_sigmaR = ROOT.RooRealVar("zjet_coef_sigmaR","zjet_coef_sigmaR",20.,10.,50)
    zjet_coef_gmean = ROOT.RooFormulaVar("zjet_coef_gmean","","@0*(1+@1*@2)",ROOT.RooArgList(zjet_coef_mean,zjet_jms_unc,CMS_jms_unc))    
    zjet_coef_gsigmaL = ROOT.RooFormulaVar("zjet_coef_gsigmaL","","@0*(1+@1*@2)",ROOT.RooArgList(zjet_coef_sigmaL,zjet_jmr_unc,CMS_jmr_unc))
    zjet_coef_gsigmaR = ROOT.RooFormulaVar("zjet_coef_gsigmaR","","@0*(1+@1*@2)",ROOT.RooArgList(zjet_coef_sigmaR,zjet_jmr_unc,CMS_jmr_unc))    
    pdf_zjet_bkg = ROOT.RooBifurGauss("pdf_zjet_bkg","pdf_zjet_bkg",obs,zjet_coef_gmean,zjet_coef_gsigmaL,zjet_coef_gsigmaR)
zjet_frac = ROOT.RooRealVar("zjet_frac","zjet_frac",0.1,0.,1.)

pdf_zjet = ROOT.RooAddPdf("pdf_zjet","pdf_zjet",ROOT.RooArgList(pdf_zjet_sig,pdf_zjet_bkg),ROOT.RooArgList(zjet_frac),True)
pdf_zjet_norm = ROOT.RooRealVar(pdf_zjet.GetName()+"_norm","",rh_zjet.sumEntries())
pdf_zjet_norm.setConstant(True)

fit_zjet_res = pdf_zjet.fitTo(rh_zjet,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"))

h_fit_zjet = pdf_zjet.createHistogram("h_fit_zjet",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_zjet.Scale(pdf_zjet_norm.getVal()*args.rebin_pdf)

h_fit_zjet_bkg = pdf_zjet_bkg.createHistogram("h_fit_zjet_bkg",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_zjet_bkg.Scale(pdf_zjet_norm.getVal()*args.rebin_pdf*(1-zjet_frac.getVal()))

h_fit_zjet_sig = pdf_zjet_sig.createHistogram("h_fit_zjet_sig",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_zjet_sig.Scale(pdf_zjet_norm.getVal()*args.rebin_pdf*zjet_frac.getVal())

h_fit_zjet_test = pdf_zjet.createHistogram("h_fit_zjet_test",obs)
h_fit_zjet_test.Scale(pdf_zjet_norm.getVal())
chi2_zjet,ndf_zjet = getChi2(h_zjet,h_fit_zjet_test)
chi2_zjet = chi2_zjet/(ndf_zjet-fit_zjet_res.floatParsFinal().getSize())
if args.cat_id == "cat0":
    max_zjet, hwhm_zjet = getMaxAndHWHM(h_fit_zjet)
else:
    max_zjet, hwhm_zjet = getMaxAndHWHM(h_fit_zjet_sig)

if args.use_poly_bkg:
    zjet_coef_1.setConstant(True)
    zjet_coef_2.setConstant(True)
else:
    zjet_coef_mean.setConstant(True)
    zjet_coef_sigmaL.setConstant(True)
    zjet_coef_sigmaR.setConstant(True)
    
zjet_frac.setConstant(True)
zjet_gmean.setConstant(True)
zjet_gsigmaL.setConstant(True)
zjet_gsigmaR.setConstant(True)
CMS_jmr_unc.setConstant(False)
CMS_jms_unc.setConstant(False)

w.Import(pdf_zjet)
w.Import(pdf_zjet_norm)

h_zjet.GetXaxis().SetTitle(args.mass_obs+" (GeV)")
h_zjet.GetYaxis().SetTitle("Events")
h_zjet.SetMarkerColor(ROOT.kBlack)
h_zjet.SetLineColor(ROOT.kBlack)
h_zjet.SetMarkerSize(0.6)
h_zjet.SetMarkerStyle(20)
h_zjet.Rebin(args.rebin_data)
h_zjet.GetYaxis().SetRangeUser(0.,h_zjet.GetMaximum()*1.25)
h_zjet.Draw("EP")
h_fit_zjet_bkg.SetLineColor(ROOT.kBlue)
h_fit_zjet_bkg.SetLineWidth(2)
h_fit_zjet_bkg.Scale(args.rebin_data)
if args.plot_bkg:
    h_fit_zjet_bkg.Draw("hist same")
h_fit_zjet.SetLineColor(ROOT.kRed)
h_fit_zjet.SetLineWidth(2)
h_fit_zjet.Scale(args.rebin_data)
h_fit_zjet.Draw("hist same")
h_zjet.Draw("EPsame")
ROOT.CMS_lumi(c,"%.1f"%(lumi_dict[args.year]),False,False,True)
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_zjet))
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_zjet))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_zjet))
c.SaveAs(args.output_directory+"/zjet_fit_"+args.cat_id+".png","png")
c.SaveAs(args.output_directory+"/zjet_fit_"+args.cat_id+".pdf","pdf")

######################
## W+jets modelling ##
######################

rh_wjet = ROOT.RooDataHist("wjet","",obs_list,h_wjet)

wjet_mw = ROOT.RooRealVar("wjet_mw","wjet_mw",80.37)
wjet_gammaw = ROOT.RooRealVar("wjet_gammaw","wjet_gammaw",2.08)
pdf_wjet_bw = ROOT.RooGenericPdf("pdf_wjet_bw","pdf_wjet_bw","@0/(pow(@0*@0 - @1*@1,2) + @2*@2*@0*@0*@0*@0/(@1*@1))",ROOT.RooArgList(obs,wjet_mw,wjet_gammaw))
wjet_mw.setConstant(True)
wjet_gammaw.setConstant(True)


wjet_gmean = ROOT.RooRealVar("wjet_gmean","wjet_gmean",zjet_gmean.getVal(),zjet_gmean.getVal()-5,zjet_gmean.getVal()+5)
wjet_gsigmaL = ROOT.RooRealVar("wjet_gsigmaL","wjet_gsigmaL",zjet_gsigmaL.getVal(),zjet_gsigmaL.getVal()*0.5,zjet_gsigmaL.getVal()*2.)
wjet_gsigmaR = ROOT.RooRealVar("wjet_gsigmaR","wjet_gsigmaR",zjet_gsigmaR.getVal(),zjet_gsigmaR.getVal()*0.5,zjet_gsigmaR.getVal()*2.)

wjet_jms_unc = ROOT.RooRealVar("wjet_jms_unc","wjet_jms_unc",args.jms_uncertainty)
wjet_jmr_unc = ROOT.RooRealVar("wjet_jmr_unc","wjet_jmr_unc",args.jmr_uncertainty)
wjet_jms_unc.setConstant(True)
wjet_jmr_unc.setConstant(True)
CMS_jms_unc.setConstant(True)
CMS_jmr_unc.setConstant(True)

wjet_gpeak = ROOT.RooFormulaVar("wjet_gpeak","","@0*(1+@1*(@0+@3)*@2)",ROOT.RooArgList(wjet_gmean,wjet_jms_unc,CMS_jms_unc,wjet_mw))
wjet_gresoL = ROOT.RooFormulaVar("wjet_gresoL","","@0*(1+@1*@2)",ROOT.RooArgList(wjet_gsigmaL,wjet_jmr_unc,CMS_jmr_unc))
wjet_gresoR = ROOT.RooFormulaVar("wjet_gresoR","","@0*(1+@1*@2)",ROOT.RooArgList(wjet_gsigmaR,wjet_jmr_unc,CMS_jmr_unc))

pdf_wjet_reso = ROOT.RooBifurGauss("pdf_wjet_reso","pdf_wjet_reso",obs,wjet_gpeak,wjet_gresoL,wjet_gresoR)

pdf_wjet = ROOT.RooFFTConvPdf("pdf_wjet","pdf_wjet",obs,pdf_wjet_bw,pdf_wjet_reso)

pdf_wjet_norm = ROOT.RooRealVar(pdf_wjet.GetName()+"_norm","",rh_wjet.sumEntries())
pdf_wjet_norm.setConstant(True)

fit_wjet_res = pdf_wjet.fitTo(rh_wjet,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"))

h_fit_wjet = pdf_wjet.createHistogram("h_fit_wjet",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_wjet.Scale(pdf_wjet_norm.getVal()*args.rebin_pdf)

h_fit_wjet_test = pdf_wjet.createHistogram("h_fit_wjet_test",obs)
h_fit_wjet_test.Scale(pdf_wjet_norm.getVal())
chi2_wjet,ndf_wjet = getChi2(h_wjet,h_fit_wjet_test)
chi2_wjet = chi2_wjet/(ndf_wjet-fit_wjet_res.floatParsFinal().getSize())
max_wjet, hwhm_wjet = getMaxAndHWHM(h_fit_wjet)

wjet_gmean.setConstant(True)
wjet_gsigmaL.setConstant(True)
wjet_gsigmaR.setConstant(True)
CMS_jms_unc.setConstant(False)
CMS_jmr_unc.setConstant(False)
    
w.Import(pdf_wjet)
w.Import(pdf_wjet_norm)

h_wjet.GetXaxis().SetTitle(args.mass_obs+" (GeV)")
h_wjet.GetYaxis().SetTitle("Events")
h_wjet.SetMarkerColor(ROOT.kBlack)
h_wjet.SetLineColor(ROOT.kBlack)
h_wjet.SetMarkerSize(0.6)
h_wjet.SetMarkerStyle(20)
h_wjet.Rebin(args.rebin_data)
h_wjet.GetYaxis().SetRangeUser(0.,h_wjet.GetMaximum()*1.25)
h_wjet.Draw("EP")
h_fit_wjet.SetLineColor(ROOT.kRed)
h_fit_wjet.SetLineWidth(2)
h_fit_wjet.Scale(args.rebin_data)
h_fit_wjet.Draw("hist same")
h_wjet.Draw("EPsame")
ROOT.CMS_lumi(c,"%.1f"%(lumi_dict[args.year]),False,False,True)
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_wjet))
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_wjet))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_wjet))
c.SaveAs(args.output_directory+"/wjet_fit_"+args.cat_id+".png","png")
c.SaveAs(args.output_directory+"/wjet_fit_"+args.cat_id+".pdf","pdf")

######################
## gj MC modelling ##
######################

rh_gj = ROOT.RooDataHist("gj","",obs_list,h_gj)

## order tune in MC
gj_coef_1 = ROOT.RooRealVar("gj_coef_1","gj_coef_1",0.001,-10,10)
gj_coef_2 = ROOT.RooRealVar("gj_coef_2","gj_coef_2",0.001,-10,10)
gj_coef_3 = ROOT.RooRealVar("gj_coef_3","gj_coef_3",0.001,-10,10)
gj_coef_4 = ROOT.RooRealVar("gj_coef_4","gj_coef_4",0.001,-10,10)

if args.tagger == "pnet" and args.cat_id == "cat3":
    pdf_gj = ROOT.RooChebychev("pdf_gj","",obs,ROOT.RooArgList(gj_coef_1,gj_coef_2,gj_coef_3,gj_coef_4))
else:
    pdf_gj = ROOT.RooChebychev("pdf_gj","",obs,ROOT.RooArgList(gj_coef_1,gj_coef_2,gj_coef_3))

fit_gj_res = pdf_gj.fitTo(rh_gj,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"))

pdf_gj_norm = ROOT.RooRealVar(pdf_gj.GetName()+"_norm","",h_gj.Integral(),0.5*h_gj.Integral(),2.0*h_gj.Integral())
pdf_gj_norm.setConstant(False)

h_fit_gj = pdf_gj.createHistogram("h_fit_gj",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_gj.Scale(pdf_gj_norm.getVal()*args.rebin_pdf)

w.Import(pdf_gj)

h_fit_gj_test = pdf_gj.createHistogram("h_fit_gj_test",obs)
h_fit_gj_test.Scale(pdf_gj_norm.getVal())
chi2_gj,ndf_gj = getChi2(h_gj,h_fit_gj_test)
chi2_gj = chi2_gj/(ndf_gj-fit_gj_res.floatParsFinal().getSize())

h_gj.GetXaxis().SetTitle(args.mass_obs+" (GeV)")
h_gj.GetYaxis().SetTitle("Events")
h_gj.SetMarkerColor(ROOT.kBlack)
h_gj.SetLineColor(ROOT.kBlack)
h_gj.SetMarkerSize(0.6)
h_gj.SetMarkerStyle(20)
h_gj.Rebin(args.rebin_data)
h_gj.GetYaxis().SetRangeUser(0.,h_gj.GetMaximum()*1.25)
h_gj.Draw("EP")
h_fit_gj.SetLineColor(ROOT.kRed)
h_fit_gj.SetLineWidth(2)
h_fit_gj.Scale(args.rebin_data)
h_fit_gj.Draw("hist same")
h_gj.Draw("EPsame")

ROOT.CMS_lumi(c,"%.1f"%(lumi_dict[args.year]),False,False,True)
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_gj))
c.SaveAs(args.output_directory+"/gj_fit_"+args.cat_id+".png","png")
c.SaveAs(args.output_directory+"/gj_fit_"+args.cat_id+".pdf","pdf")

####################
## Zγ modelling ##
####################

rh_zg = ROOT.RooDataHist("zg","",obs_list,h_zg)

## Relativistic Breit-Wigner for Z boson (same parameters as Z+jets)
zg_mz = ROOT.RooRealVar("zg_mz","zg_mz",91.18)
zg_gammaz = ROOT.RooRealVar("zg_gammaz","zg_gammaz",2.49)
pdf_zg_bw = ROOT.RooGenericPdf("pdf_zg_bw","pdf_zg_bw","@0/(pow(@0*@0 - @1*@1,2) + @2*@2*@0*@0*@0*@0/(@1*@1))",ROOT.RooArgList(obs,zg_mz,zg_gammaz))
zg_mz.setConstant(True)
zg_gammaz.setConstant(True)

## Gaussian resolution
zg_gmean = ROOT.RooRealVar("zg_gmean","zg_gmean",0.,-20.,20.)
zg_gsigmaL = ROOT.RooRealVar("zg_gsigmaL","zg_gsigmaL",5.,1.,15.)
zg_gsigmaR = ROOT.RooRealVar("zg_gsigmaR","zg_gsigmaR",5.,1.,15.)

zg_jms_unc = ROOT.RooRealVar("zg_jms_unc","zg_jms_unc",args.jms_uncertainty)
zg_jmr_unc = ROOT.RooRealVar("zg_jmr_unc","zg_jmr_unc",args.jmr_uncertainty)
zg_jms_unc.setConstant(True)
zg_jmr_unc.setConstant(True)
CMS_jms_unc.setConstant(True)
CMS_jmr_unc.setConstant(True)

zg_gpeak = ROOT.RooFormulaVar("zg_gpeak","","@0*(1+@1*(@0+@3)*@2)",ROOT.RooArgList(zg_gmean,zg_jms_unc,CMS_jms_unc,zg_mz))
zg_gresoL = ROOT.RooFormulaVar("zg_gresoL","","@0*(1+@1*@2)",ROOT.RooArgList(zg_gsigmaL,zg_jmr_unc,CMS_jmr_unc))
zg_gresoR = ROOT.RooFormulaVar("zg_gresoR","","@0*(1+@1*@2)",ROOT.RooArgList(zg_gsigmaR,zg_jmr_unc,CMS_jmr_unc))

pdf_zg_reso = ROOT.RooBifurGauss("pdf_zg_reso","pdf_zg_reso",obs,zg_gpeak,zg_gresoL,zg_gresoR)

pdf_zg_sig = ROOT.RooFFTConvPdf("pdf_zg_sig","pdf_zg_sig",obs,pdf_zg_bw,pdf_zg_reso)
    
## combinatorial background
if args.use_poly_bkg:
    zg_coef_1 = ROOT.RooRealVar("zg_coef_1","zg_coef_1",0.001,-10,10)
    zg_coef_2 = ROOT.RooRealVar("zg_coef_2","zg_coef_2",0.001,-10,10)
    pdf_zg_bkg = ROOT.RooChebychev("pdf_zg_bkg","pdf_zg_bkg",obs,ROOT.RooArgList(zg_coef_1,zg_coef_2))
else:
    zg_coef_mean = ROOT.RooRealVar("zg_coef_mean","zg_coef_mean",91.,80.,110.)
    zg_coef_sigmaL = ROOT.RooRealVar("zg_coef_sigmaL","zg_coef_sigmaL",20.,10.,50)
    zg_coef_sigmaR = ROOT.RooRealVar("zg_coef_sigmaR","zg_coef_sigmaR",20.,10.,50)
    zg_coef_gmean = ROOT.RooFormulaVar("zg_coef_gmean","","@0*(1+@1*@2)",ROOT.RooArgList(zg_coef_mean,zg_jms_unc,CMS_jms_unc))    
    zg_coef_gsigmaL = ROOT.RooFormulaVar("zg_coef_gsigmaL","","@0*(1+@1*@2)",ROOT.RooArgList(zg_coef_sigmaL,zg_jmr_unc,CMS_jmr_unc))
    zg_coef_gsigmaR = ROOT.RooFormulaVar("zg_coef_gsigmaR","","@0*(1+@1*@2)",ROOT.RooArgList(zg_coef_sigmaR,zg_jmr_unc,CMS_jmr_unc))    
    pdf_zg_bkg = ROOT.RooBifurGauss("pdf_zg_bkg","pdf_zg_bkg",obs,zg_coef_gmean,zg_coef_gsigmaL,zg_coef_gsigmaR)
zg_frac = ROOT.RooRealVar("zg_frac","zg_frac",0.1,0.,1.)

pdf_zg = ROOT.RooAddPdf("pdf_zg","pdf_zg",ROOT.RooArgList(pdf_zg_sig,pdf_zg_bkg),ROOT.RooArgList(zg_frac),True)
pdf_zg_norm = ROOT.RooRealVar(pdf_zg.GetName()+"_norm","",rh_zg.sumEntries())
pdf_zg_norm.setConstant(True)

fit_zg_res = pdf_zg.fitTo(rh_zg,ROOT.RooFit.Save(),ROOT.RooFit.Optimize(1),ROOT.RooFit.SumW2Error(True),ROOT.RooFit.Minimizer("Minuit2"))

h_fit_zg = pdf_zg.createHistogram("h_fit_zg",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_zg.Scale(pdf_zg_norm.getVal()*args.rebin_pdf)

h_fit_zg_bkg = pdf_zg_bkg.createHistogram("h_fit_zg_bkg",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_zg_bkg.Scale(pdf_zg_norm.getVal()*args.rebin_pdf*(1-zg_frac.getVal()))

h_fit_zg_sig = pdf_zg_sig.createHistogram("h_fit_zg_sig",obs,ROOT.RooFit.Binning(obs.getBins()*args.rebin_pdf))
h_fit_zg_sig.Scale(pdf_zg_norm.getVal()*args.rebin_pdf*zg_frac.getVal())

h_fit_zg_test = pdf_zg.createHistogram("h_fit_zg_test",obs)
h_fit_zg_test.Scale(pdf_zg_norm.getVal())
chi2_zg,ndf_zg = getChi2(h_zg,h_fit_zg_test)
chi2_zg = chi2_zg/(ndf_zg-fit_zg_res.floatParsFinal().getSize())
if args.cat_id == "cat0":
    max_zg, hwhm_zg = getMaxAndHWHM(h_fit_zg)
else:
    max_zg, hwhm_zg = getMaxAndHWHM(h_fit_zg_sig)

if args.use_poly_bkg:
    zg_coef_1.setConstant(True)
    zg_coef_2.setConstant(True)
else:
    zg_coef_mean.setConstant(True)
    zg_coef_sigmaL.setConstant(True)
    zg_coef_sigmaR.setConstant(True)
    
zg_frac.setConstant(True)
zg_gmean.setConstant(True)
zg_gsigmaL.setConstant(True)
zg_gsigmaR.setConstant(True)
CMS_jmr_unc.setConstant(False)
CMS_jms_unc.setConstant(False)

w.Import(pdf_zg)
w.Import(pdf_zg_norm)

h_zg.GetXaxis().SetTitle(args.mass_obs+" (GeV)")
h_zg.GetYaxis().SetTitle("Events")
h_zg.SetMarkerColor(ROOT.kBlack)
h_zg.SetLineColor(ROOT.kBlack)
h_zg.SetMarkerSize(0.6)
h_zg.SetMarkerStyle(20)
h_zg.Rebin(args.rebin_data)
h_zg.GetYaxis().SetRangeUser(0.,h_zg.GetMaximum()*1.25)
h_zg.Draw("EP")
h_fit_zg_bkg.SetLineColor(ROOT.kBlue)
h_fit_zg_bkg.SetLineWidth(2)
h_fit_zg_bkg.Scale(args.rebin_data)
if args.plot_bkg:
    h_fit_zg_bkg.Draw("hist same")
h_fit_zg.SetLineColor(ROOT.kRed)
h_fit_zg.SetLineWidth(2)
h_fit_zg.Scale(args.rebin_data)
h_fit_zg.Draw("hist same")
h_zg.Draw("EPsame")
ROOT.CMS_lumi(c,"%.1f"%(lumi_dict[args.year]),False,False,True)
label.DrawLatex(0.65,0.8,"#chi^{2}/ndf=%.2f"%(chi2_zg))
label.DrawLatex(0.65,0.75,"Peak=%.2f"%(max_zg))
label.DrawLatex(0.65,0.7,"HWHM=%.2f"%(hwhm_zg))
c.SaveAs(args.output_directory+"/zg_fit_"+args.cat_id+".png","png")
c.SaveAs(args.output_directory+"/zg_fit_"+args.cat_id+".pdf","pdf")

#### prefit distributions
pdf_gj_norm.setVal(rh_data.sumEntries()-pdf_wjet_norm.getVal()-pdf_zjet_norm.getVal()-pdf_zg_norm.getVal())
w.Import(pdf_gj_norm)
h_fit_gj.Scale(pdf_gj_norm.getVal()*args.rebin_pdf/h_fit_gj.Integral())

### Save output files
f_out = ROOT.TFile(args.output_directory+"/workspace_"+args.cat_id+".root","RECREATE")
f_out.cd()
w.Write("w")
f_out.Close()

### write datacard
os.chdir(args.output_directory)
def try_get_ewk(fhandle, summed_files, hist_base, suffix=''):
    # if fhandle available try direct Get
    if fhandle:
        key = 'hist_ak8_{}{}'.format(hist_base, suffix)
        if fhandle.GetListOfKeys().Contains(key):
            return fhandle.Get(key)
    # else try to find in summed files (e.g. ZGto2QG) using helper
    if 'files' in locals() and summed_files:
        return sum_histograms_from_files(summed_files, 'hist_ak8_{}{}'.format(hist_base, suffix))
    # fallback to the nominal hist if nothing else
    if 'z' in hist_base.lower():
        return locals().get('h_zjet')
    elif 'w' in hist_base.lower():
        return locals().get('h_wjet')
    elif 'zg' in hist_base.lower():
        return locals().get('h_zg')
    else:
        return None

h_zjet_ewk = try_get_ewk(f_zjet, z_files if 'z_files' in locals() else None, args.mass_obs, '_ewk')
h_wjet_ewk = try_get_ewk(f_wjet, w_files if 'w_files' in locals() else None, args.mass_obs, '_ewk')
h_zg_ewk = try_get_ewk(f_zg, zg_files if 'zg_files' in locals() else None, args.mass_obs, '_ewk')

def safe_corr(h_ewk, h_nom):
    try:
        if h_nom is None:
            return 1.0
        nom = h_nom.Integral()
        if nom == 0:
            return 1.0
        if h_ewk is None:
            return 1.0
        return 0.5*((h_ewk.Integral() + nom) / nom)
    except Exception:
        return 1.0

# fallback: if no ewk histograms found, use nominal histograms to avoid None
if h_zjet_ewk is None:
    h_zjet_ewk = h_zjet
if h_wjet_ewk is None:
    h_wjet_ewk = h_wjet
if h_zg_ewk is None:
    h_zg_ewk = h_zg


sys = {"NLOEWK_corr_zjet": safe_corr(h_zjet_ewk, h_zjet), "NLOEWK_corr_wjet": safe_corr(h_wjet_ewk, h_wjet), "NLOEWK_corr_zg": safe_corr(h_zg_ewk, h_zg)}
dc = createCardTemplate(args.cat_id,"workspace_"+args.cat_id+".root",sys,"w")
