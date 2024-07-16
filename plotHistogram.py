import ROOT
import os
import sys
import CMS_lumi
import tdrstyle
#####################
### Configuration ###
#####################

histBins = {
"lep1pt": [20, 0, 250],
"MET": [20, 0, 300]
}


# TDR
tdrstyle.setTDRStyle()

# Change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.cmsText = "CMS"
CMS_lumi.writeExtraText = True
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_7TeV = "5.1 fb^{-1}"
CMS_lumi.lumi_8TeV = "19.7 fb^{-1}"
CMS_lumi.lumi_13TeV = "20.1 fb^{-1}"
CMS_lumi.lumi_sqrtS = "13 TeV"  # Used with iPeriod = 0, e.g., for simulation-only plots

iPos = 0
if iPos == 0:
    CMS_lumi.relPosX = 0.13

H_ref = 600
W_ref = 800
W = W_ref
H = H_ref

# Define the period for lumi
iPeriod = 4  # 13 TeV only

# Iterate through histogram bins and create efficiency plots
inF = ROOT.TFile(sys.argv[1], "READ")

for var in histBins:
    teff = ROOT.TEfficiency(inF.Get(var + "_num"), inF.Get(var + "_den"))
    c = ROOT.TCanvas(var + "C", var + "C", 800, 600)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.18)
    c.SetTickx(0)
    c.SetTicky(0)

    teff.SetTitle(";%s;%s" % (var, "Efficiency"))
    teff.Draw("AP")
    c.Update()
    gr = teff.GetPaintedGraph()
    gr.SetMinimum(0)
    gr.GetXaxis().SetLimits(histBins[var][1], histBins[var][2])

    # Apply the CMS lumi style
    CMS_lumi.CMS_lumi(c, iPeriod, iPos)

    c.SaveAs(sys.argv[2] + "/%s.pdf" % (var))



inF.Close()
