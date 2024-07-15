import ROOT
import os
import sys
import cmsstyle as CMS
import tdrstyle
#####################
### Configuration ###
#####################

# TDR
tdrstyle.setTDRStyle()

# Set ROOT to batch mode to suppress plots
ROOT.gROOT.SetBatch(True)

histBins = {
  "lep1pt": [25, 0, 250],
  "met"  : [35, 0, 350],
}

### Plotting ###

inF = ROOT.TFile(sys.argv[1], "READ")

for var in histBins:
  teff = ROOT.TEfficiency(inF.Get(var +"_num"), inF.Get(var+"_den"))
  CMS.SetExtraText("Simulation Preliminary")
  CMS.SetLumi("")
#  c = CMS.cmsCanvas('', 0, 1, 0, 1, '', '', square = CMS.kSquare, extraSpace=0.01, iPos=11)
  c = ROOT.TCanvas(var + "C", var + "C", 800, 600)
  teff.SetTitle(";%s;%s"%(var,"Efficiency"))
  teff.Draw("AP")
  c.Update()
  gr = teff.GetPaintedGraph()
  gr.SetMinimum(0)
  gr.GetXaxis().SetLimits(histBins[var][1], histBins[var][2])
  c.SaveAs(sys.argv[2] + "/%s.pdf"%(var))

inF.Close()
