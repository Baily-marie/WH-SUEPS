import ROOT
import os, sys

#Configurations
data = False
inputFiles = "WHleptonicpythia_leptonic_M125.0_MD3.00_T3.00_HT-1_UL17_NANOAOD.root"
outputHistos = "MC_mu_efficiencies.root"

#HLT Path(s) to measure #
hlt = ["HLT_IsoMu27", "HLT_Mu50"]

# Definition of offline cuts
offlineCuts = {
    "lep1pt": 40,
    "met": 40,
}

# Histogram configuration
histBins = {
    "lep1pt" : [20, 0, 200],
    "met"    : [20, 0, 200]
}
# Create the actual histograms for saving #
histos = {}

for var in histBins:
    histos[var + "_num"] = ROOT.TH1F(var + "_num", var+"_num", histBins[var][0], histBins[var][1], histBins[var][2])
    histos[var + "_den"] = ROOT.TH1F(var + "_den", var+"_den", histBins[var][0], histBins[var][1], histBins[var][2])

# Now loop over the events #
print("Starting %s"%(inputFiles))

tf = ROOT.TFile(inputFiles, "READ")
events = tf.Get("Events")

#event counter
iEv = 0
#total number of events
nEv = events.GetEntries()


for ev in events:
    if(iEv%1000 == 0): print("%i/%i events in file done"%(iEv, nEv))
    iEv += 1
    if(iEv%10000 == 0): break

    #first make sure everything passes the reference cuts. This is important for data.
    if data:
        if not(passRefCuts(ev)): continue

    passHLT = False
    # Check if we pass numerator
    for hltpath in hlt:
        if getattr(ev, hltpath, False): passHLT = True

    # Find the lepton with the highest pT
    highest_pt = -1
    highest_pt_lepton_index = -1

    #all events pass the good lepton definition
    if not(ev.nMuon > 0): continue

    passgoodLepCut  = (ev.Muon_tightId
                       and (abs(ev.Muon_eta) < 2.4)
                       and (abs(ev.Muon_dz) <= 0.05)
                       and (abs(ev.Muon_dxy) <= 0.02)
                       and (ev.Muon_pfIsoId >= 5))

    for muonIndex in range(ev.passgoodLepCut):
        if ev.Muon_pt[muonIndex] > highest_pt:
            highest_pt = ev.Muon_pt[muonIndex]
            highest_pt_lepton_index = muonIndex

    # If no valid muon found, continue to the next event
    if highest_pt_lepton_index == -1:
        continue
    muonIndex = highest_pt_lepton_index

    # Variables you want to study:
    passmetCut = ev.MET_pt >= offlineCuts["met"]
    passlepCut = ev.Muon_pt[muonIndex] >= offlineCuts["lep1pt"]

    # Then save denominator and numerator
    for var in histBins:
        passDen = False
        if var == "lep1pt":
            passDen = passgoodLepCut and passmetCut
            fillvar = ev.Muon_pt[muonIndex]
        if var == "met":
            passDen = passgoodLepCut and passlepCut
            fillvar = ev.MET_pt
        if passDen:
            histos[var + "_den"].Fill(fillvar)
            if passHLT:
                histos[var + "_num"].Fill(fillvar)
tf.Close()


outF = ROOT.TFile(outputHistos, "RECREATE")
for h in histos:
  histos[h].Write()
