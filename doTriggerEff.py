import ROOT
import os, sys

# Define lepton type
lepton = "Muon"  # Change to "Electron" for electrons

# Lepton-specific configurations
if lepton == "Muon":
    inputFiles = "WHleptonicpythia_leptonic_M125.0_MD3.00_T3.00_HT-1_UL17_NANOAOD.root"
    outputHistos = "muon_efficiencies.root"
    hlt = ["HLT_IsoMu27", "HLT_Mu50"]
    offlineCuts = {
        "lep1pt": 40,
        "MET": 40,
    }
    histBins = {
        "lep1pt": [20, 0, 200],
        "MET": [20, 0, 200]
    }
    pt_label = "Muon pT [GeV]"
else:  # Electron-specific configurations
    inputFiles = "WHleptonicpythia_leptonic_M125.0_MD3.00_T3.00_HT-1_UL17_NANOAOD.root"
    outputHistos = "electron_efficiencies.root"
    if "UL17" in inputFiles or "UL18" in inputFiles:
        hlt = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200"]
    else:
        hlt = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175"]
    offlineCuts = {
        "lep1pt": 35,
        "MET": 40,
    }
    histBins = {
        "lep1pt": [20, 0, 200],
        "MET": [20, 0, 200]
    }
    pt_label = "Electron pT [GeV]"

# Create the actual histograms for saving
histos = {}

for var in histBins:
    x_label = pt_label if var == "lep1pt" else var
    histos[var + "_num"] = ROOT.TH1F(var + "_num", var+"_num", histBins[var][0], histBins[var][1], histBins[var][2])
    histos[var + "_den"] = ROOT.TH1F(var + "_den", var+"_den", histBins[var][0], histBins[var][1], histBins[var][2])

def passes_lepton_cuts(ev, lepton, leptonIndex):
    if lepton == "Muon":
        return (
            getattr(ev, lepton + "_tightId")[leptonIndex]
            and abs(getattr(ev, lepton + "_eta")[leptonIndex]) < 2.4
            and abs(getattr(ev, lepton + "_dz")[leptonIndex]) <= 0.05
            and abs(getattr(ev, lepton + "_dxy")[leptonIndex]) <= 0.02
            and getattr(ev, lepton + "_pfIsoId")[leptonIndex] >= 5
        )
    else:
        return (
            ev.Electron_cutBased[leptonIndex] >= 2
            and ev.Electron_mvaFall17V2Iso_WP80[leptonIndex]
            and abs(ev.Electron_dxy[leptonIndex]) < 0.05 + 0.05 * (abs(ev.Electron_eta[leptonIndex]) > 1.479)
            and abs(ev.Electron_dz[leptonIndex]) < 0.10 + 0.10 * (abs(ev.Electron_eta[leptonIndex]) > 1.479)
            and ((abs(ev.Electron_eta[leptonIndex]) < 1.444) or (abs(ev.Electron_eta[leptonIndex]) > 1.566))
            and abs(ev.Electron_eta[leptonIndex]) < 2.5
        )

# Now loop over the events
print("Starting %s" % inputFiles)

tf = ROOT.TFile(inputFiles, "READ")
events = tf.Get("Events")

# Event counter
iEv = 0
# Total number of events
nEv = events.GetEntries()

for ev in events:
    if iEv % 1000 == 0:
        print("%i/%i events in file done" % (iEv, nEv))
    iEv += 1
#    if(iEv % 500 == 0): break

    passHLT = False
    # Check if we pass numerator
    for hltpath in hlt:
        if getattr(ev, hltpath, False): passHLT = True

    # Find the lepton with the highest pT
    highest_pt = -1
    highest_pt_lepton_index = -1

    if not getattr(ev, "n" + lepton) > 0:
        continue

    for leptonIndex in range(getattr(ev, "n" + lepton)):
        if passes_lepton_cuts(ev, lepton, leptonIndex):
            if getattr(ev, lepton + "_pt")[leptonIndex] > highest_pt:
                highest_pt = getattr(ev, lepton + "_pt")[leptonIndex]
                highest_pt_lepton_index = leptonIndex

    # If no valid lepton found, continue to the next event
    if highest_pt_lepton_index == -1:
        continue
    leptonIndex = highest_pt_lepton_index

    # Variables you want to study
    passmetCut = ev.MET_pt >= offlineCuts["MET"]
    passlepCut = getattr(ev, lepton + "_pt")[leptonIndex] >= offlineCuts["lep1pt"]

    # Then save denominator and numerator
    for var in histBins:
        passDen = False
        fillvar = None
        if var == "lep1pt":
            passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passmetCut
            fillvar = getattr(ev, lepton + "_pt")[leptonIndex]
        elif var == "MET":
            passDen = passes_lepton_cuts(ev, lepton, leptonIndex) and passlepCut
            fillvar = ev.MET_pt
        if passDen and fillvar is not None:
            histos[var + "_den"].Fill(fillvar)
            if passHLT:
                histos[var + "_num"].Fill(fillvar)

tf.Close()

outF = ROOT.TFile(outputHistos, "RECREATE")
for h in histos:
    histos[h].Write()
outF.Close()

print("Processing complete. Output written to %s" % outputHistos)
