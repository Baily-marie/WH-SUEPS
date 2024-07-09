import uproot
import argparse
import numpy as np
import ROOT
import awkward as ak
from array import array

# Sets batch mode so no popup window
ROOT.gROOT.SetBatch(True)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("--input", help="Name of input file", type=str)
args = vars(parser.parse_args())

# Name of sample
sample_name = args["input"]

output_file = "MC_muon_efficiencies.root"
input_file = "/eos/user/j/jreicher/SUEP/WH_private_signals/merged/" + sample_name + ".root"

# suep decay type
decay_type = ""
if "generic" in sample_name:
    decay_type = "generic"
elif "hadronic" in sample_name:
    decay_type = "hadronic"
else:
    decay_type = "leptonic"

# conditions for what year
if "UL18" in sample_name:
    year="2018 conditions"
    folder = "MC_eff_outputs_2018/"
elif "UL17" in sample_name:
    year = "2017 conditions"
    folder = "MC_eff_outputs_2017/"
elif "UL16APV" in sample_name:
    year = "2016 APV conditions"
    folder = "MC_eff_outputs_2016APV/"
else:
    year = "2016 conditions"
    folder = "MC_eff_outputs_2016/"

# dark meson (phi) mass
md = ""
if "MD2.00" in sample_name:
    md = "2.00 [GeV]"
elif "MD4.00" in sample_name:
    md = "4.00 [GeV]"
elif "MD3.00" in sample_name:
    md = "3.00 [GeV]"
elif "MD8.00" in sample_name:
    md = "8.00 [GeV]"
elif "MD1.00" in sample_name:
    md = "1.00 [GeV]"
else:
    md = "1.40 [GeV]"

# temperature
temp = ""
if "T0.25" in sample_name:
    temp = "0.25"
if "T0.35" in sample_name:
    temp = "0.35"
if "T0.50" in sample_name:
    temp = "0.50"
elif "T0.75" in sample_name:
    temp = "0.75"
elif "T1.00" in sample_name:
    temp = "1.00"
elif "T1.50" in sample_name:
    temp = "1.50"
elif "T2.00" in sample_name:
    temp = "2.00"
elif "T3.00" in sample_name:
    temp = "3.00"
elif "T4.00" in sample_name:
    temp = "4.00"
elif "T8.00" in sample_name:
    temp = "8.00"
elif "T12.00" in sample_name:
    temp = "12.00"
elif "T16.00" in sample_name:
    temp = "16.00"
elif "T32.00" in sample_name:
    temp = "32.00"
else:
    temp = "6.00"

# Gets relevant events from file
def Events(f):
    evs=f['Events'].arrays(['HLT_IsoMu27',
                'HLT_Mu50',
                'Muon_pt',
                'Muon_eta',
                'Muon_phi',
                'Muon_dz',
                'Muon_dxy',
                'Muon_pfIsoId',
                'Muon_tightId',
                'Muon_pdgId',
                'TrigObj_pt',
                'TrigObj_eta',
                'TrigObj_phi',
                'TrigObj_id',
                'TrigObj_filterBits'])

    return evs

with uproot.open(input_file) as f:
    evs = Events(f)

# Defining a offline Muon
muons = ak.zip({
    "pt": evs["Muon_pt"],
    "eta": evs["Muon_eta"],
    "phi": evs["Muon_phi"],
    "mass": 0.0,
    "charge": evs["Muon_pdgId"]/(-13),
    "pdgId": evs["Muon_pdgId"],
    "isTight": evs["Muon_tightId"],
    "isTightIso": evs["Muon_pfIsoId"] >= 4
    }, with_name="Momentum4D")

cutMuons = ((evs["Muon_tightId"])
                & (abs(evs["Muon_eta"]) < 2.4)
                & (abs(evs["Muon_dz"]) <= 0.05)
                & (abs(evs["Muon_dxy"]) <= 0.02)
                & (evs["Muon_pfIsoId"] >= 5))
offlineMuons = muons[cutMuons]

# Gets matched online/offline muons from file
def isHLTMatched(events, offlineMuons):
    trigObj = ak.zip({
            "pt": events["TrigObj_pt"],
            "eta": events["TrigObj_eta"],
            "phi": events["TrigObj_phi"],
            "mass": 0.0,
            "id": events["TrigObj_id"],
            "filterBits": events["TrigObj_filterBits"]
            }, with_name = "Momentum4D")

    # Defining the conditions for filtering triggers
    filterbits1 = (((events['TrigObj_filterBits'] & 2) == 2) & ((events['TrigObj_filterBits'] & 8) == 8) & (events['HLT_IsoMu27']))
    filterbits2 = (((events['TrigObj_filterBits'] & 1024) == 1024) & (events['HLT_Mu50']))

    trigObjSingleMu = trigObj[((abs(trigObj.id) == 13)
                              & (trigObj.pt >= 10)
                              & (abs(trigObj.eta) < 2.4)
                              & (filterbits1 | filterbits2))]

    # Computes deltaR2
    def deltaR2(eta1, phi1, eta2, phi2):
            deta = eta1 - eta2
            dphi = phi1 - phi2
            dphi = np.arccos(np.cos(dphi))

            return deta**2 + dphi**2

    toMatch1Mu, trigObjSingleMu = ak.unzip(ak.cartesian([offlineMuons, trigObjSingleMu], axis = 1, nested = True))
    alldr2 = deltaR2(toMatch1Mu.eta, toMatch1Mu.phi, trigObjSingleMu.eta, trigObjSingleMu.phi)
    min_dr2 = ak.min(alldr2, axis = 2)
    match1Mu = (min_dr2 < 0.1)

    return match1Mu, min_dr2

# Define eta bins
eta_bins = [[0.0, 0.9], [0.9, 2.1], [2.1, 2.4]]


# Create TH2F


th2_hist = ROOT.TH2F("Muon_pt_vs_minDeltaR", "Muon pT vs Min DeltaR;Muon pT [GeV];Min DeltaR",
                     100, 0, 200, 100, 0, 0.5)

# Fill the histogram


match1Mu, min_dr2 = isHLTMatched(evs, offlineMuons)

matched_pts = ak.flatten(offlineMuons[match1Mu].pt).to_numpy()
min_deltaRs = (ak.flatten(min_dr2[match1Mu]).to_numpy())

for pt, deltaR in zip(matched_pts, min_deltaRs):
    th2_hist.Fill(pt, deltaR)

# Draw and save as PNG


c1 = ROOT.TCanvas("c1", "Muon pT vs Min DeltaR", 800, 600)
th2_hist.Draw("COLZ")
th2_hist.SetStats(0)
c1.SaveAs(sample_name + "_pt_vs_minDeltaR.pdf")

# Save to ROOT file


output_root_file = ROOT.TFile(output_file, "RECREATE")
th2_hist.Write()
output_root_file.Close()

print("Sample " + sample_name + " processed and histogram saved to " + output_file)
print("Histogram also saved as " + sample_name + "_pt_vs_minDeltaR.pdf")
