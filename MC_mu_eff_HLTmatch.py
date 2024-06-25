#CURRENTLY ADDING HLT ONLINE/OFFLINE MATCHING, SEE ORIGINAL IN GIT OR TERMINAL
# to run: $ python2 MC_muon_efficiency.py --input "name_of_folder"
import os
import argparse
import uproot
import numpy as np
import ROOT
import awkward as ak
from array import array
# Sets batch mode so no popup window
ROOT.gROOT.SetBatch(True)

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("--input_folder", help="Name of input folder", type=str)
args = vars(parser.parse_args())

# Name of folder containing samples
input_folder = args["input_folder"]
output_file= "MC_mu_efficiencies.root"


# Defines binning and histograms
mu_bin_edges=array('d',[0,2,4,6,8,10,12,
                         14,16,18,20,22,
                         24,26,28,30,32,
                         34,36,38,40,50,
                         60,70,80,90,100,
                         120,140,160,180,200])

# Histograms for overall efficiency
mu_totalhist=ROOT.TH1D("total_events","Total Events",len(mu_bin_edges)-1,mu_bin_edges)
mu_filthist=ROOT.TH1D("filt_events","Filtered Events",len(mu_bin_edges)-1,mu_bin_edges)


# Split into three regions of eta
eta1_mu_totalhist=ROOT.TH1D("eta1_total_events","Total Events |eta| < 0.9",len(mu_bin_edges)-1,mu_bin_edges)
eta1_mu_filthist=ROOT.TH1D("eta1_filt_events","Filtered Events |eta| < 0.9",len(mu_bin_edges)-1,mu_bin_edges)
eta2_mu_totalhist=ROOT.TH1D("eta2_total_events","Total Events 0.9 < |eta| < 2.1",len(mu_bin_edges)-1,mu_bin_edges)
eta2_mu_filthist=ROOT.TH1D("eta2_filt_events","Filtered Events 0.9 < |eta| < 2.1",len(mu_bin_edges)-1,mu_bin_edges)
eta3_mu_totalhist=ROOT.TH1D("eta3_total_events","Total Events 2.1 < |eta| < 2.4",len(mu_bin_edges)-1,mu_bin_edges)
eta3_mu_filthist=ROOT.TH1D("eta3_filt_events","Filtered Events 2.1 < |eta| < 2.4",len(mu_bin_edges)-1,mu_bin_edges)

#Gets Relevent Events from file
def Events(f):
	evs=f['Events'].arrays(['HLT_IsoMu27',
                'HLT_Mu50',
                'Muon_pt',
                'Muon_eta',
                'Muon_dz',
                'Muon_dxy',
                'Muon_pfRelIso03_all',
                'Muon_pfRelIso03_chg',
                'Muon_mediumId'])
	return evs



def isHLTMatched(events, muons):
        trigObj = ak.zip({
            "pt": events.TrigObj_pt,
            "eta": events.TrigObj_eta,
            "phi": events.TrigObj_phi,
            "mass": 0.,
            "id": events.TrigObj_id,
            "filterBits": events.TrigObj_filterBits
        },with_name = "Momentum4D")
 
        trigObjSingleMu = trigObj[(abs(trigObj.id) == 13)
                                  & ((trigObj.filterBits & 2) 
                                  | (trigObj.filterBits & 8)
                                  | (trigObj.filterBits & 1024))]

        #Calculating DeltaR
        def deltaR2(eta1, phi1, eta2, phi2):
	    deta = eta1 - eta2
            dphi = phi1- phi2
            dphi = np.mod(dphi + np.pi, 2*np.pi) - np.pi
        return deta**2 + dphi**2

        toMatch1mu, trigObjSingleMu = ak.unzip(ak.cartesian([muons, trigObjSingleMu], axis = 1, nested = True))
        alldr2 = deltaR2(toMatch1mu.eta, toMatch1mu.phi, trigObjSingleMu.eta, trigObjSingleMu.phi)
        min_dr2 = ak.min(alldr2, axis = 2)
        match1El = ak.any(min_dr2 < 0.1, axis = 1)
        return match1mu


# Function for filling the histograms
def muon_hists(events,etas,hists):
	totalhist=hists[0]
	filthist=hists[1]
	eta_min=etas[0]
	eta_max=etas[1]
   
	# trigger
	triggerSingleMuon = (events["HLT_IsoMu27"] | events["HLT_Mu50"])

	# quality requirements for muons
	muon_quality_check = ((events["Muon_mediumId"])
		& (events["Muon_pt"] >= 10)
                & (np.abs(events["Muon_eta"]) < 2.4)
                & (np.abs(events["Muon_dz"]) <= 0.1)
                & (np.abs(events["Muon_dxy"]) <= 0.02)
                & (events["Muon_pfRelIso03_all"] < 0.25)
                & (events["Muon_pfRelIso03_chg"] < 0.25)

	# cut on eta
	eta_split=((np.abs(events["Muon_eta"]) >= eta_min) & (np.abs(events["Muon_eta"]) < eta_max ))

	# Select based on quality and eta
	muons = events[muon_quality_check & eta_split]
	
	
	#Fill histograms
	for muon in muons:
            totalhist.Fill(muon.pt)

        #HLT Matching
        hltMatched = isHLTMatched(events, muons)
        tr_muons = muons[hltMatched & triggerSingleMuon]

        #Fill filtered histograms with HLT matched muons
        for tr_muon in tr_muons:
            filthiest.Fill(tr_muon.pt)
        return 0


#go through each file in folder
for file_name in os.listdir(input_folder):
	if file_name.endswith(".root"): 
		input_file = os.path.join(input_folder, file_name)
                        print(f"Processing File: {file_name}")
		
		with uproot.open(input_file) as f:
			evs=Events(f)
			eta_split=[[0.0,2.4],[0.0,0.9],[0.9,2.1],[2.1,2.4]]
			eta_hists=[[mu_totalhist,mu_filthist],[eta1_mu_totalhist,eta1_mu_filthist],[eta2_mu_totalhist,eta2_mu_filthist],[eta3_mu_totalhist,eta3_mu_filthist]]
			for (etas,hists) in zip(eta_split,eta_hists):
				muon_hists(evs,etas,hists)
	



# Fills efficiency plots
eta1_effs=ROOT.TEfficiency(eta1_mu_filthist,eta1_mu_totalhist)
eta2_effs=ROOT.TEfficiency(eta2_mu_filthist,eta2_mu_totalhist)
eta3_effs=ROOT.TEfficiency(eta3_mu_filthist,eta3_mu_totalhist)
c1 = ROOT.TCanvas ("canvas","",800,600)

# Get overall Efficiency:
mu_eff=mu_filthist.Clone("Overall Efficiency")

# Creates Efficiency Plot w legend
eta1_effs.SetTitle("Muon Trigger Efficiency in bins of pT;Muon pT [GeV];Efficiency")
legend=ROOT.TLegend(0.5,0.1,0.9,0.4)
legend.AddEntry(eta1_effs,"|#eta|<0.9","l")
legend.AddEntry(eta2_effs,"0.9<|#eta|<2.1","l")
legend.AddEntry(eta3_effs,"2.1<|#eta|<2.4","l")
legend.SetTextColor(ROOT.kBlack)
legend.SetTextFont(42)
legend.SetTextSize(0.03)

# Draw plot	
eta1_effs.Draw()
eta2_effs.SetLineColor(ROOT.kRed)
eta2_effs.Draw("same")
eta3_effs.SetLineColor(ROOT.kBlue)
eta3_effs.Draw("same")
legend.Draw("same")
c1.Update()

# Creates ouput file if needed and saves plot to pdf
os.makedirs("MC_eff_outputs_combined", exist_ok = True)
c1.SaveAs("MC_eff_outputs_combined/Overall_Efficiency.pdf")
	
	
# Saves overall efficiency
root_file = ROOT.TFile(output_file,"UPDATE")
root_file.cd()
	
eff_dir=root_file.Get("Efficiencies")
if not eff_dir:
        eff_dir=root_file.mkdir("Efficiencies")
        eff_dir.cd()
mu_eff.Write()

root_file.Close()


print("Samples' Efficiency Complete")
