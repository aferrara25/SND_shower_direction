import csv
import ROOT

def read_csv_into_matrix(filename):
    matrix = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)
        for row in csv_reader:
            # Convert row elements to integers
            row = [float(item) for item in row]
            matrix.append(row)
    return matrix

filename = 'shower_info.csv'
matrix = read_csv_into_matrix(filename)
#print(matrix)
tag_corr = ROOT.TH2D("Paper_vs_Density_tagging", "Paper_vs_Density_tagging; paper_tagging; density_tagging", 7, -1.5, 5.5, 7, -1.5, 5.5)
en_dist = ROOT.TH1D("E_reco", "E_reco; E (GeV); entries", 50, 0, 300)
scifi_corr = ROOT.TH2D("Paper_Wall_vs_E_SciFi", "Paper_Wall_vs_E_SciFi; paper_tagging; SciFi_energy (GeV)", 5, 0.5, 5.5, 50, -50, 150)
us_corr = ROOT.TH2D("Paper_Wall_vs_E_US", "Paper_Wall_vs_E_US; paper_tagging; US_energy (GeV)", 5, 0.5, 5.5, 50, 0, 200)
ratio = ROOT.TH2D("Scifi/E_reco vs starting wall", "Scifi/E_reco vs starting wall; paper_tagging; Scifi/E_reco (a.u.)", 5, 0.5, 5.5, 50, -1, 1)

for row in matrix:
    tag_corr.Fill(row[2],row[7])
    en_dist.Fill(row[5])
    scifi_corr.Fill(row[2],row[3])
    us_corr.Fill(row[2],row[4])
    ratio.Fill(row[2],row[3]/row[5])

# Create a ROOT file
root_file = ROOT.TFile("plots_shower.root", "RECREATE")

canvases = [ROOT.TCanvas(f"c{i}",f"c{i}", 800, 600) for i in range(5)]

# Write the histogram to the ROOT file
canvases[0].cd()
tag_corr.Draw("colz")
ROOT.gStyle.SetOptStat("ne")
canvases[0].Write()

canvases[1].cd()
en_dist.Draw("hist")
ROOT.gStyle.SetOptStat("ne")
canvases[1].Write()

canvases[2].cd()
scifi_corr.Draw("colz")
ROOT.gStyle.SetOptStat("ne")
canvases[2].Write()

canvases[3].cd()
us_corr.Draw("colz")
ROOT.gStyle.SetOptStat("ne")
canvases[3].Write()

canvases[4].cd()
ratio.Draw("colz")
ROOT.gStyle.SetOptStat("ne")
canvases[4].Write()

# Close the ROOT file
root_file.Close()