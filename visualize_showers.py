import ROOT

def plot_histograms_root(root_files, output_root, runs, filen):

    # Create a ROOT file to save histograms
    output_file = ROOT.TFile(output_root, "RECREATE")

    # Loop through each root file
    for j, root_file in enumerate(root_files):
        file = ROOT.TFile(root_file)
        canvas = ROOT.TCanvas(f"Run_{runs[j]}_{filen[j]}", f"Run_{runs[j]}_{filen[j]}", 800, 600)

        # Get histograms from the file
        hist_names = ["NoCut_Shower_development_X", "NoCut_Shower_development_Y", "GuilCut_Shower_development_X", "GuilCut_Shower_development_Y"]

        # Plot histograms
        for i, hist_name in enumerate(hist_names):
            hist = file.Get(hist_name)
            hist.GetZaxis().SetTitle("qdc (a.u.)")
            if i % 4 == 0:
                canvas.Clear()
                canvas.Divide(2, 2)

            canvas.cd(i % 4 + 1)
            output_file.cd()
            hist.Draw("COLZ0")
            ROOT.gStyle.SetOptStat("ne")
            stats = hist.FindObject("stats")
            stats.SetX1NDC(0.1)
            stats.SetX2NDC(0.3)
            stats.SetY1NDC(0.0)
            stats.SetY2NDC(0.2)

            # Draw horizontal lines at specific Y-values
            y_values = [2, 3, 4, 5]  # Example Y-values
            if i == 0:
                lines = [ROOT.TLine(hist.GetXaxis().GetXmin(), y, hist.GetXaxis().GetXmax(), y) for y in y_values]
            for line in lines:
                line.SetLineColor(ROOT.kRed)
                line.Draw()

            if (i + 1) % 4 == 0 or i == len(hist_names) - 1:
                canvas.Write()


        file.Close()

    # Close the output ROOT file
    output_file.Close()

# List of root files to open
runs = [4752, 4809, 4815, 4819, 4976, 4992, 5013, 5056, 5099, 5120, 5132, 5152, 5171, 5180, 5239, 5389, 5981, 6050, 6069, 6250, 6252, 6268, 6279, 6286, 6290, 6295, 6296, 6296, 6568, 6590, 6596, 6610, 6640]
filen = [51, 6, 3, 49, 64, 10, 42, 101, 75, 99, 49, 40, 18, 30, 174, 50, 18, 177, 73, 41, 142, 38, 63, 178, 116, 119, 4, 9, 256, 49, 192, 59, 72]
if len(runs) != len(filen):
    print("Different lenght!")

root_files = [f"TI18_output/TI18_outputRun_{r}_{f}.root" for r, f in zip(runs, filen)]

# Output PDF file
output = "shower_display.root"

# Plot histograms and save to PDF
plot_histograms_root(root_files, output, runs, filen)
