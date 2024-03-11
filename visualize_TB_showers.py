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
            y_values = [2, 3, 4]  # Example Y-values
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
runs = [100639, 100647, 100672, 100673, 100631]
filen = [10, 11, 8, 9, 9]
if len(runs) != len(filen):
    print("Different lenght!")

root_files = [f"output/TB_outputRun_{r}_{f}.root" for r, f in zip(runs, filen)]

# Output PDF file
output = "shower_display_TB.root"

# Plot histograms and save to PDF
plot_histograms_root(root_files, output, runs, filen)
