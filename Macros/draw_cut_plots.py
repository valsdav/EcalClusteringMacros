import ROOT as R 
import argparse
R.gStyle.SetPalette(R.kLightTemperature)

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--plots", type=int, nargs="+", default=list(range(50)))
args = parser.parse_args()

file = R.TFile("output.root")

cuts_en = ["cut_100KeV_50ns","cut_500KeV_50ns", "cut_1MeV_50ns", "cut_5MeV_50ns",
            "cut_10MeV_50ns","cut_100MeV_50ns"]
cuts_time = [ "cut_0KeV_50ns", "cut_0KeV_30ns","cut_0KeV_20ns","cut_0KeV_15ns",
               "cut_0KeV_10ns", "cut_0KeV_5ns" ]
simhit_cuts = ["cut_0KeV","cut_100KeV", "cut_1MeV", "cut_10MeV", "cut_100MeV"]

histos = []
canvas = []

'''
0. PFCluster / GenParticle energy
1. PCalohit timing spectrum
2. PCalohit energy spectrum
3. Pcalohits / GenParticle energy (energy cut)
4. Pcalohits / GenParticle energy (time cut)
5. Pcalohits / Simhits energy (energy cut)
6. Pcalohits / Simhits energy (time cut)
7. Simhits / GenParticle energy
8. SimHit energy spectrum
9. Simhit / GenParticle energy (energy cut)
10. GenParticle energy spectrum

'''

outdir="/eos/user/d/dvalsecc/ECAL/EcalClustering/Cut_Pcalohit_plots/"


if 0 in args.plots:
    c_recogen = R.TCanvas("c_recogen")
    h = file.Get("reco_gen")
    h.SetTitle("PFCluster / GenParticle energy")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_recogen.Draw()
    canvas.append(c_recogen)

if 1 in args.plots:
    c_calotime = R.TCanvas("c_calotime")
    legt = R.TLegend(0.5, 0.5, 0.85,0.7)
    h = file.Get("pcalohit_time")
    h2 = file.Get("pcalohit_time_encut")
    h.SetTitle("PCaloHit timing")
    legt.AddEntry(h, "PCaloHit time - no cut")
    legt.AddEntry(h2, "PCaloHit time - > 1 MeV")
    histos.append(h)
    histos.append(h2)
    h.SetLineWidth(2)
    h2.SetLineWidth(2)
    h.GetXaxis().SetTitle("Time [ns]")
    h.Draw("HIST PLC")
    h2.Draw("HIST PLC SAME")
    legt.Draw("same")
    c_calotime.SetLogy()
    c_calotime.Draw()
    canvas.append(c_calotime)


if 2 in args.plots:
    c_caloenergy = R.TCanvas("c_caloenergy")
    h = file.Get("pcalohit_energy")
    h.SetTitle("PCaloHit energy")
    h.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_caloenergy.SetLogy()
    c_caloenergy.Draw()
    canvas.append(c_caloenergy)

    c_caloenergy_zoom = R.TCanvas("c_caloenergy_zoom")
    h = file.Get("pcalohit_energy_zoom")
    h.SetTitle("PCaloHit energy")
    h.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_caloenergy_zoom.SetLogy()
    c_caloenergy_zoom.Draw()
    canvas.append(c_caloenergy_zoom)

if 3 in args.plots:
    histos3 = []
    max3 = -1
    c1 = R.TCanvas("c_calo_gen_en","Energy PCalohit / GenParticle energy  (energy cut)")
    leg = R.TLegend(0.15, 0.5, 0.4,0.8)
    for cut in cuts_en:
        h = file.Get("calo_gen_"+ cut)
        h.SetTitle("Energy PCalohit / GenParticle energy  (energy cut)")
        histos3.append(h)
        if h.GetMaximum() > max3: max3 = h.GetMaximum()
        h.SetLineWidth(2)
        leg.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    for h in histos3: h.GetYaxis().SetRangeUser(0, max3*1.25)
    leg.Draw("same")
    c1.Draw()
    canvas.append(c1)

if 4 in args.plots:
    histos4 = []
    max4 = -1
    c2 = R.TCanvas("c_calo_gen_time", "Energy PCalohit / GenParticle energy  (time cut)")
    leg2 = R.TLegend(0.15, 0.5, 0.4,0.8)
    for cut in cuts_time:
        h = file.Get("calo_gen_"+ cut)
        h.SetTitle("Energy PCalohit / GenParticle energy  (time cut)")
        histos4.append(h)
        if h.GetMaximum() > max4: max4 = h.GetMaximum()
        h.SetLineWidth(2)
        leg2.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    for h in histos4: h.GetYaxis().SetRangeUser(0, max4*1.25)
    leg2.Draw("same")
    c2.Draw()
    canvas.append(c2)


if 5 in args.plots:
    histos5 = []
    max5 = -1
    c1_simhit = R.TCanvas("c_calo_simhit_en","Energy PCalohits / SimHits  (energy cut)")
    leg3 = R.TLegend(0.15, 0.5, 0.4,0.8)
    for cut in cuts_en:
        h = file.Get("calo_simhit_"+ cut)
        h.SetTitle("Energy PCalohits / SimHits  (energy cut)")
        histos5.append(h)
        if h.GetMaximum() > max5: max5 = h.GetMaximum()
        h.SetLineWidth(2)
        leg3.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    for h in histos5: h.GetYaxis().SetRangeUser(0, max5*1.25)
    leg3.Draw("same")
    c1_simhit.Draw()
    canvas.append(c1_simhit)

if 6 in args.plots:
    histos6 = []
    max6 = -1
    c2_simhit = R.TCanvas("c_calo_simhit_time","Energy PCalohits / SimHits  (time cut)")
    leg6 = R.TLegend(0.15, 0.5, 0.4,0.8)
    maxi = -1
    for cut in cuts_time:
        h = file.Get("calo_simhit_"+ cut)
        h.SetTitle("Energy PCalohits / SimHits  (time cut)")
        histos6.append(h)
        if h.GetMaximum() > max6: max6 = h.GetMaximum()
        h.SetLineWidth(2)
        leg6.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    for h in histos6: h.GetYaxis().SetRangeUser(0, max6*1.25)
    leg6.Draw("same")
    c2_simhit.Draw()
    canvas.append(c2_simhit)


if 7 in args.plots:
    c_simhitgen = R.TCanvas("c_simhitgen")
    h = file.Get("simhits_gen")
    h.SetTitle("Simhit energy / GenParticle energy")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_simhitgen.SetLogy()
    c_simhitgen.Draw()
    canvas.append(c_simhitgen)

if 8 in args.plots:
    c_simhiten = R.TCanvas("c_simhiten")
    h = file.Get("simhit_energy")
    h.SetTitle("Simhits energy")
    h.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_simhiten.SetLogy()
    c_simhiten.Draw()
    canvas.append(c_simhiten)


    c_simhiten_zoom = R.TCanvas("c_simhiten_zoom")
    h = file.Get("simhit_energy_zoom")
    h.SetTitle("Simhits energy")
    h.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_simhiten_zoom.SetLogy()
    c_simhiten_zoom.Draw()
    canvas.append(c_simhiten_zoom)

if 9 in args.plots:
    histos9 = []
    max9 = -1
    c_simhit_cut = R.TCanvas("c_simhitcut","Energy SimHits / GenParticle (energy cut)")
    leg9 = R.TLegend(0.15, 0.5, 0.4,0.8)
    maxi = -1
    for cut in simhit_cuts:
        h = file.Get("simhit_gen_"+ cut)
        h.SetTitle("Energy SimHits / GenParticle (energy cut)")
        histos9.append(h)
        if h.GetMaximum() > max9: max9 = h.GetMaximum()
        h.SetLineWidth(2)
        leg9.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    for h in histos9: h.GetYaxis().SetRangeUser(0, max9*1.25)
    leg9.Draw("same")
    c_simhit_cut.Draw()
    canvas.append(c_simhit_cut)


if 10 in args.plots:
    c_genenergy = R.TCanvas("c_genenergy")
    h = file.Get("gen_energy")
    h.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_genenergy.SetLogy()
    c_genenergy.Draw()
    canvas.append(c_genenergy)

if 11 in args.plots:
    # c_simhiten_cumul = R.TCanvas("c_simhiten_cumul")
    # h = file.Get("simhit_energy")
    # bins = h.GetNbinsX()
    # h.SetBinContent(bins, h.GetBinContent(bins)+ h.GetBinContent(bins+1))
    # h.Scale(1/h.Integral())
    # h2 = h.GetCumulative(False)
    # h2.SetTitle("SimHit energy > E")
    # h2.GetXaxis().SetTitle("Energy [GeV]")
    # histos.append(h2)
    # h2.SetLineWidth(2)
    # h2.Draw("HIST PLC")
    # c_simhiten_cumul.SetLogy()
    # c_simhiten_cumul.Draw()
    # canvas.append(c_simhiten_cumul)

    # c_simhiten_cumul_2 = R.TCanvas("c_simhiten_cumul_2")
    # h = file.Get("simhit_energy")
    # bins = h.GetNbinsX()
    # h.SetBinContent(bins, h.GetBinContent(bins)+ h.GetBinContent(bins+1))
    # h.Scale(1/h.Integral())
    # h2 = h.GetCumulative(True)
    # h2.SetTitle("SimHit energy < E")
    # histos.append(h2)
    # h2.SetLineWidth(2)
    # h2.GetXaxis().SetTitle("Energy [GeV]")
    # h2.Draw("HIST PLC")
    # #c_simhiten_cumul_2.SetLogy()
    # c_simhiten_cumul_2.Draw()
    # canvas.append(c_simhiten_cumul_2)

    c_simhiten_cumul_zoom = R.TCanvas("c_simhiten_cumul_zoom")
    hsimhit_energy_zoom = file.Get("simhit_energy_zoom")
    bins = hsimhit_energy_zoom.GetNbinsX()
    hsimhit_energy_zoom.SetBinContent(bins, 
                hsimhit_energy_zoom.GetBinContent(bins)+ hsimhit_energy_zoom.GetBinContent(bins+1))
    hsimhit_energy_zoom.SetBinContent(bins+1, 0)
    hsimhit_energy_zoom.Scale(1/hsimhit_energy_zoom.Integral())
    histos.append(hsimhit_energy_zoom)
    h3 = hsimhit_energy_zoom.GetCumulative(False)
    print(h3.GetBinContent(h3.FindBin(0.001)))
    h3.SetTitle("SimHit energy > E")
    h3.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h3)
    h3.SetLineWidth(2)
    h3.Draw("HIST PLC")
    c_simhiten_cumul_zoom.Draw()
    canvas.append(c_simhiten_cumul_zoom)

    c_simhiten_cumul_zoom2 = R.TCanvas("c_simhiten_cumul_zoom2")
    h3 = hsimhit_energy_zoom.GetCumulative(True)
    print(h3.GetBinContent(h3.FindBin(0.001)))
    h3.SetTitle("SimHit energy < E")
    h3.GetXaxis().SetTitle("Energy [GeV]")
    histos.append(h3)
    h3.SetLineWidth(2)
    h3.Draw("HIST PLC")
    #c_simhiten_cumul.SetLogy()
    c_simhiten_cumul_zoom2.Draw()
    canvas.append(c_simhiten_cumul_zoom2)



for c in canvas:
    c.SaveAs(outdir + c.GetName() +".png")
    c.SaveAs(outdir + c.GetName() +".C")
