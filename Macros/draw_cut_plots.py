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
simhit_cuts = ["cut_100KeV", "cut_1MeV", "cut_10MeV", "cut_100MeV"]

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
    h = file.Get("pcalohit_time")
    h.SetTitle("PCaloHit timing")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_calotime.SetLogy()
    c_calotime.Draw()
    canvas.append(c_calotime)

if 2 in args.plots:
    c_caloenergy = R.TCanvas("c_caloenergy")
    h = file.Get("pcalohit_energy")
    h.SetTitle("PCaloHit energy")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_caloenergy.SetLogy()
    c_caloenergy.Draw()
    canvas.append(c_caloenergy)

    c_caloenergy_zoom = R.TCanvas("c_caloenergy_zoom")
    h = file.Get("pcalohit_energy_zoom")
    h.SetTitle("PCaloHit energy")
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
    c_simhiten = R.TCanvas("c_simhitem")
    h = file.Get("simhit_energy")
    h.SetTitle("Simhits energy")
    histos.append(h)
    h.SetLineWidth(2)
    h.Draw("HIST PLC")
    c_simhiten.SetLogy()
    c_simhiten.Draw()
    canvas.append(c_simhiten)


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


for c in canvas:
    c.SaveAs(outdir + c.GetName() +".png")
    c.SaveAs(outdir + c.GetName() +".C")