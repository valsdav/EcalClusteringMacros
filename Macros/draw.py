import ROOT as R 
import argparse
R.gStyle.SetPalette(R.kLightTemperature)

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--plots", type=int, nargs="+", default=list(range(50)))
args = parser.parse_args()

file = R.TFile("output.root")

cuts_en = ["cut000MeV_50ns", "cut030MeV_50ns","cut050MeV_50ns", "cut100MeV_50ns", "cut200MeV_50ns"]
cuts_time = ["cut000MeV_50ns", "cut000MeV_30ns","cut000MeV_20ns","cut000MeV_15ns",
                "cut000MeV_10ns","cut000MeV_5ns" ]

histos = []
canvas = []

'''
1. PCalohit timing
2. PCalohit energy
3. Pcalohits / GenParticle energy (energy cut)
4. Pcalohits / GenParticle energy (time cut)
5. Pcalohits / Simhits energy (energy cut)
6. Pcalohits / Simhits energy (time cut)
7. Simhits / GenParticle energy

'''

outdir="/eos/user/d/dvalsecc/ECAL/EcalClustering/Cut_Pcalohit_plots/"


if 0 in args.plots:
    c_recogen = R.TCanvas("c_recogen")
    h = file.Get("reco_gen")
    h.SetTitle("super")
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
    canvas.append(c_caloenergy)

if 3 in args.plots:
    c1 = R.TCanvas("c_calo_gen_en","Energy PCalohit / GenParticle energy  (energy cut)")
    leg = R.TLegend(0.15, 0.5, 0.4,0.8)
    for cut in cuts_en:
        h = file.Get("calo_gen_"+ cut)
        h.SetTitle("Energy PCalohit / GenParticle energy  (energy cut)")
        histos.append(h)
        h.SetLineWidth(2)
        leg.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    leg.Draw("same")
    c1.Draw()
    canvas.append(c1)

if 4 in args.plots:
    c2 = R.TCanvas("c_calo_gen_time", "Energy PCalohit / GenParticle energy  (time cut)")
    leg2 = R.TLegend(0.15, 0.5, 0.4,0.8)
    for cut in cuts_time:
        h = file.Get("calo_gen_"+ cut)
        h.SetTitle("Energy PCalohit / GenParticle energy  (time cut)")
        histos.append(h)
        h.SetLineWidth(2)
        leg2.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")

    leg2.Draw("same")
    c2.Draw()
    canvas.append(c2)


if 5 in args.plots:
    c1_simhit = R.TCanvas("c_calo_simhit_en","Energy PCalohits / SimHits  (energy cut)")
    leg3 = R.TLegend(0.15, 0.5, 0.4,0.8)
    for cut in cuts_en:
        h = file.Get("calo_simhit_"+ cut)
        h.SetTitle("Energy PCalohits / SimHits  (energy cut)")
        histos.append(h)
        h.SetLineWidth(2)
        leg3.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
    leg3.Draw("same")
    c1_simhit.Draw()
    canvas.append(c1_simhit)

if 6 in args.plots:
    c2_simhit = R.TCanvas("c_calo_simhit_time","Energy PCalohits / SimHits  (time cut)")
    leg6 = R.TLegend(0.15, 0.5, 0.4,0.8)
    maxi = -1
    for cut in cuts_time:
        h = file.Get("calo_simhit_"+ cut)
        h.SetTitle("Energy PCalohits / SimHits  (time cut)")
        histos.append(h)
        if h.GetMaximum() > maxi: maxi = h.GetMaximum()
        h.SetLineWidth(2)
        leg6.AddEntry(h, cut, "l")
        h.Draw("HIST PLC SAME")
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


for c in canvas:
    c.SaveAs(outdir + c.GetName() +".png")
    c.SaveAs(outdir + c.GetName() +".C")
