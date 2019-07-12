import ROOT as R 
import sys 
import os
from array import array
from collections import defaultdict

R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)

def getSigma(inputfile, en,eta, ofile):

    f = R.TFile(inputfile)
    tree = f.Get("recosimdumper/caloTree")
    Npfclusters = 0 

    cprof = R.TH2F("cprofile_en{}_eta{}".format(en,eta), "pfCluster profile",7,-3.5,3.5,7,-3.5,3.5);

    Npfclusters = 0

    for event in tree:
        hasPfCluster = False
        maxhit_ieta = -999
        maxhit_iphi = -999
        maxhit_energy = -1

        for ieta,iphi,clusterhit, rechit, rechit_ismatched in zip(event.simHit_ieta,
            event.simHit_iphi, event.pfClusterHit_energy, 
            event.recHit_energy, event.pfRecHit_isMatched):

            if (clusterhit>0 and rechit_ismatched):
                hasPfCluster = True
                if (rechit > maxhit_energy):
                    maxhit_ieta = ieta 
                    maxhit_iphi = iphi 
                    maxhit_energy = rechit
            
        if hasPfCluster: Npfclusters+=1
        
        for ieta,iphi,clusterhit, rechit, rechit_ismatched in zip(event.simHit_ieta,
            event.simHit_iphi, event.pfClusterHit_energy, 
            event.recHit_energy, event.pfRecHit_isMatched):
            if (clusterhit>0 and rechit_ismatched):
                cprof.Fill(ieta-maxhit_ieta, iphi-maxhit_iphi, rechit)

        # end of loop on enets

    print("Tot pfClusters: ", Npfclusters)
    if Npfclusters == 0:
        print("No clusters")
        return -1,-1
    #print(hits)
    # for (ieta,iphi), n in hits.items():
    #     bin = cprof.FindBin(ieta,iphi)
    #     print(ieta, iphi, cprof.GetBinContent(bin))
    #     cprof.SetBinContent(bin, cprof.GetBinContent(bin) / n )

    for ix in range(cprof.GetNbinsX()+1):
        for iy in range(cprof.GetNbinsY()+1):
            bin = cprof.GetBin(ix, iy)
            cprof.SetBinContent(bin, cprof.GetBinContent(bin)/ Npfclusters)
        

    print("Fitting")
    f2 = R.TF2("xygaus","xygaus",-3.5,3.5,-3.5,3.5);
    f2.SetParameter(3, 0.5)
    f2.SetParameter(5, 0.5)
    f2.SetParLimits(3, -1,1)
    f2.SetParLimits(5, -1,1)
    results = cprof.Fit(f2,"SN");
    print(results.Chi2())

    ofile.cd()
    cprof.Write()
    f2.Write()
    

    return f2.GetParameter("SigmaX"), f2.GetParameter("SigmaY")


if __name__ == "__main__":

    etas = [0.2,0.5, 1, 1.2, 1.8, 2.5]
    ens = [30,50,100]

    results= {} 

    outputfile = R.TFile("output.root", "recreate")

    inputdir = sys.argv[1]

    for eta in etas:
        for en in ens:
            print("Eta", eta, "en", en)
            file = inputdir+ "/cluster_en{:.1f}_eta{:.1f}.root".format(en,eta)
            if not os.path.exists(file):
                print("Not exists")
                results[(eta,en)] = (0,0)
                continue
            results[(eta,en)] = getSigma(file,en,eta, outputfile)



    graphsx_eta = {}
    graphsy_eta = {}

    for eta in etas:
        x = []
        sigmax = []
        sigmay = []
        for en in ens:
            sx, sy = results[(eta,en)]
            x.append(en)
            sigmax.append(sx)
            sigmay.append(sy)
        
        grx = R.TGraph(len(x), array("d", x), array("d", sigmax))
        gry = R.TGraph(len(x), array("d", x), array("d", sigmay))

        graphsx_eta[eta] = grx
        graphsy_eta[eta] = gry


    cx = R.TCanvas("cx", "", 800, 600)
    legx = R.TLegend(0.8,0.8,0.95,0.95)
    mgx = R.TMultiGraph()
    for eta, g in graphsx_eta.items():
        mgx.Add(g)
        legx.AddEntry(g, "Eta: {}".format(eta), "lp")

    mgx.Draw("APL0 PMC PLC")
    legx.Draw("SAME")
    cx.Draw()

    cy = R.TCanvas("cy", "", 800, 600)
    legy = R.TLegend(0.8,0.8,0.95,0.95)
    mgy = R.TMultiGraph()
    for eta, g in graphsy_eta.items():
        mgy.Add(g)
        legy.AddEntry(g, "Eta: {}".format(eta), "lp")

    mgy.Draw("APL0 PMC PLC")
    legy.Draw("SAME")
    cy.Draw()
    
    outputfile.Close()
    print(results)



