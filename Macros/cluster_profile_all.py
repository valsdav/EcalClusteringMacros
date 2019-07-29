import sys
sys.argv.append("-b")
import ROOT as R  
import os
from tqdm import tqdm
from array import array
from collections import defaultdict
import CMS_lumi
import tdrstyle
from math import sqrt

R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)

tdrstyle.setTDRStyle()
#CMS_lumi.lumiText = "Train: {}, eta ring interval: {}".format(train, ring.replace("_","-"))
CMS_lumi.lumiText = ""
CMS_lumi.extraText="   Simulation Preliminary"

def getSigma(inputfile, en,eta, ofile, outputdir):

    f = R.TFile(inputfile)
    tree = f.Get("recosimdumper/caloTree")
    
    Npfclusters = 0 

    # if eta < 2:
    cprof = R.TH2F("cprofile_en{}_eta{}".format(en,eta), "pfCluster profile",25,-12.5,12.5,25,-12.5,12.5)
    xproj = R.TH1F("cprofx_en{}_eta{}".format(en,eta), "",25,-12.5,12.5)
    yproj = R.TH1F("cprofy_en{}_eta{}".format(en,eta), "",25,-12.5,12.5)
    # else:
    #     cprof = R.TH2F("cprofile_en{}_eta{}".format(en,eta), "pfCluster profile",13,-6.5,6.5,13,-6.5,6.5)
    #     xproj = R.TH1F("cprofx_en{}_eta{}".format(en,eta), "",13,-6.5,6.5)
    #     yproj = R.TH1F("cprofy_en{}_eta{}".format(en,eta), "",13,-6.5,6.5)

    Npfclusters = 0
    Nmulticlusterpercalo = 0
    pbar = tqdm(total=tree.GetEntries())
    for event in tree:
        pbar.update()
        ncalo = event.caloParticle_pt.size()
        #print("...............................................")
        for icalo in range(ncalo):  
            #print("----")
            hasPfCluster = False
            maxhit_ieta = -999
            maxhit_iphi = -999
            maxhit_energy = -1
            clusters_N = []
            ihit = -1

            for ieta,iphi,clusterhits in zip(event.simHit_ieta[icalo],
                event.simHit_iphi[icalo], event.pfClusterHit_energy[icalo]):
                ihit += 1

                if clusterhits.size() >= 1:
                    hasPfCluster = True
                    rechit = 0.
                    for clhit in clusterhits:
                        #print(ieta, iphi, clhit.first, clhit.second )
                        # there should be only 1 cluster since we are not overlapping them
                        cluster_n = clhit.first
                        if cluster_n not in clusters_N: clusters_N.append(cluster_n)
                        rechit += clhit.second
                    if (rechit > maxhit_energy):
                        maxhit_ieta = ieta 
                        maxhit_iphi = iphi 
                        maxhit_energy = rechit
                
            if hasPfCluster: Npfclusters+=1
            
            for ieta,iphi,clusterhits in zip(event.simHit_ieta[icalo],
                event.simHit_iphi[icalo], event.pfClusterHit_energy[icalo]):

                if clusterhits.size() >= 1:
                    rechit = 0.
                    for clhit in clusterhits:
                        # there should be only 1 cluster since we are not overlapping them
                        cluster_n = clhit.first
                        if cluster_n not in clusters_N: clusters_N.append(cluster_n)
                        rechit += clhit.second
                    cprof.Fill(ieta-maxhit_ieta, iphi-maxhit_iphi, rechit)

            # Adding the rechits of noise for each pfCluster associated to the caloparticle
            for nc in clusters_N:
                for ieta, iphi, ien in zip(event.pfClusterHit_noCaloPart_ieta[nc], 
                                event.pfClusterHit_noCaloPart_iphi[nc], 
                                event.pfClusterHit_noCaloPart_energy[nc]):
                    cprof.Fill(ieta-maxhit_ieta, iphi-maxhit_iphi, ien)
                    print(ieta, iphi, ien)
        

    pbar.close()
        # end of loop on enets

    print("Tot pfClusters: ", Npfclusters)
    if Npfclusters == 0:
        print("No clusters")
        return None, None
    #print(hits)
    # for (ieta,iphi), n in hits.items():
    #     bin = cprof.FindBin(ieta,iphi)
    #     print(ieta, iphi, cprof.GetBinContent(bin))
    #     cprof.SetBinContent(bin, cprof.GetBinContent(bin) / n )
    
    for ix in range(cprof.GetNbinsX()+1):
        ixsum = 0.
        ixerr = 0.
        for iy in range(cprof.GetNbinsY()+1):
            bin = cprof.GetBin(ix, iy)
            cprof.SetBinContent(bin, cprof.GetBinContent(bin)/ Npfclusters)
            cprof.SetBinError(bin, cprof.GetBinError(bin)/ sqrt(Npfclusters))
            ixsum += cprof.GetBinContent(bin)
            ixerr += cprof.GetBinError(bin)**2
        xproj.SetBinContent(ix, ixsum)
        xproj.SetBinError(ix, sqrt(ixerr) )
 
    for iy in range(cprof.GetNbinsY()+1):
        iysum = 0.
        iyerr = 0.
        for ix in range(cprof.GetNbinsX()+1):
            bin = cprof.GetBin(ix, iy)
            iysum += cprof.GetBinContent(bin)
            iyerr += cprof.GetBinError(bin)**2
    
        yproj.SetBinContent(iy, iysum)
        yproj.SetBinError(iy, sqrt(iyerr))


    gx = R.TF1("gausx", "gaus")
    gy = R.TF1("gausy", "gaus")
    xproj.SetTitle(";ieta/ix;Mean energy [GeV]")
    yproj.SetTitle(";iphi/iy;Mean energy [GeV]")
    # Projection X and Y
    # xproj = cprof.ProjectionX()
    # yproj = cprof.ProjectionY()
        
    ofile.cd()
    cprof.Write()

    label = R.TPaveText(0.75, 0.75, 0.95, 0.8, "NB NDC" )
    text = label.AddText("{} GeV - eta {:.1f}".format(en,eta))
    text.SetLineColor(R.kBlack)
    label.SetFillColor(R.kWhite)
   

    c = R.TCanvas("c_en{}_et{}".format(en, eta), "", 900, 800)
    c.SetRightMargin(0.2)
    cprof.Draw("COLZ")
    cprof.SetTitle("Energy: {}, eta: {}".format(en,eta))
    cprof.GetZaxis().SetTitle("Mean Energy (GeV)")
    cprof.GetXaxis().SetTitle("ieta/ix")
    cprof.GetYaxis().SetTitle("iphi/iy")
    c.Update()
    CMS_lumi.lumiText = "{} GeV - eta {:.1f}".format(en,eta)
    CMS_lumi.CMS_lumi(c, 0, 0)
    c.SaveAs(outputdir+"/c_en{}_et{}.png".format(en, eta))
    c.SaveAs(outputdir+"/c_en{}_et{}.C".format(en, eta))

    CMS_lumi.lumiText = ""

    c2x = R.TCanvas("cprojx_en{}_et{}".format(en, eta), "", 800, 800)
    xproj.Fit("gausx", "WL")
    xproj.Draw()
    c2x.SetLogy()
    c2x.Update()
    CMS_lumi.CMS_lumi(c2x, 0, 0)
    label.Draw("same")
    c2x.SaveAs(outputdir+"/cprojx_en{}_et{}.png".format(en, eta))

    c2y = R.TCanvas("cprojy_en{}_et{}".format(en, eta), "", 800, 800)
    yproj.Fit("gausy", "WL")
    yproj.Draw()
    c2y.SetLogy()
    c2y.Update()
    label.Draw("same")
    CMS_lumi.CMS_lumi(c2y, 0, 0)
    c2y.SaveAs(outputdir+"/cprojy_en{}_et{}.png".format(en, eta))
       
    return gx.GetParameter("Sigma"), gy.GetParameter("Sigma"), xproj.GetRMS(), yproj.GetRMS()


if __name__ == "__main__":

    etas = [0.2,0.5, 1, 1.2, 1.8]
    ens = [5, 10, 30,50,75, 100]
    #etas = [1.2]
    #ens = [5]

    results= {} 
    results_rms = {}

    outputfile = R.TFile("output.root", "recreate")

    inputdir = sys.argv[1]
    outputdir = sys.argv[2]

    for eta in etas:
        for en in ens:
            print("Eta", eta, "en", en)
            file = inputdir+ "/cluster_en{:.1f}_eta{:.1f}.root".format(en,eta)
            if not os.path.exists(file):
                print("Not exists")
                results[(eta,en)] = (0,0)
                continue
            sx, sy, rmsx, rmsy = getSigma(file,en,eta, outputfile, outputdir)
            if sx != None and sy != None:
                results[(eta, en)] = sx, sy
                results_rms[(eta, en)] = rmsx, rmsy



    graphsx_eta = {}
    graphsy_eta = {}
    graphrmsx_eta = {}
    graphrmsy_eta = {}

    for eta in etas:
        x = []
        sigmax = []
        sigmay = []
        rmsx = []
        rmsy = []
        for en in ens:
            if (eta, en) in results:
                sx, sy = results[(eta,en)]
                x.append(en)
                sigmax.append(sx)
                sigmay.append(sy)
                rx, ry = results_rms[(eta,en)]
                rmsx.append(rx)
                rmsy.append(ry)
                print(sx, sy)
                
        grx = R.TGraph(len(x), array("d", x), array("d", sigmax))
        gry = R.TGraph(len(x), array("d", x), array("d", sigmay))
        gr_rmsx = R.TGraph(len(x), array("d", x), array("d", rmsx))
        gr_rmsy = R.TGraph(len(x), array("d", x), array("d", rmsy))

        graphsx_eta[eta] = grx
        graphsy_eta[eta] = gry
        graphrmsx_eta[eta] = gr_rmsx
        graphrmsy_eta[eta] = gr_rmsy


    outputfile.cd()
    cx = R.TCanvas("cx", "", 800, 600)
    legx = R.TLegend(0.2,0.63,0.36,0.9)
    mgx = R.TMultiGraph()
    for eta, g in graphsx_eta.items():
        mgx.Add(g)
        g.SetLineWidth(2)
        legx.AddEntry(g, "Eta: {}".format(eta), "lp")

    mgx.Draw("APL0 PMC PLC")
    mgx.SetTitle("#sigma ieta;En [GeV];#sigma ieta")
    legx.Draw("SAME")
    cx.Update()
    CMS_lumi.CMS_lumi(cx, 0, 0)
    cx.Draw()
    cx.SaveAs(outputdir+"/Sigma_ieta.png")
    cx.SaveAs(outputdir+"/Sigma_ieta.C")

    cy = R.TCanvas("cy", "", 800, 600)
    legy = R.TLegend(0.2,0.63,0.36,0.9)
    mgy = R.TMultiGraph()
    for eta, g in graphsy_eta.items():
        mgy.Add(g)
        g.SetLineWidth(2)
        legy.AddEntry(g, "Eta: {}".format(eta), "lp")

    mgy.Draw("APL0 PMC PLC")
    mgy.SetTitle("#sigma iphi;En [GeV];#sigma iphi")
    legy.Draw("SAME")
    cy.Update()
    CMS_lumi.CMS_lumi(cy, 0, 0)
    cy.Draw()
    cy.SaveAs(outputdir+"/Sigma_iphi.png")
    cy.SaveAs(outputdir+"/Sigma_iphi.C")
    
    cx.Write()
    cy.Write()
    mgx.Write()
    mgy.Write()


    ## RMS graphs
    cx_rms = R.TCanvas("cx_rms", "", 800, 600)
    legxrms = R.TLegend(0.2,0.63,0.36,0.9)
    mgx_rms = R.TMultiGraph()
    for eta, g in graphrmsx_eta.items():
        mgx_rms.Add(g)
        g.SetLineWidth(2)
        legxrms.AddEntry(g, "Eta: {}".format(eta), "lp")

    mgx_rms.Draw("APL0 PMC PLC")
    mgx_rms.SetTitle("#sigma ieta;En [GeV];#sigma ieta")
    legxrms.Draw("SAME")
    cx_rms.Update()
    CMS_lumi.CMS_lumi(cx_rms, 0, 0)
    cx_rms.Draw()
    cx_rms.SaveAs(outputdir+"/Sigma_ieta_rms.png")
    cx_rms.SaveAs(outputdir+"/Sigma_ieta_rms.C")

    cy_rms = R.TCanvas("cy_rms", "", 800, 600)
    legy_rms = R.TLegend(0.2,0.63,0.36,0.9)
    mgy_rms = R.TMultiGraph()
    for eta, g in graphrmsy_eta.items():
        mgy_rms.Add(g)
        g.SetLineWidth(2)
        legy_rms.AddEntry(g, "Eta: {}".format(eta), "lp")

    mgy_rms.Draw("APL0 PMC PLC")
    mgy_rms.SetTitle("#sigma iphi;En [GeV];#sigma iphi")
    legy_rms.Draw("SAME")
    cy_rms.Update()
    CMS_lumi.CMS_lumi(cy_rms, 0, 0)
    cy_rms.Draw()
    cy_rms.SaveAs(outputdir+"/Sigma_iphi_rms.png")
    cy_rms.SaveAs(outputdir+"/Sigma_iphi_rms.C")
    
    cx_rms.Write()
    cy_rms.Write()
    mgx_rms.Write()
    mgy_rms.Write()

    outputfile.Close()
    print(results)



