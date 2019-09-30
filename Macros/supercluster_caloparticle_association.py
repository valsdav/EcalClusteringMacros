from __future__ import print_function
import ROOT as R 
import sys 
import os
from tqdm import tqdm
from collections import defaultdict
from math import cosh
from itertools import islice, chain
from numpy import mean
from array import array
from math import sqrt
from pprint import pprint
import argparse
import association_strategies

R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)
#R.gStyle.SetOptStat(0)
R.gStyle.SetPalette(R.kLightTemperature)

'''
This script analyse the pfclusters and caloparticle and perform
an association based on the caloparticle simEnergy fraction. 
'''

debug = False

f = R.TFile(sys.argv[1]);
tree = f.Get("recosimdumper/caloTree")

pbar = tqdm(total=tree.GetEntries())

if len(sys.argv)>2:
    nevent = int(sys.argv[2])
    if len(sys.argv) > 3:
        nevent2 = int(sys.argv[3])
    else:
        nevent2 = nevent+1
    tree = islice(tree, nevent, nevent2)


def DeltaR(phi1, eta1, phi2, eta2):
        dphi = phi1 - phi2
        if dphi > R.TMath.Pi(): dphi -= 2*R.TMath.Pi()
        if dphi < -R.TMath.Pi(): dphi += 2*R.TMath.Pi()
        deta = eta1 - eta2
        deltaR = (deta*deta) + (dphi*dphi)
        return deltaR

totevents = 0

histos = defaultdict(dict)
histos2d = defaultdict(dict)

for strategy in association_strategies.strategies.keys():
    h_ereco_etrue = R.TH1F("ereco_etrue_{}".format(strategy), "", 50, 0, 2)
    h_etareco_etatrue = R.TH1F("etareco_etatrue{}".format(strategy),"", 60, -2, 2)
    h_ncl = R.TH1F("ncluster_{}".format(strategy), "", 10, 0, 10)
    h2_etareco_etatrue = R.TH2F("2D_etareco_etatrue{}".format(strategy),"", 40, -2.5, 2.5, 40, -2.5, 2.5)
    h2_etaratio_ncl = R.TH2F("etaratio_ncl{}".format(strategy),"", 40, -2, 2, 4, 1, 4)
    h2_etaratio_etacalo = R.TH2F("etaratio_etacalo{}".format(strategy),"", 40, 0, 2.5, 40, -2, 2)
    h2_etaratio_ecalo = R.TH2F("etaratio_ecalo{}".format(strategy),"", 40, 0,200, 40, -2, 2 )
    h_deltaR = R.TH1F("deltaR_{}".format(strategy),"", 40, 0, 1)
    
    h_ereco_etrue.SetTitle("E reco/E true;E_reco/E_true;")
    h_ncl.SetTitle("N SuperClusters per caloparticle;N SuperClusters;")
    h_etareco_etatrue.SetTitle("Eta SCluster / Eta calo;Eta SCluster / Eta calo")
    h2_etareco_etatrue.SetTitle("Eta calo:Eta SCluster (algo:{});Eta calo;Eta SCluster".format(strategy))
    h_deltaR.SetTitle("#DeltaR caloparticle - SCluster; #DeltaR")
    h2_etaratio_ncl.SetTitle("Eta calo/eta SCluster (algo:{}):  N. SClusters;Eta calo / eta SCluster;N. SClusters".format(strategy))
    h2_etaratio_etacalo.SetTitle("Eta calo/eta SCluster (algo:{}): Eta calo;Eta calo;Eta calo / eta SCluster".format(strategy))
    h2_etaratio_ecalo.SetTitle("Eta calo/eta SCluster (algo:{}): Energy calo;Energy calo [GeV];Eta calo / eta SCluster".format(strategy))
    
    histos["ereco_etrue"][strategy] = h_ereco_etrue
    histos["etareco_etatrue"][strategy] = h_etareco_etatrue
    histos["deltaR"][strategy] = h_deltaR
    histos["ncl"][strategy] = h_ncl
    histos2d["etareco_etatrue"][strategy] = h2_etareco_etatrue
    histos2d["etaratio_ncl"][strategy] = h2_etaratio_ncl
    histos2d["etaratio_etacalo"][strategy] = h2_etaratio_etacalo
    histos2d["etaratio_ecalo"][strategy] = h2_etaratio_ecalo


for iev, event in enumerate(tree):
    totevents+=1
    pbar.update()
    if debug: print( '---', iev)

    ncalo = event.caloParticle_pt.size()
    calo_simE = event.caloParticle_simEnergy
    calo_eta = event.caloParticle_eta
    calo_phi = event.caloParticle_phi
    sCluster_energy = event.superCluster_energy
    sCluster_eta = event.superCluster_eta
    sCluster_phi = event.superCluster_phi
    
    assoc, _ = association_strategies.get_all_associations(event, cluster_type="superCluster",debug=debug)
    if debug: print(assoc)
    
    printout = False
    for strategy, (cluster_calo_assoc, calo_cluster_assoc) in assoc.items():
        
        for calo, clusters in calo_cluster_assoc.items():
            trueE = calo_simE[calo]
            recoE = 0.
            
            for icl, clid in enumerate(clusters):
                recoE += sCluster_energy[clid]
            
            histos["ncl"][strategy].Fill(len(clusters))
            histos["ereco_etrue"][strategy].Fill(recoE/trueE)

            etadiff = sCluster_eta[clusters[0]]/ calo_eta[calo]
            etaratio = sCluster_eta[clusters[0]]/ calo_eta[calo]

            histos["etareco_etatrue"][strategy].Fill(etaratio)
            histos2d["etareco_etatrue"][strategy].Fill(calo_eta[calo], sCluster_eta[clusters[0]])
            histos["deltaR"][strategy].Fill(DeltaR(sCluster_phi[clusters[0]], sCluster_eta[clusters[0]],
                                                    calo_phi[calo], calo_eta[calo]))
            histos2d["etaratio_ncl"][strategy].Fill(etaratio, len(clusters))
            histos2d["etaratio_etacalo"][strategy].Fill(abs(calo_eta[calo]), etaratio)
            histos2d["etaratio_ecalo"][strategy].Fill(trueE, etaratio)
            
            if abs(1-etaratio)> 0.8 :
                printout = True
        
    if printout:
        print("------------------------------------")
        for strategy, (cluster_calo_assoc, calo_cluster_assoc) in assoc.items():
            print(">>>",iev, "|", strategy)
            print("Cluster-calo map\n" , cluster_calo_assoc)
            print("Calo-cluster map\n" , calo_cluster_assoc)


cache = []
for var, hs in histos.items():
    c1 = R.TCanvas("c_{}".format(var), var, 800, 800)
    leg = R.TLegend()
    for strategy, h in hs.items():
        h.Draw("hist same PLC ")
        leg.AddEntry(h, strategy, "l")
        h.SetLineWidth(3)
    leg.Draw("same")
    cache.append(leg)
    c1.Draw()
    cache.append(c1)
    

for var, hs in histos2d.items():
    for strategy, h in hs.items():
        c1 = R.TCanvas("c_{}_{}".format(var, strategy), var, 800, 800)
        h.Draw("COLZ")
        cache.append(c1)
        c1.SetLogz()
        c1.Draw()
    


