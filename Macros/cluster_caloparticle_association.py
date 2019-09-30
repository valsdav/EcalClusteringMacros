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
    h_ncl = R.TH1F("ncluster_{}".format(strategy), "", 10, 0, 10)
    h_enfrac_nclu = R.TProfile("enfrac_ncluster_{}".format(strategy),"", 10, 0, 10)
    h2_ncl_eta = R.TH2F("ncl_eta_{}".format(strategy), "", 20, 0, 3,10, 1, 11)
    
    h_ereco_etrue.SetTitle("E reco/E true;E_reco/E_true;")
    h_ncl.SetTitle("N Clusters per caloparticle;N Clusters;")
    h_enfrac_nclu.SetTitle("Algo fraction for associated cluster;Nth associated cluster")
    h2_ncl_eta.SetTitle("N. clusters: eta gen (algo:{});eta calo;N. clusters".format(strategy))

    histos["ereco_etrue"][strategy] = h_ereco_etrue
    histos["ncl"][strategy] = h_ncl 
    histos["enfrac_nclu"][strategy] = h_enfrac_nclu
    histos2d["ncl_eta"][strategy] = h2_ncl_eta


for iev, event in enumerate(tree):
    totevents+=1
    pbar.update()
    if debug: print( '---', iev)

    ncalo = event.caloParticle_pt.size()
    calo_simE = event.caloParticle_simEnergy
    calo_eta = event.caloParticle_eta
    pfCluster_energy = event.pfCluster_energy
    
    assoc, _ = association_strategies.get_all_associations(event, cluster_type="pfCluster",debug=debug)
    if debug: print(assoc)
    
    printout = False
    for strategy, (cluster_calo_assoc, calo_cluster_assoc) in assoc.items():
        for calo, clusters in calo_cluster_assoc.items():
            trueE = calo_simE[calo]
            recoE = 0.
            for icl, clid in enumerate(clusters):
                recoE += pfCluster_energy[clid]
                # If it is associated the calo is the first fraction in the cluster_calo list
                histos["enfrac_nclu"][strategy].Fill(icl+1, cluster_calo_assoc[clid][0][1])
            
            histos["ereco_etrue"][strategy].Fill(recoE/trueE)
            histos["ncl"][strategy].Fill(len(clusters))
            histos2d["ncl_eta"][strategy].Fill( abs(calo_eta[calo]), len(clusters))

            if debug and len(clusters)> 5:
                printout = True
    
        if printout:
            print("\n>>>",iev, "|", strategy, "| Cluster-calo map\n" , cluster_calo_assoc)
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
    


