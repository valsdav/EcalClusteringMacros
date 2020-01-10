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
R.gStyle.SetPalette(R.kViridis)

'''
This script analyse the pfclusters and caloparticle and perform
an association based on the caloparticle simEnergy fraction. 
'''

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputfile", type=str, help="inputfile", required=True)
parser.add_argument("-o","--outputdir", type=str, help="outputdir", required=True)
parser.add_argument("-n","--nevents", type=int,nargs="+", help="n events iterator", required=False)
parser.add_argument("-d","--debug", action="store_true",  help="debug", default=False)
parser.add_argument("-b","--batch", action="store_true",  help="debug", default=False)
args = parser.parse_args()

debug = args.debug

if args.batch: R.gROOT.SetBatch(True)

f = R.TFile(args.inputfile);
tree = f.Get("recosimdumper/caloTree")
pbar = tqdm(total=tree.GetEntries())

if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)

if args.nevents and len(args.nevents) >= 1:
    nevent = args.nevents[0]
    if len(args.nevents) == 2:
        nevent2 = args.nevents[1]
    else:
        nevent2 = nevent+1
    tree = islice(tree, nevent, nevent2)


def DeltaR(phi1, eta1, phi2, eta2):
        dphi = phi1 - phi2
        if dphi > R.TMath.Pi(): dphi -= 2*R.TMath.Pi()
        if dphi < -R.TMath.Pi(): dphi += 2*R.TMath.Pi()
        deta = eta1 - eta2
        deltaR = (deta*deta) + (dphi*dphi)
        return sqrt(deltaR)

def get_etaphi_clusters(clusters, xtal_cluster, xtal_cluster_noise):
    mean_eta = 0.
    mean_phi = 0.
    sum_w = 0.
    for (eta, phi, clid, clhit) in filter( lambda (ieta,iphi,clid,clhit): clid in clusters, xtal_cluster.keys()):
        mean_eta += eta*clhit 
        mean_phi += phi*clhit
        sum_w += clhit
    for (eta, phi, clid, clhit) in filter( lambda (ieta,iphi,clid,clhit): clid in clusters, xtal_cluster_noise):
        mean_eta += eta*clhit
        mean_phi += phi*clhit
        sum_w += clhit
    return mean_eta / sum_w, mean_phi / sum_w

totevents = 0

histos = defaultdict(dict)
histos2d = defaultdict(dict)
effplots = defaultdict(dict)

for strategy in association_strategies.strategies.keys():
    h_etacalo = R.TH1F("etacalo_{}".format(strategy), "", 40, -3, 3)
    h_etacluster = R.TH1F("etacluster_{}".format(strategy), "", 40, -3,3)
    h_encalo = R.TH1F("energy_calo_{}".format(strategy), "", 40, 0, 400)
    h_encluster = R.TH1F("energy_cluster_{}".format(strategy), "", 40, 0,400)
    h_ereco_etrue_tot = R.TH1F("ereco_etrue_tot_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue = R.TH1F("ereco_etrue_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue_EB_tot = R.TH1F("ereco_etrue__EB_tot_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue_EB = R.TH1F("ereco_etrue_EB_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue_EE_tot = R.TH1F("ereco_etrue_EE_tot_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue_EE = R.TH1F("ereco_etrue_EE_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue_EB_best = R.TH1F("ereco_etrue_EB_best_{}".format(strategy), "", 25, 0, 1.5)
    h_ereco_etrue_EE_best = R.TH1F("ereco_etrue_EE_best_{}".format(strategy), "", 25, 0, 1.5)

    h_nxtal_clu = R.TH1F("nxtals_clu_{}".format(strategy), "", 50, 0, 50)
    h_deltaR_recotrue = R.TH1F("deltaR_recotrue_{}".format(strategy), "", 30, 0, 2)
    h_deltaR_recotrue_tot = R.TH1F("deltaR_recotrue_tot_{}".format(strategy), "", 30, 0, 2)

    h_etaratio = R.TH1F("etaratio_{}".format(strategy), "", 40, -2, 2)
    h_etaratio_tot = R.TH1F("etaratio_tot_{}".format(strategy), "", 40, -2, 2)
    h_ncl = R.TH1F("ncluster_{}".format(strategy), "", 20, 0, 20)
    h2_ncl_eta = R.TH2F("ncl_eta_{}".format(strategy), "", 20, 0, 3,20, 1, 21)
    h2_phireco_phitrue = R.TH2F("phireco_phitrue{}".format(strategy),"", 40, -3.14,3.14, 40, -3.14,3.14)
    h2_phireco_phitrue_tot = R.TH2F("phireco_phitrue_tot_{}".format(strategy),"", 40, -3.14,3.14, 40, -3.14,3.14)
    h2_etareco_etatrue = R.TH2F("etareco_etatrue{}".format(strategy),"", 40, -3, 3, 40, -3, 3)
    h2_etaratio_etacalo = R.TH2F("etaratio_etacalo{}".format(strategy),"",40, -3, 3, 40, -2, 2)
    h2_etaratio_encalo =  R.TH2F("etaratio_ecalo{}".format(strategy),"", 40, 0,500, 40, -2,2)
    h2_enratio_etacalo = R.TH2F("enratio_etacalo{}".format(strategy),"", 40, -3,3, 40, 0,2)
    h2_enratio_encalo =   R.TH2F("enratio_ecalo{}".format(strategy),"",  40, 0,500,  40, 0,2)
    h2_deltaR_encalo =  R.TH2F("deltaR_ecalo{}".format(strategy),"", 40, 0,500, 40, 0,1.5, )
    h2_deltaR_etacalo =  R.TH2F("deltaR_etacalo{}".format(strategy),"", 40, -3,3, 40, 0,1.5, )

    h2_etareco_etatrue_tot = R.TH2F("etareco_etatrue_tot_{}".format(strategy),"", 40, -3, 3, 40, -3, 3, )
    h2_etaratio_etacalo_tot = R.TH2F("etaratio_etacalo_tot_{}".format(strategy),"", 40, -3, 3, 40, -2, 2 )
    h2_etaratio_encalo_tot =  R.TH2F("etaratio_ecalo_tot_{}".format(strategy),"",40, 0,500, 40, -2,2 )
    h2_enratio_etacalo_tot = R.TH2F("enratio_etacalo_tot_{}".format(strategy),"", 40, -3,3, 40, 0,2)
    h2_enratio_encalo_tot =   R.TH2F("enratio_ecalo_tot_{}".format(strategy),"",  40, 0,500, 40, 0,2)
    h2_deltaR_encalo_tot =  R.TH2F("deltaR_ecalo_tot{}".format(strategy),"", 40, 0,500, 15, 0,0.3, )
    h2_deltaR_etacalo_tot =  R.TH2F("deltaR_etacalo_tot{}".format(strategy),"", 40, -3,3, 15, 0,0.3 )

    e_eff_eta = R.TEfficiency("eff_eta_{}".format(strategy), "", 30, -3, 3)
    e_eff_pt_EB = R.TEfficiency("eff_pt_{}".format(strategy), "", 25, 1, 100)
    e_eff_pt_EE = R.TEfficiency("eff_pt_{}".format(strategy), "", 25, 1, 100)
    e2_eff_eta_pt = R.TEfficiency("eff_eta_pt_{}".format(strategy), "", 25, -3, 3, 5, 1, 100)

    #######################
    ### Titles
    
    h_etacalo.SetTitle("#eta of calo particles; #eta calo")
    h_etacluster.SetTitle("#eta of pfClusters; #eta pfCluster")
    h_encalo.SetTitle("Energy of calo particles; Energy calo (GeV)")
    h_encluster.SetTitle("Energy of pfClusters; Energy pfCluster (GeV)")
    h_ereco_etrue.SetTitle("E reco/E true;E_reco/E_true;")
    h_ereco_etrue_tot.SetTitle("E reco (all pfClusters)/E true;E_reco (all pfCl.)/E_true;")
    h_ereco_etrue_EB.SetTitle("E reco/E true (EB);E_reco/E_true;")
    h_ereco_etrue_EB_tot.SetTitle("E reco (all pfClusters)/E true (EB);E_reco (all pfCl.)/E_true;")
    h_ereco_etrue_EE.SetTitle("E reco/E true (EE);E_reco/E_true;")
    h_ereco_etrue_EE_tot.SetTitle("E reco (all pfClusters)/E true (EE);E_reco (all pfCl.)/E_true;")
    h_ereco_etrue_EB_best.SetTitle("E reco (best score pfCluster)/ E true (EB); E_reco (best pfCluster)/E_true;")
    h_ereco_etrue_EE_best.SetTitle("E reco (best score pfCluster)/ E true (EE); E_reco (best pfCluster)/E_true;")
    h_etaratio.SetTitle("|#eta reco/ #eta true|;|#eta_reco/#eta_true|")
    h_etaratio_tot.SetTitle("|Mean #eta reco / #eta true|;|Mean #eta_reco/#eta_true|")
    h_ncl.SetTitle("N Clusters per caloparticle;N Clusters;")
    h_nxtal_clu.SetTitle("N xtals in cluster; N. xtals")
    h_deltaR_recotrue.SetTitle("#Delta R Calo-pfCluster; #DeltaR")
    h_deltaR_recotrue_tot.SetTitle("#Delta R Calo- All associated pfClusters; #DeltaR")
    h2_ncl_eta.SetTitle("N. clusters: eta gen (algo:{});eta calo;N. clusters".format(strategy))
    
    h2_phireco_phitrue.SetTitle("Phi calo: Phi pfCluster (algo: {}); Phi calo; Phi pfCluster".format(strategy))
    h2_phireco_phitrue_tot.SetTitle("Phi calo: Mean Phi pfClusters (algo: {}); Phi calo; Mean Phi pfClusters".format(strategy))
    h2_etareco_etatrue.SetTitle("Eta calo:Eta pfCluster (algo:{});Eta calo;Eta pfCluster".format(strategy))
    h2_etaratio_etacalo.SetTitle("Eta calo/eta pfCluster (algo:{}): Eta calo;Eta calo;Eta calo / eta pfCluster".format(strategy))
    h2_etaratio_encalo.SetTitle("Eta calo/eta pfCluster (algo:{}): Energy calo (GeV);Energy calo [GeV];Eta calo / eta pfCluster".format(strategy))
    h2_enratio_etacalo.SetTitle("En calo/En pfCluster (algo:{}): Eta calo;Eta calo; En calo / En pfCluster".format(strategy))
    h2_enratio_encalo.SetTitle("En calo/En pfCluster (algo:{}): Energy calo (GeV);Energy calo [GeV];En calo / En pfCluster".format(strategy))
    h2_deltaR_encalo.SetTitle("#Delta R Calo-pfCluster (algo:{}): Energy calo;Energy calo (GeV);#Delta R calo-cluster".format(strategy))
    h2_deltaR_etacalo.SetTitle("#Delta R Calo-pfCluster (algo:{}): Eta calo;Eta calo;#Delta R calo-cluster".format(strategy))

    h2_etareco_etatrue_tot.SetTitle("Eta calo:Mean Eta pfClusters (algo:{});Eta calo;Mean Eta pfCluster".format(strategy))
    h2_etaratio_etacalo_tot.SetTitle("Eta calo/Mean Eta pfClusters (algo:{}):Eta calo; Eta calo;Eta calo /Mean Eta pfCluster".format(strategy))
    h2_etaratio_encalo_tot.SetTitle("Eta calo/Mean Eta pfClusters (algo:{}): Energy calo;Energy calo [GeV];Eta calo /Mean Eta pfCluster".format(strategy))
    h2_enratio_etacalo_tot.SetTitle("En calo/Tot En pfClusters (algo:{}): Eta calo;Eta calo;En calo / Tot En pfCluster".format(strategy))
    h2_enratio_encalo_tot.SetTitle("En calo/Tot En pfClusters (algo:{}): Energy calo;Energy calo [GeV];En calo /Tot En pfCluster".format(strategy))
    h2_deltaR_encalo_tot.SetTitle("#Delta R Calo- All pfClusters (algo:{}): Energy calo;Energy calo (GeV);#Delta R calo- all clusters".format(strategy))
    h2_deltaR_etacalo_tot.SetTitle("#Delta R Calo- All pfClusters (algo:{}): Eta calo;Eta calo (GeV);#Delta R calo- all clusters".format(strategy))

    e_eff_eta.SetTitle("Eff. pfCluster (algo:{});#eta calo".format(strategy))
    e_eff_pt_EB.SetTitle("Eff. pfCluster (algo:{}) (EB);Pt calo".format(strategy))
    e_eff_pt_EE.SetTitle("Eff. pfCluster (algo:{}) (EE);Pt calo".format(strategy))
    e2_eff_eta_pt.SetTitle("Eff. pfCluster (algo:{});#eta calo; Pt calo".format(strategy))
    
    histos["etacalo"][strategy] = h_etacalo
    histos["etacluster"][strategy] = h_etacluster
    histos["encalo"][strategy] = h_encalo
    histos["encluster"][strategy] = h_encluster
    histos["ereco_etrue"][strategy] = h_ereco_etrue
    histos["ereco_etrue_tot"][strategy] = h_ereco_etrue_tot
    histos["ereco_etrue_EB"][strategy] = h_ereco_etrue_EB
    histos["ereco_etrue_EB_tot"][strategy] = h_ereco_etrue_EB_tot
    histos["ereco_etrue_EB_best"][strategy] = h_ereco_etrue_EB_best
    histos["ereco_etrue_EE"][strategy] = h_ereco_etrue_EE
    histos["ereco_etrue_EE_tot"][strategy] = h_ereco_etrue_EE_tot
    histos["ereco_etrue_EE_best"][strategy] = h_ereco_etrue_EE_best
    histos["etaratio"][strategy] = h_etaratio
    histos["etaratio_tot"][strategy] = h_etaratio_tot
    histos["ncl"][strategy] = h_ncl 
    histos["deltaR_recotrue"][strategy] = h_deltaR_recotrue
    histos["deltaR_recotrue_tot"][strategy] = h_deltaR_recotrue_tot
    histos["nxtal_clu"][strategy] = h_nxtal_clu

    histos2d["ncl_eta"][strategy] = h2_ncl_eta

    histos2d["etareco_etatrue"][strategy] = h2_etareco_etatrue
    histos2d["phireco_phitrue"][strategy] = h2_phireco_phitrue
    histos2d["etaratio_etacalo"][strategy] = h2_etaratio_etacalo
    histos2d["etaratio_encalo"][strategy] = h2_etaratio_encalo
    histos2d["deltaR_encalo"][strategy] = h2_deltaR_encalo
    histos2d["enratio_etacalo"][strategy] = h2_enratio_etacalo
    histos2d["enratio_encalo"][strategy] = h2_enratio_encalo
    histos2d["deltaR_etacalo"][strategy] = h2_deltaR_etacalo

    histos2d["etareco_etatrue_tot"][strategy] = h2_etareco_etatrue_tot    
    histos2d["phireco_phitrue_tot"][strategy] = h2_phireco_phitrue_tot
    histos2d["etaratio_etacalo_tot"][strategy] = h2_etaratio_etacalo_tot
    histos2d["etaratio_encalo_tot"][strategy] = h2_etaratio_encalo_tot
    histos2d["deltaR_encalo_tot"][strategy] = h2_deltaR_encalo_tot
    histos2d["enratio_etacalo_tot"][strategy] = h2_enratio_etacalo_tot
    histos2d["enratio_encalo_tot"][strategy] = h2_enratio_encalo_tot
    histos2d["deltaR_etacalo_tot"][strategy] = h2_deltaR_etacalo_tot

    effplots["eff_eta"][strategy] = e_eff_eta
    effplots["eff_pt_EB"][strategy] = e_eff_pt_EB
    effplots["eff_pt_EE"][strategy] = e_eff_pt_EE
    effplots["eff2_eta_pt"][strategy] = e2_eff_eta_pt


for iev, event in enumerate(tree):
    totevents+=1
    if not args.batch: pbar.update()
    if debug: print( '---', iev)

    ncalo = event.caloParticle_simEnergy.size()
    calo_simE = event.caloParticle_simEnergy
    calo_pt = event.caloParticle_simPt
    calo_simE = event.caloParticle_simEnergy
    calo_eta = event.caloParticle_simEta
    calo_phi = event.caloParticle_simPhi
    calo_geneta = event.caloParticle_genEta
    calo_genphi = event.caloParticle_genPhi
    pfCluster_energy = event.pfCluster_energy
    pfCluster_eta = event.pfCluster_eta
    pfCluster_phi = event.pfCluster_phi
    gen_eta = event.genParticle_eta
    gen_phi = event.genParticle_phi
    gen_pt = event.genParticle_pt
    
    assoc, (xtal_cluster, xtal_calo, xtal_cluster_noise) = \
            association_strategies.get_all_associations(event, cluster_type="pfCluster",debug=debug)
    if debug: print(assoc)
    
    printout = False
    for strategy, (cluster_calo_assoc, calo_cluster_assoc) in assoc.items():
        has_clusters_calo = {}
        for calo, clusters in calo_cluster_assoc.items():
            trueE = calo_simE[calo]
            histos["encalo"][strategy].Fill(trueE)
            histos["etacalo"][strategy].Fill(calo_eta[calo])
            recoE = 0.  
            # Save energy of first cluster
            if abs(pfCluster_eta[clusters[0]]) > 1.5:
                histos["ereco_etrue_EE_best"][strategy].Fill(pfCluster_energy[clusters[0]]/trueE)
            else:
                histos["ereco_etrue_EB_best"][strategy].Fill(pfCluster_energy[clusters[0]]/trueE)

            # Check efficiency
            hasCluster = False
            
            for clid in clusters:
                deltaR_clcalo = DeltaR(calo_genphi[calo], calo_geneta[calo], pfCluster_phi[clid], pfCluster_eta[clid])
                recoE += pfCluster_energy[clid]
                etaratio = pfCluster_eta[clid] / calo_eta[calo]
                enratio = pfCluster_energy[clid] / calo_simE[calo]
                
                if enratio > 0.4 and enratio < 1.4:
                    hasCluster = True

                histos["nxtal_clu"][strategy].Fill(
                    len(filter(lambda (ieta,iphi,icl,clhit):icl==clid, chain(xtal_cluster.keys(), xtal_cluster_noise) ))
                    )
                histos["etacluster"][strategy].Fill(pfCluster_eta[clid])
                histos["encluster"][strategy].Fill(pfCluster_energy[clid])
                histos["deltaR_recotrue"][strategy].Fill(deltaR_clcalo)

                histos["ereco_etrue"][strategy].Fill(pfCluster_energy[clid]/trueE)
                if abs(pfCluster_eta[clid]) > 1.5:
                    histos["ereco_etrue_EE"][strategy].Fill(pfCluster_energy[clid]/trueE)
                else:
                    histos["ereco_etrue_EB"][strategy].Fill(pfCluster_energy[clid]/trueE)

                histos["etaratio"][strategy].Fill(etaratio)
                # plots by single cluster
                histos2d["etareco_etatrue"][strategy].Fill(calo_eta[calo],pfCluster_eta[clid])
                histos2d["phireco_phitrue"][strategy].Fill(calo_phi[calo],pfCluster_phi[clid])

                histos2d["etaratio_etacalo"][strategy].Fill(calo_eta[calo],etaratio)
                histos2d["etaratio_encalo"][strategy].Fill(trueE, etaratio)
                histos2d["enratio_etacalo"][strategy].Fill(calo_eta[calo],enratio)
                histos2d["enratio_encalo"][strategy].Fill( trueE, enratio)
                histos2d["deltaR_encalo"][strategy].Fill(calo_simE[calo], deltaR_clcalo)
                histos2d["deltaR_etacalo"][strategy].Fill(calo_eta[calo],deltaR_clcalo)

            # Save if the calo has clusters that pass the efficiency criteria
            has_clusters_calo[calo] = hasCluster

            # plot with tot cluter
            mean_eta_cls, mean_phi_cls = get_etaphi_clusters(clusters, xtal_cluster, xtal_cluster_noise)
            etaratio_tot = mean_eta_cls / calo_eta[calo]
            enratio_tot = recoE / trueE
            deltaR_cls_calo = DeltaR(calo_genphi[calo], calo_geneta[calo], mean_phi_cls, mean_eta_cls)
            histos2d["etareco_etatrue_tot"][strategy].Fill(calo_eta[calo],mean_eta_cls)
            histos2d["phireco_phitrue_tot"][strategy].Fill(calo_phi[calo], mean_phi_cls)
            histos2d["etaratio_etacalo_tot"][strategy].Fill(calo_eta[calo],etaratio_tot)
            histos2d["etaratio_encalo_tot"][strategy].Fill(trueE,etaratio_tot)
            histos2d["enratio_etacalo_tot"][strategy].Fill(calo_eta[calo],enratio_tot)
            histos2d["enratio_encalo_tot"][strategy].Fill(trueE, enratio_tot )
            histos2d["deltaR_encalo_tot"][strategy].Fill( calo_simE[calo], deltaR_cls_calo)
            histos2d["deltaR_etacalo_tot"][strategy].Fill(calo_eta[calo], deltaR_cls_calo )
            histos["etaratio_tot"][strategy].Fill(etaratio_tot)
            histos["ereco_etrue_tot"][strategy].Fill(enratio_tot)
            if abs(mean_eta_cls) > 1.5:
                histos["ereco_etrue_EE_tot"][strategy].Fill(enratio_tot)
            else:
                histos["ereco_etrue_EB_tot"][strategy].Fill(enratio_tot)
            histos["deltaR_recotrue_tot"][strategy].Fill(deltaR_cls_calo)
            histos["ncl"][strategy].Fill(len(clusters))
            histos2d["ncl_eta"][strategy].Fill( abs(calo_eta[calo]), len(clusters))

            if debug and len(clusters)> 5:
                printout = True

        # Cluster efficiency
        for ngen, (g_eta, g_phi) in enumerate(zip(gen_eta, gen_phi)):
            calogen = None
            for calo in calo_cluster_assoc.keys():
                #remove crack
                # comparing gen with caloparticle_gen info, should be the same
                dr = DeltaR(g_phi, g_eta, calo_genphi[calo], calo_geneta[calo])
                if dr < 1e-3:
                    calogen = calo
                    break
            if calogen != None: 
                # Saving caloParticle simInfo, not gen one
                effplots["eff_eta"][strategy].Fill(has_clusters_calo[calogen], calo_eta[calo])
                if abs(calo_eta[calo])< 1.492 and abs(calo_eta[calo]) > 1.442: continue
                if abs(calo_eta[calo])> 1.5:
                    effplots["eff_pt_EE"][strategy].Fill(has_clusters_calo[calogen], calo_pt[calo])
                else:
                    effplots["eff_pt_EB"][strategy].Fill(has_clusters_calo[calogen], calo_pt[calo])
                effplots["eff2_eta_pt"][strategy].Fill(has_clusters_calo[calogen], calo_eta[calo], calo_pt[calo])
            else:
                effplots["eff_eta"][strategy].Fill(False, g_eta)
                effplots["eff2_eta_pt"][strategy].Fill(False, g_eta, gen_pt[ngen])
                #print("Lost gen",strategy, g_eta)
                if abs(calo_eta[calo])< 1.492 and abs(calo_eta[calo]) > 1.442: continue
                if abs(g_eta) > 1.5:
                    effplots["eff_pt_EE"][strategy].Fill(False, gen_pt[ngen])
                else:
                    effplots["eff_pt_EB"][strategy].Fill(False, gen_pt[ngen])

                

        if printout:
            print("\n>>>",iev, "|", strategy, "| Cluster-calo map\n" , cluster_calo_assoc)
            print("Calo-cluster map\n" , calo_cluster_assoc)

cache = []
R.gStyle.SetOptStat(0)
for var, hs in histos.items():
    c1 = R.TCanvas("c_{}".format(var), var, 800, 800)
    leg = R.TLegend(0.2, 0.638, 0.53, 0.873)
    m = -1
    i = 0
    for strategy, h in hs.items():
        for i in range(h.GetNbinsX()+1):
            if h.GetBinContent(i) == 0:
                h.SetBinContent(i, 10e-4)
               
        h.Draw("hist same PLC")
        leg.AddEntry(h, strategy, "l")
        h.SetLineWidth(3)
        mm = h.GetMaximum()
        if mm > m: m = mm

    leg.Draw("same")

    for _, h in hs.items():
        h.GetYaxis().SetRangeUser(0.8, m*100)
    c1.SetLogy()
    c1.Draw()
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+"_log.png")
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+"_log.C")
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+"_log.root")

    c1.SetLogy(False)
    for _, h in hs.items():
        h.GetYaxis().SetRangeUser(0, m*1.2)
    c1.Draw()
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+".png")
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+".C")
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+".root")
   
    

for var, hs in histos2d.items():
    for strategy, h in hs.items():
        c1 = R.TCanvas("c_{}_{}".format(var, strategy), var, 800, 800)
        h.Draw("COLZ")
        cache.append(c1)
        c1.SetLogz()
        c1.Draw()
        c1.SaveAs(args.outputdir+"/"+c1.GetName()+".png")
        c1.SaveAs(args.outputdir+"/"+c1.GetName()+".C")
        c1.SaveAs(args.outputdir+"/"+c1.GetName()+".root")


colors= {
    "sim_fraction" : R.TColor.GetColorPalette(1), 
    "nxtals":R.TColor.GetColorPalette(50), 
    "deltaR":R.TColor.GetColorPalette(100), 
    "sim_fraction_min1":R.TColor.GetColorPalette(150),
    "sim_fraction_min3":R.TColor.GetColorPalette(190),
    "sim_rechit_diff": R.TColor.GetColorPalette(220),
    "sim_rechit_fractions": R.TColor.GetColorPalette(250),
}
for var, hs in effplots.items():
    c1 = R.TCanvas("c_{}".format(var), var, 800, 800)
    leg = R.TLegend(0.35, 0.14, 0.68, 0.40)
    i = 0

    for strategy, h in hs.items():
        h.SetMarkerStyle(10)
        h.SetLineWidth(2)
        h.SetMarkerSize(1)
        h.SetLineColor(colors[strategy])
        h.SetMarkerColor(colors[strategy])

        if i == 0:
            h.Draw("APL")
            h.SetTitle("Eff. pfCluster")
            i+=1
        else:
            h.Draw("same PL")
        leg.AddEntry(h, strategy, "lp")

    leg.Draw("same")
    c1.Draw()
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+".png")
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+".C")
    c1.SaveAs(args.outputdir+"/"+c1.GetName()+".root")