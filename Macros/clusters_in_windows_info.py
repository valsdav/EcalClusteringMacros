from __future__ import print_function
import ROOT as R 
R.gROOT.SetBatch(True)
#R.PyConfig.IgnoreCommandLineOptions = True

import sys 
import os
from tqdm import tqdm
from collections import defaultdict
from math import cosh
from itertools import islice, chain
from numpy import mean
from operator import itemgetter, attrgetter
from array import array
from math import sqrt
from pprint import pprint as pp
import numpy as np
import argparse
import association_strategies

'''
This script analyse the overlapping of two caloparticles
'''

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputfile", type=str, help="inputfile", required=True)
parser.add_argument("-o","--outputdir", type=str, help="outputdir", required=True)
parser.add_argument("-n","--nevents", type=int,nargs="+", help="n events iterator", required=False)
parser.add_argument("-d","--debug", action="store_true",  help="debug", default=False)
args = parser.parse_args()

debug = args.debug

if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)


f = R.TFile(args.inputfile);
tree = f.Get("recosimdumper/caloTree")

#pbar = tqdm(total=tree.GetEntries())


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
        return deltaR

# maximum baffo
window_ieta = 10
window_iphi = 20

totevents = 0

energies_maps = []


for iev, event in enumerate(tree):
    totevents+=1
    #pbar.update()
    print ('---', iev)
    # if iev % 10 == 0: print(".",end="")

    # Branches
    pfCluster_energies = event.pfCluster_energy
    pfCluster_eta = event.pfCluster_eta
    pfCluster_ieta = event.pfCluster_ieta
    pfCluster_iphi = event.pfCluster_iphi
    pfCluster_iz = event.pfCluster_iz
    calo_ieta = event.caloParticle_ieta
    calo_iphi = event.caloParticle_iphi
    calo_iz = event.caloParticle_iz
    calo_simE = event.caloParticle_simEnergy
    
    rechits_E = event.recHit_energy
    rechits_ieta = event.simHit_ieta #using simHit ieta structure
    rechits_iphi = event.simHit_iphi
    rechits_iz = event.simHit_iz
    pfrechits_filter = event.pfRecHit_isMatched
    rechits_nocalo_E = event.pfRecHit_unMatched_energy
    rechits_nocalo_ieta = event.pfRecHit_unMatched_ieta
    rechits_nocalo_iphi = event.pfRecHit_unMatched_iphi
    rechits_nocalo_iz = event.pfRecHit_unMatched_iz


    # xtal_calo = map (ieta,iphi,icalo, simhit):(iclu, clhit)
    # xtal_cluster =  map (ieta, iphi, iclu, clhit):{icalo: simhit}
    # xtal_cluster_noise = list (ieta, iphi, iclu, noisehit)

    (cluster_calo, calo_cluster), (xtal_cluster, xtal_calo, xtal_cluster_noise) =   \
                association_strategies.get_association(event, "sim_fraction", cluster_type="pfCluster",debug=debug)

   
    print(">>>>> Cluster_calo association")
    pp(cluster_calo)
    print("\n\n>>> Calo cluster association")
    pp(calo_cluster)


    # 1) Look for highest energy cluster
    clenergies_ordered = sorted([ (ic , en) for ic, en in enumerate(pfCluster_energies)], key=itemgetter(1), reverse=True)
    print ("biggest cluster", clenergies_ordered)

    used_pfclusters = []

    # Now iterate over clusters in order of energies
    for iw, (cli, clen) in enumerate(clenergies_ordered):
        if abs(pfCluster_eta[cli] > 1.479): continue

        if cli in used_pfclusters: continue
        
        sc_histo = R.TH2F("supercl_map_iev{}_iw{}".format(iev, iw), "", window_iphi*2, 0, window_iphi*2, window_ieta*2, 0, window_ieta*2)
        print("------ WINDOW {} -- cluster {}".format(iw, cli))

        pfClHit_numpy_window = []
        energy_dict = defaultdict(float)
        pfclusters_in_window = []
        caloparticle_in_window = []

        w_ieta = (pfCluster_ieta[cli] - window_ieta, pfCluster_ieta[cli] + window_ieta)
        w_iphi = (pfCluster_iphi[cli] - window_iphi, pfCluster_iphi[cli] + window_iphi)
        #print(pfCluster_eta[cli], w_ieta, pfCluster_phi[cli], w_iphi)
        
        # Take all clhits in the box, 
        # take w
        clhits_in_window = filter( 
                 lambda ((ieta,iphi,icl,clhit), l): icl not in used_pfclusters and \
                         ieta>=w_ieta[0] and ieta < w_ieta[1] and iphi>=w_iphi[0] and iphi<w_iphi[1] ,
                    xtal_cluster.items() )

        clhits_noise_in_window = filter( 
                        lambda (ieta,iphi, icln, clnhit) : icln not in used_pfclusters and \
                            ieta>=w_ieta[0] and ieta < w_ieta[1] and iphi>=w_iphi[0] and iphi<w_iphi[1] ,
                 xtal_cluster_noise)

        # Fraction
        caloenergy_in_window = defaultdict(float)
        calofraction_in_window = defaultdict(float)          

        for (ieta,iphi, icl, iclhit), caloinfo in clhits_in_window:
            energy_dict[(ieta, iphi)] += iclhit
            pfclusters_in_window.append(icl)
            caloparticle_in_window.append(cluster_calo[icl][0][0])
            if debug: sc_histo.Fill(iphi - w_iphi[0], ieta- w_ieta[0], iclhit)

            # filling calo fraction
            # N.B there could be two cluster in the same xtal-- > fraction not exact
            for icalo , calosim in caloinfo.items():
                caloenergy_in_window[icalo] += calosim
                calofraction_in_window[icalo] += calosim / calo_simE[icalo]

        
        for (ieta,iphi, iclnoise, clnoisehit) in clhits_noise_in_window:
            pfclusters_in_window.append(iclnoise)
            energy_dict[(ieta, iphi)] += clnoisehit

        caloparticle_in_window = list(set(caloparticle_in_window))
        pfclusters_in_window = list(set(pfclusters_in_window))
        
        print("PfClusters in window:")
        print(pfclusters_in_window)
        print("All caloparticle in window (at least one xtal in window)")
        print(caloenergy_in_window.keys())
        print("Caloparticle associated to pfclusters in window (winning the score):")
        for icalo in caloparticle_in_window:
            print("\tcalo: {} | clusters: {}".format(icalo, calo_cluster[icalo]))

        print("-- Calo fractions in this windows:")
        for icalo in caloenergy_in_window.keys():
            print("\tcalo: {} | fraction: {} | energy: {} ".format(icalo, calofraction_in_window[icalo], caloenergy_in_window[icalo], ))
        
        print("-- Calo particle scores for pfcluster")
        for icl in pfclusters_in_window:
            if icl in cluster_calo:
                print("\tcluster: {} | calo_associated: True".format(icl))
                for calo, caloscore in cluster_calo[icl]:
                    print("\t\tcalo: {} | score: {}".format(calo, caloscore))
            else:
                print("\tcluster: {} | calo_associated: False".format(icl))
        
        if cli in cluster_calo:
            #main caloparticle
            # Windows with associated calo
            main_calo = cluster_calo[cli] 

        used_pfclusters += pfclusters_in_window

        if len(energy_dict)> 0 and len(caloparticle_in_window)>0:
            c = R.TCanvas("c_{}_{}".format(iev, iw))
            sc_histo.Draw("COLZ")
            c.SaveAs(args.outputdir+"/c_{}_{}.png".format(iev, iw))

    

                    
        