from __future__ import print_function
import ROOT as R 
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
from pprint import pprint
import argparse
R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)
R.gStyle.SetOptStat(0)

'''
This script analyse the pfclusters and caloparticle and perform
an association based on the caloparticle simEnergy fraction. 
'''

debug = True

f = R.TFile(sys.argv[1]);
tree = f.Get("recosimdumper/caloTree")


if len(sys.argv)>2:
    nevent = int(sys.argv[2])
    if len(sys.argv) > 3:
        nevent2 = int(sys.argv[3])
    else:
        nevent2 = nevent+1
    tree = islice(tree, nevent, nevent2)


nbadevents = {"noPfClusters":0, "noGamma1Cluster": 0, "noGamma2Cluster":0}


def DeltaR(phi1, eta1, phi2, eta2):
        dphi = phi1 - phi2
        if dphi > R.TMath.Pi(): dphi -= 2*R.TMath.Pi()
        if dphi < -R.TMath.Pi(): dphi += 2*R.TMath.Pi()
        deta = eta1 - eta2
        deltaR = (deta*deta) + (dphi*dphi)
        return deltaR

totevents = 0

for iev, event in enumerate(tree):
    totevents+=1
    if debug: print( '---', iev)

    ncalo = event.caloParticle_pt.size()
    calo_genE = event.caloParticle_energy
    calo_simE = event.caloParticle_simEnergy
    calo_eta  = event.caloParticle_eta
    calo_phi  = event.caloParticle_phi
    pfCluster_eta = event.pfCluster_eta
    pfCluster_phi = event.pfCluster_phi
    pfCluster_energy = event.pfCluster_energy

    # map (ieta,iphi,icalo, simhit):(iclu, clhit)
    xtal_calo = defaultdict(list)
    # map (ieta, iphi, iclu, clhit):(icalo, simhit)
    xtal_cluster = defaultdict(list)
    # map iclu: (energy, ieta, iphi)
    # xtal_cluster_noise = defaultdict(list)

    for icalo in range(ncalo):
        print( "--- icalo ", icalo)
        # analysis cluster hit to get the cluster number
        if debug: print ("ieta iphi simhit [ pfcluster index , pfcluster hit]")
        for i, (ieta, iphi, simhit, clhit) in enumerate(zip(event.simHit_ieta[icalo],  event.simHit_iphi[icalo],
                        event.simHit_energy[icalo], event.pfClusterHit_energy[icalo])):
            if clhit.size() > 0:
                if debug:   print( ieta, iphi, "{:.5f}".format(simhit), [(hit.first, '{:.5f}'.format(hit.second)) for hit in clhit] )
                for chit in clhit:
                    xtal_cluster[(ieta, iphi, chit.first, chit.second)].append((icalo, simhit))
                    xtal_calo[(ieta, iphi, icalo, simhit)].append((chit.first, chit.second))
            #else:
                #print ieta, iphi, "{:.5f}".format(simhit), "!!! No Cluster hits"


    ####################################################
    # Analysis on the clusterhits for each caloparticle
    # Associate to each pfCluster the calo with the greatest fraction of 
    # its simEnergy
    #######################################################
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo:
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            cluster_calo_fraction[clid][icalo] += (simhit / calo_simE[icalo])

    print(cluster_calo_fraction)
    cluster_calo_assoc = {}
    for clid, calofrac in cluster_calo_fraction.items():
        caloid = sorted(calofrac.items(), key=itemgetter(1), reverse=True)[0][0]
        cluster_calo_assoc[clid] = caloid

    print ("cluster-calo association: ", cluster_calo_assoc)
    

