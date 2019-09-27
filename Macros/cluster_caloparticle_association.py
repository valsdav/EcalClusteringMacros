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
    
    assoc = association_strategies.get_all_associations(event, debug)
    print(assoc)
    

    # print("cluster calo fractions")
    # for clind , caloinfo in cluster_calo_fraction.items():
    #     print(">> cluster ", clind)
    #     for caloind, frac in caloinfo.items():
    #         print("     > calo: ", caloind, "  fraction: {:.4f}".format(frac))
    

    # print ("cluster-calo association: ", cluster_calo_assoc)
    # print("\n\n\nCalo-clusters association", sorted_calo_cluster_assoc)
    

