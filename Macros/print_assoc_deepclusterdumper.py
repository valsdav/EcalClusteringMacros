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
parser.add_argument("-o","--outputfile", type=str, help="outputfile", required=True)
parser.add_argument("-n","--nevents", type=int,nargs="+", help="n events iterator", required=False)
parser.add_argument("-d","--debug", action="store_true",  help="debug", default=False)
parser.add_argument("-b","--batch", action="store_true",  help="batch", default=False)
parser.add_argument("-s","--strategy", type=str, help="Strategy to output")
args = parser.parse_args()

debug = args.debug

if args.batch: R.gROOT.SetBatch(True)

f = R.TFile(args.inputfile);
tree = f.Get("deepclusteringdumper/caloTree")
pbar = tqdm(total=tree.GetEntries())


if args.nevents and len(args.nevents) >= 1:
    nevent = args.nevents[0]
    if len(args.nevents) == 2:
        nevent2 = args.nevents[1]
    else:
        nevent2 = nevent+1
    tree = islice(tree, nevent, nevent2)


totevents = 0

out = open(args.outputfile, "w")

for iev, event in enumerate(tree):
    totevents+=1
    if not args.batch: pbar.update()
    if debug: print( '---', iev)

    cl_assoc = getattr(event, "caloParticle_pfCluster_{}_MatchedIndex".format(args.strategy))
   
    for calo, clusters in enumerate(cl_assoc):
        #print(calo, clusters)
        if clusters[0] != -1:
            out.write("events:{} | calo: {} | clusters:{}\n".format(iev, calo, sorted(clusters)))

out.close()