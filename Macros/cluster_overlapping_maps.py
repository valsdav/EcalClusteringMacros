from __future__ import print_function
import ROOT as R 
import sys 
from tqdm import tqdm
from collections import defaultdict
from math import cosh
from itertools import islice
from numpy import mean
from operator import itemgetter
from pprint import pprint
from array import array
import association_strategies

R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)
R.gStyle.SetOptStat(0)

'''
This script plots a ieta:iphi map for cluster position and caloparticle position. 
It is useful for debugging purposes.
'''

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

for event in tree:
    print("---")
    pbar.update()
    calo_simE = event.caloParticle_simEnergy
    pfCluster_energy = event.pfCluster_energy
    
    # Get all associations
    assoc, (xtal_cluster, xtal_calo, xtal_cluster_noise) = association_strategies.get_all_associations(event, debug=True)


    all_calo_clusters = list(set(map( itemgetter(2), xtal_cluster.keys())))
   

    # now the plotting
    mean_ieta_cl =   mean([ ieta for (ieta, iphi,_,_) in xtal_cluster.keys()]).round()
    mean_iphi_cl =   mean([ iphi for (ieta, iphi,_,_) in xtal_cluster.keys()]).round()
    print(mean_ieta_cl, mean_iphi_cl)

    hratio_gamma1 = R.TH1F("eratio_g1", "", 100, -1, 2)
    hxtal_cluster       = R.TH2F("xtal_cluster", "xtal_cluster", 31, -15.5, +15.5, 31, -15.5, 15.5)
    hxtal_calo          = R.TH2F("xtal_calo", "xtal_caloparticle", 31, -15.5, +15.5, 31, -15.5, 15.5)
    contours = array("d",[1,2,3,4])
    hxtal_calo.SetContour(len(contours), contours )
    hxtal_calo.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_calo.GetZaxis().SetNdivisions(3, False)
    hxtal_cluster.SetContour(len(contours), contours )
    hxtal_cluster.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_cluster.GetZaxis().SetNdivisions(3, False)

    for (ieta, iphi, iclu, clhit), calohits in xtal_cluster.items():
        hxtal_cluster.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, iclu+1)
    for (ieta, iphi, icalo, simhit), clhits in xtal_calo.items():
            hxtal_calo.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, icalo+1)

    histos = {}
    for strategy, (cluster_calo_assoc, calo_cluster_assoc) in assoc.items():

        print("cluster calo fractions: >>>> ", strategy)
        for clind , caloinfo in cluster_calo_assoc.items():
            print(">> cluster ", clind)
            for caloind, frac in caloinfo:
                print("     > calo: ", caloind, "  fraction: {:.4f}".format(frac))

        print ("cluster-calo association: ", cluster_calo_assoc)
        print("calo-clusters association", calo_cluster_assoc)


        hxtal_calo_fraction1A = R.TH2F("calo0_cluster0_{}".format(strategy), "calo0_cluster0{}".format(strategy), 31, -15.5, +15.5, 31, -15.5, 15.5)
        hxtal_calo_fraction2A = R.TH2F("calo1_cluster0_{}".format(strategy), "calo1_cluster0{}".format(strategy), 31, -15.5, +15.5, 31, -15.5, 15.5)
        hxtal_calo_fraction1B = R.TH2F("calo0_cluster1_{}".format(strategy), "calo0_cluster1{}".format(strategy), 31, -15.5, +15.5, 31, -15.5, 15.5)
        hxtal_calo_fraction2B = R.TH2F("calo1_cluster0_{}".format(strategy), "calo1_cluster1{}".format(strategy), 31, -15.5, +15.5, 31, -15.5, 15.5)
        hxtal_calo_fractions = [[hxtal_calo_fraction1A,hxtal_calo_fraction1B],[hxtal_calo_fraction2A,hxtal_calo_fraction2B]]    
        histos[strategy] = hxtal_calo_fractions
    

        for icalo in range(2):
            for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
                for ical, simhit in caloinfo:
                    if ical != icalo: continue
                    if strategy == "sim_fraction":
                        hxtal_calo_fractions[icalo][clid].Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, simhit/calo_simE[icalo])
                    elif strategy == "sim_rechit_fractions":
                        hxtal_calo_fractions[icalo][clid].Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, 
                                                                abs( (simhit/calo_simE[icalo]) - (clhit/pfCluster_energy[clid])))
                    if strategy == "sim_rechit_global_fraction":
                        hxtal_calo_fractions[icalo][clid].Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, 
                                                               (simhit/calo_simE[icalo]) - (clhit/pfCluster_energy[clid]))
                

    c1 = R.TCanvas("c1", "", 800,800)
    hxtal_cluster.Draw("COLZ")
    c1.Draw()               
            
    c2 = R.TCanvas("c2", "", 800,800)
    hxtal_calo.Draw("COLZ")
    c2.Draw()

    for strategy, hs in histos.items():
        print(strategy, hs)
        c4 = R.TCanvas("c4_{}".format(strategy), strategy, 1000,1000)
        c4.Divide(2,2)
        c4.cd(1)
        hs[0][0].Draw("COLZ")
        c4.cd(2)
        hs[0][1].Draw("COLZ")
        c4.cd(3)
        hs[1][0].Draw("COLZ")
        c4.cd(4)
        hs[1][1].Draw("COLZ")
        #c4.SetLogz()
        if strategy== "sim_rechit_fractions":
            hs[0][0].GetZaxis().SetRangeUser(0, 0.1)
            hs[0][1].GetZaxis().SetRangeUser(0, 0.1)
            hs[1][0].GetZaxis().SetRangeUser(0, 0.1)
            hs[1][1].GetZaxis().SetRangeUser(0, 0.1)
        elif strategy == "sim_fraction":
            hs[0][0].GetZaxis().SetRangeUser(1e-5, 1)
            hs[0][1].GetZaxis().SetRangeUser(1e-5, 1)
            hs[1][0].GetZaxis().SetRangeUser(1e-5, 1)
            hs[1][1].GetZaxis().SetRangeUser(1e-5, 1)
            for i in range(1,5):
                c4.cd(i)
                c4.SetLogz()

        c4.cd()
        c4.Draw()   

                           
                        
            
    raw_input("next?")