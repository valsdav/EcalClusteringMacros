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
    ncalo = event.caloParticle_pt.size()
    calo_trueE = event.caloParticle_simEnergy
    
    # map (ieta,iphi,icalo, simhit):(iclu, clhit)
    xtal_calo = defaultdict(list)
    # map (ieta, iphi, iclu, clhit):(icalo, simhit)
    xtal_cluster = defaultdict(list)
    # map (ieta, iphi):(iclu, noisehit)
    xtal_cluster_noise = defaultdict(list)
  
    for icalo in range(ncalo):
        print "--- icalo: ", icalo
        print "ieta iphi simhit [ pfcluster index , pfcluster hit]"
        for i, (ieta, iphi, simhit, clhit) in enumerate(zip(event.simHit_ieta[icalo],  event.simHit_iphi[icalo],
                        event.simHit_energy[icalo], event.pfClusterHit_energy[icalo])):
            if clhit.size() > 0:
                print ieta, iphi, "{:.5f}".format(simhit), [(hit.first, '{:.5f}'.format(hit.second)) for hit in clhit] 
                for chit in clhit:
                    xtal_cluster[(ieta, iphi, chit.first, chit.second)].append((icalo, simhit))
                    xtal_calo[(ieta, iphi, icalo, simhit)].append((chit.first, chit.second))
        #
    # Check the noise hits (not overlapping with caloparticle sihimt)
    for nclus , (energys, ietas,iphis) in enumerate(zip(event.pfClusterHit_noCaloPart_energy,   
                    event.pfClusterHit_noCaloPart_ieta, event.pfClusterHit_noCaloPart_iphi)):
        
        #print nclus , [(ieta, iphi, en) for ieta,iphi, en in zip(energys, ietas, iphis)]
        for en, ieta, iphi in zip(energys, ietas, iphis):
            xtal_cluster_noise[(ieta, iphi)].append((nclus, en))

    print "XTAL_cluster"
    pprint(xtal_cluster)
    print "XTAL_cluster_noise" 
    pprint(xtal_cluster_noise)
    print "XTAL_calo"
    pprint(xtal_calo)
    

    
    # now the plotting
    mean_ieta_cl =   mean([ ieta for (ieta, iphi,_,_) in xtal_cluster.keys()]).round()
    mean_iphi_cl =   mean([ iphi for (ieta, iphi,_,_) in xtal_cluster.keys()]).round()
    print(mean_ieta_cl, mean_iphi_cl)

    hratio_gamma1 = R.TH1F("eratio_g1", "", 100, -1, 2)
    hxtal_cluster       = R.TH2F("xtal_cluster", "xtal_cluster", 41, -20.5, +20.5, 41, -20.5, 20.5)
    hxtal_calo          = R.TH2F("xtal_calo", "xtal_caloparticle", 41, -20.5, +20.5, 41, -20.5, 20.5)
    contours = array("d",[1,2,3,4])
    hxtal_calo.SetContour(len(contours), contours )
    hxtal_calo.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_calo.GetZaxis().SetNdivisions(3, False)
    hxtal_cluster.SetContour(len(contours), contours )
    hxtal_cluster.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_cluster.GetZaxis().SetNdivisions(3, False)
   
    for (ieta, iphi, iclu, clhit), calohits in xtal_cluster.items():
        # print(ieta, iphi, iclu, clhit)
        # print(calohits)
        hxtal_cluster.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, iclu+1)

    for (ieta, iphi, icalo, simhit), clhits in xtal_calo.items():
        print(ieta, iphi, icalo)
        hxtal_calo.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, icalo+1)


    c1 = R.TCanvas("c1", "", 800,800)
    hxtal_cluster.Draw("COLZ")
    c1.Draw()               
            
    c2 = R.TCanvas("c2", "", 800,800)
    hxtal_calo.Draw("COLZ")
    c2.Draw()               
            
    raw_input("next?")