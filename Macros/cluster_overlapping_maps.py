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
    calo_simE = event.caloParticle_simEnergy
    
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

    ####################################################
    # Analysis of the caloparticle, pfcluster association
    # - Associate to each pfCluster the calo with the greatest fraction of 
    #   its simEnergy. (save the ordered list of them by fraction)
    # - Save a list of pfcluster associated with the method before to 
    #   a caloparticle (there can be more than one)
    #######################################################
    all_calo_clusters = list(set(map( itemgetter(2), xtal_cluster.keys())))

    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo:
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            cluster_calo_fraction[clid][icalo] += (simhit / calo_simE[icalo])
    
    # # Filter out calo with less than 5% fraction 
    # clean_cluster_calo_fraction = defaultdict(dict)
    # for clid, calos in cluster_calo_fraction.items():
    #     for icalo, frac in calos.items():
    #         if frac > 0.05:
    #             clean_cluster_calo_fraction[clid][icalo] = frac

    # print(clean_cluster_calo_fraction)

    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in cluster_calo_fraction.items():
        # order caloparticle by fraction
        caloids =  list(sorted(calofrac.items(), key=itemgetter(1), reverse=True))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction 
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1), reverse=True)))
    

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
    hxtal_cluster_assoc = R.TH2F("xtal_cluster_assoc", "xtal_cluster_assoc", 41, -20.5, +20.5, 41, -20.5, 20.5)
    hxtal_calo          = R.TH2F("xtal_calo", "xtal_caloparticle", 41, -20.5, +20.5, 41, -20.5, 20.5)
    hxtal_calo_fraction1 = R.TH2F("xtal_calo_fraction1", "xtal_caloparticle_frac1", 41, -20.5, +20.5, 41, -20.5, 20.5)
    hxtal_calo_fraction2 = R.TH2F("xtal_calo_fraction2", "xtal_caloparticle_frac2", 41, -20.5, +20.5, 41, -20.5, 20.5)
    hxtal_calo_fractions = [hxtal_calo_fraction1,hxtal_calo_fraction2]    

    contours = array("d",[1,2,3,4])
    hxtal_calo.SetContour(len(contours), contours )
    hxtal_calo.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_calo.GetZaxis().SetNdivisions(3, False)
    hxtal_cluster.SetContour(len(contours), contours )
    hxtal_cluster.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_cluster.GetZaxis().SetNdivisions(3, False)
    hxtal_cluster_assoc.SetContour(len(contours), contours )
    hxtal_cluster_assoc.GetZaxis().SetRangeUser(contours[0], contours[-1])
    hxtal_cluster_assoc.GetZaxis().SetNdivisions(3, False)
   
    for (ieta, iphi, iclu, clhit), calohits in xtal_cluster.items():
        # print(ieta, iphi, iclu, clhit)
        # print(calohits)
        hxtal_cluster.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, iclu+1)
        hxtal_cluster_assoc.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, cluster_calo_assoc[iclu][0][0]+1)

    for (ieta, iphi, icalo, simhit), clhits in xtal_calo.items():
        print(ieta, iphi, icalo)
        hxtal_calo.Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, icalo+1)

    for icalo in range(2):
        for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
            for ical, simhit in caloinfo:
                if ical != icalo: continue
                hxtal_calo_fractions[icalo].Fill(ieta-mean_ieta_cl, iphi-mean_iphi_cl, simhit / calo_simE[icalo])
            

    c1 = R.TCanvas("c1", "", 800,800)
    hxtal_cluster.Draw("COLZ")
    c1.Draw()               
            
    c2 = R.TCanvas("c2", "", 800,800)
    hxtal_calo.Draw("COLZ")
    c2.Draw()

    c3 = R.TCanvas("c3", "", 800,800)
    hxtal_cluster_assoc.Draw("COLZ")
    c3.Draw()   


    c4 = R.TCanvas("c4", "", 800,800)
    hxtal_calo_fractions[0].Draw("COLZ")
    c4.SetLogz()
    hxtal_calo_fractions[0].GetZaxis().SetRangeUser(1e-5, 1)
    c4.Draw()   

    c5 = R.TCanvas("c5", "", 800,800)
    hxtal_calo_fractions[1].Draw("COLZ")
    c5.SetLogz()
    hxtal_calo_fractions[1].GetZaxis().SetRangeUser(1e-5, 1)
    c5.Draw()                        
            
    raw_input("next?")