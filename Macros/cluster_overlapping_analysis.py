import ROOT as R 
import sys 
from tqdm import tqdm
from collections import defaultdict
from math import cosh
from itertools import islice, chain
from numpy import mean
from operator import itemgetter, attrgetter
from pprint import pprint
from array import array

R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)

'''
This script plots a ieta:iphi map for cluster position and caloparticle position. 
It is useful for debugging purposes.
'''


f = R.TFile(sys.argv[1]);
tree = f.Get("recosimdumper/caloTree")

#pbar = tqdm(total=tree.GetEntries())

calo1_en = float(sys.argv[2])

if len(sys.argv)>3:
    nevent = int(sys.argv[3])
    if len(sys.argv) > 4:
        nevent2 = int(sys.argv[4])
    else:
        nevent2 = nevent+1
    tree = islice(tree, nevent, nevent2)


#plots
hgamma1_Eratio = R.TH1F("hgamma1_Eratio", "Epf/Etrue #gamma_1", 100, 0.7, 1.3)
hgamma2_Eratio = R.TH1F("hgamma2_Eratio", "Epf/Etrue #gamma_2", 100, 0.7, 1.3)

debug = False


for iev, event in enumerate(tree):
    print '---', iev
    #pbar.update()

    ncalo = event.caloParticle_pt.size()
    calo_genE = event.caloParticle_energy
    calo_simE = event.caloParticle_simEnergy
    calo_eta  = event.caloParticle_eta
    calo_phi  = event.caloParticle_phi

    # map (ieta,iphi,icalo, simhit):(iclu, clhit)
    xtal_calo = defaultdict(list)
    # map (ieta, iphi, iclu, clhit):(icalo, simhit)
    xtal_cluster = defaultdict(list)
    # map iclu: (energy, ieta, iphi)
    xtal_cluster_noise = defaultdict(list)
    
    nclus_hits = defaultdict(int)
    gamma1 = -1   #gamma 1 is fixed in the overlapping scan
    gamma2 = -1   #gamma 2 vary

    for icalo in range(ncalo):
        # Analysis gamma1, the fixed one
        if abs(calo_genE[icalo] - calo1_en) < 0.001:
            gamma1 = icalo
        else: 
            gamma2 = icalo

        # analysis cluster hit to get the cluster number
        if debug: print "ieta iphi simhit [ pfcluster index , pfcluster hit]"
        for i, (ieta, iphi, simhit, clhit) in enumerate(zip(event.simHit_ieta[icalo],  event.simHit_iphi[icalo],
                        event.simHit_energy[icalo], event.pfClusterHit_energy[icalo])):
            if clhit.size() > 0:
                if debug:   print ieta, iphi, "{:.5f}".format(simhit), [(hit.first, '{:.5f}'.format(hit.second)) for hit in clhit] 
                for chit in clhit:
                    xtal_cluster[(ieta, iphi, chit.first, chit.second)].append((icalo, simhit))
                    xtal_calo[(ieta, iphi, icalo, simhit)].append((chit.first, chit.second))
            #else:
                #print ieta, iphi, "{:.5f}".format(simhit), "!!! No Cluster hits"
                    
      
    for iclu in [0,1]: 
        if event.pfClusterHit_noCaloPart_energy[iclu].size() > 0:
            #print nclus , [(ieta, iphi, en) for ieta,iphi, en in zip(energys, ietas, iphis)]
            for en, ieta, iphi in zip(event.pfClusterHit_noCaloPart_energy[iclu],
                    event.pfClusterHit_noCaloPart_ieta[iclu], event.pfClusterHit_noCaloPart_iphi[iclu]):
                xtal_cluster_noise[iclu].append((en,ieta, iphi))


    ####################################################
    # Analysis on the clusterhits for each caloparticle
    #######################################################à
    cluster_energies = {}
    calo_cluster_assoc = {}

    for gammaX in [gamma1, gamma2]:
        # cluster ids associated with caloparticle
        clids =  list( 
            map( itemgetter(0),        
                    chain.from_iterable( 
                        # chain clusters hits
                        map( lambda (k,clusters): clusters,
                            # select gamma1
                            filter( lambda (k, v): k[2] == gammaX, xtal_calo.items() )
            ) ) ) )

        if len(clids) > 0:
            # get the cluster id with more hits associated with the caloparticle
            gammaX_iclu = max(clids,  key = clids.count  ) 
            calo_cluster_assoc[gammaX] = gammaX_iclu
                    
            # sum the cluster hit associated with the cluster with more hits in the gammaX calo xtals
            cluster_energies[gammaX]  =  sum( 
                    map(itemgetter(1), 
                        # filter cluster index
                        filter( lambda v: v[0] == gammaX_iclu, 
                            chain.from_iterable( 
                                # chain clusters hits
                                map( lambda (k,clusters): clusters,
                                    # select gammaX
                                    filter( lambda (k, v): k[2] == gammaX, xtal_calo.items() )
                ) ) ) ) 
                )  + sum(
                        map(itemgetter(1), #summing the noise energy in that cluster 
                            chain.from_iterable( 
                                map(lambda (k, v): v,
                                    filter(lambda (k,v): k==gammaX_iclu, xtal_cluster_noise.items()) 
                ) ) ) ) 
        else:
            # calohit associated
            calo_cluster_assoc[gammaX] = -1
            cluster_energies[gammaX]  = 0

    
    print("Gamma1) cluster energy: ", cluster_energies[gamma1], " sum of simhits calo: ", calo_simE[gamma1], 
        " calo truth: ", calo_genE[gamma1])
    hgamma1_Eratio.Fill(cluster_energies[gamma1]/calo_simE[gamma1])
    print("Gamma2) cluster energy: ", cluster_energies[gamma2], " sum of simhits calo: ", calo_simE[gamma2],
       " calo truth: ", calo_genE[gamma2])
    hgamma2_Eratio.Fill(cluster_energies[gamma2]/calo_simE[gamma2])


c1 = R.TCanvas("c1")
hgamma1_Eratio.Draw("hist")
c1.Draw()  
c2 = R.TCanvas("c2")
hgamma2_Eratio.Draw("hist")
c2.Draw()              
            
    # c2 = R.TCanvas("c2")
    # hxtal_calo.Draw("COLZ")
    # c2.Draw()               
            
    # raw_input("next?")