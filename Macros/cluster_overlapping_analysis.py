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
import argparse
R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)
R.gStyle.SetOptStat(0)

'''
This script analyse the overlapping of two caloparticles
'''

parser = argparse.ArgumentParser()
parser.add_argument("--energy", type=float, help="energy gamma1", required=True)
parser.add_argument("--eta", type=float, help="eta gamma1", required=True)
parser.add_argument("--maxR", type=float, help="maxR plot", required=True)
parser.add_argument("--inputfile", type=str, help="inputfile", required=True)
parser.add_argument("--outputdir", type=str, help="outputdir", required=True)
parser.add_argument("--nevents", type=int,nargs="+", help="n events iterator", required=False)
parser.add_argument("--debug", type=bool,  help="debug", default=False)
args = parser.parse_args()

debug = args.debug

if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)


f = R.TFile(args.inputfile);
tree = f.Get("recosimdumper/caloTree")

pbar = tqdm(total=tree.GetEntries())

calo1_en = args.energy

if args.nevents and len(args.nevents) >= 1:
    nevent = args.nevents[0]
    if len(args.nevents) == 2:
        nevent2 = args.nevents[1]
    else:
        nevent2 = nevent+1
    tree = islice(tree, nevent, nevent2)


#plots Pf- true
hgamma1_Eratio = R.TH1F("hgamma1_Eratio", "Epf/Ecalo #gamma1", 100, 0.7, 1.3)
hgamma1_Eratio_gen = R.TH1F("hgamma1_Eratio_gen", "Ecalo/Egen #gamma1", 100, 0.7, 1.3)
hgamma2_Eratio = R.TH1F("hgamma2_Eratio", "Epf/Ecalo #gamma2", 100, 0.7, 1.3)

hgamma1_dEta = R.TH1F("gamma1_deltaEta", "#Delta#eta #gamma1 (PF - true)", 50, -0.2,-0.2)
hgamma1_dPhi = R.TH1F("gamma1#DeltaPhi", "#Delta#phi #gamma1 (PF - true)", 50, -0.1,-0.1)
hgamma1_dR   = R.TH1F("gamma1#DeltaR", "#DeltaR #gamma1 (PF - true)", 30, 0,0.03)

#overlap plot
hgamma12_dEta = R.TH1F("gamma12#DeltaEta", "#Delta#eta PF #gamma1-2", 50, -1,1)
hgamma12_dPhi = R.TH1F("gamma12#DeltaPhi", "#Delta#phi PF #gamma1-2", 50, -1,1)
hgamma12_dR   = R.TH1F("gamma12#DeltaR", "#DeltaR PF #gamma1-2", 30, 0,0.1)
hscan_Egamma1 = R.TProfile2D("scan_Egamma1", "En PF #gamma1 - E true #gamma1",25, 0.5, 100, 25, 0, args.maxR)
hscan_Egamma2 = R.TProfile2D("scan_Egamma2", "En PF #gamma2 - E true #gamma2",25, 0.5, 100, 25, 0, args.maxR)
hbadevent_dR    = R.TH1F("hbadevent_dR", "#DeltaR #gamma1-2 bad enents",  30, 0,args.maxR)
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
    pbar.update()
    if debug: print '---', iev

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
                    
    ####################################################
    # Analysis on the clusterhits for each caloparticle
    #######################################################
    cluster_energies = {}
    calo_cluster_assoc = {}

    # the event is bad if one of the caloparticle
    # doesn't have any associated pfClusterHit
    good_event = True

    all_calo_clusters = set(map( itemgetter(2), xtal_cluster.keys()))
    
    # cluster ids associated with caloparticle
    clids =  list( 
        map( itemgetter(0),        
                chain.from_iterable( 
                    # chain clusters hits
                    map( lambda (k,clusters): clusters,
                        # select gamma1
                        filter( lambda (k, v): k[2] == gamma1, xtal_calo.items() )
        ) ) ) )

    if len(all_calo_clusters) == 0:
        # no calo clusters
        good_event = False
        nbadevents["noPfClusters"] +=1
    else:
        if len(clids) > 0:
            # get the cluster id with more hits associated with the caloparticle
            gamma1_iclu = max(clids,  key = clids.count  ) 
            all_calo_clusters.remove(gamma1_iclu)
            calo_cluster_assoc[gamma1] = gamma1_iclu
            # use the cluster raw energy saved by the dumper (in includes already noise
            cluster_energies[gamma1] = pfCluster_energy[gamma1_iclu]
        else:
            # no calohit for gamma1
            good_event = False
            nbadevents["noGamma1Cluster"] += 1 
            
        if len(all_calo_clusters) == 1: 
            gamma2_iclu = list(all_calo_clusters)[0]
            calo_cluster_assoc[gamma2] = gamma2_iclu
            # use the cluster raw energy saved by the dumper (in includes already noise
            cluster_energies[gamma2] = pfCluster_energy[gamma2_iclu]
        else:
            # No clusters remained, can be that there is only 1 cluster associated to both the 
            # caloparticles, or that one caloparticle is non associates 
            calo_cluster_assoc[gamma2] = -1
            cluster_energies[gamma2]  = 0
            good_event = False
            nbadevents["noGamma2Cluster"] += 1
            
    if not good_event: 
        hbadevent_dR.Fill(DeltaR(calo_phi[gamma2], calo_eta[gamma2],
                             calo_phi[gamma1], calo_eta[gamma1] ))
        continue

    if debug: 
        print "Gamma1) cluster energy: ", cluster_energies[gamma1], " sum of simhits calo: ", calo_simE[gamma1], \
            " calo truth: ", calo_genE[gamma1]
        print "Gamma2) cluster energy: ", cluster_energies[gamma2], " sum of simhits calo: ", calo_simE[gamma2], \
            " calo truth: ", calo_genE[gamma2]

    

    hgamma1_Eratio.Fill(cluster_energies[gamma1]/calo_simE[gamma1])
    hgamma1_Eratio_gen.Fill(calo_simE[gamma1]/calo_genE[gamma1])    

    # Gamma 1 pfCluster vs true caloparticle (simhit)
    hgamma1_dEta.Fill(pfCluster_eta[ calo_cluster_assoc[gamma1] ] - calo_eta[gamma1])
    hgamma1_dPhi.Fill(pfCluster_phi[ calo_cluster_assoc[gamma1] ] - calo_phi[gamma1])
    hgamma1_dR.Fill(DeltaR(pfCluster_phi[ calo_cluster_assoc[gamma1] ],
                             pfCluster_eta[ calo_cluster_assoc[gamma1] ] ,
                             calo_phi[gamma1], calo_eta[gamma1] ))

    

    hgamma2_Eratio.Fill(cluster_energies[gamma2]/calo_simE[gamma2])
    #Overlap plots
    deltaR_clusters = DeltaR(pfCluster_phi[ calo_cluster_assoc[gamma1] ],
                             pfCluster_eta[ calo_cluster_assoc[gamma1] ] ,
                             pfCluster_phi[ calo_cluster_assoc[gamma2] ],
                             pfCluster_eta[ calo_cluster_assoc[gamma2] ]  )
    
    hgamma12_dEta.Fill(pfCluster_eta[ calo_cluster_assoc[gamma1] ] - pfCluster_eta[ calo_cluster_assoc[gamma2] ])
    hgamma12_dPhi.Fill(pfCluster_phi[ calo_cluster_assoc[gamma1] ] - pfCluster_phi[ calo_cluster_assoc[gamma2] ])
    hgamma12_dR.Fill( deltaR_clusters)

    # plot: x = gamma2 true energy, y = deltaR clusters, Z = pf energy - true energy
    hscan_Egamma1.Fill(calo_simE[gamma2],deltaR_clusters, cluster_energies[gamma1] - calo_simE[gamma1] )
    hscan_Egamma2.Fill(calo_simE[gamma2],deltaR_clusters, cluster_energies[gamma2] - calo_simE[gamma2] )


c1 = R.TCanvas("c1", "", 800, 600)
hgamma1_Eratio.Draw("hist")
hgamma1_Eratio.GetXaxis().SetTitle("Epf /Ecalo")
c1.Draw()  
c1.SaveAs(args.outputdir+ "/Gamma1_Eratio_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c1.SaveAs(args.outputdir+ "/Gamma1_Eratio_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))


c2 = R.TCanvas("c2", "", 800, 600)
hgamma2_Eratio.Draw("hist")
hgamma2_Eratio.GetXaxis().SetTitle("Epf /Ecalo")
c2.Draw()
c2.SaveAs(args.outputdir+ "/Gamma2_Eratio_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c2.SaveAs(args.outputdir+ "/Gamma2_Eratio_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))

c3 = R.TCanvas("c3", "", 800, 600)
hgamma1_dEta.Draw("hist")
hgamma1_dEta.GetXaxis().SetTitle("#Delta#eta")
c3.SaveAs(args.outputdir+ "/Gamma1_dEta_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c3.SaveAs(args.outputdir+ "/Gamma1_dEta_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c3.Draw()  

c4 = R.TCanvas("c4", "", 800, 600)
hgamma1_dPhi.Draw("hist")
hgamma1_dPhi.GetXaxis().SetTitle("#Delta#phi")
c4.SaveAs(args.outputdir+ "/Gamma1_dPhi_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c4.SaveAs(args.outputdir+ "/Gamma1_dPhi_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c4.Draw()  

c5 = R.TCanvas("c5", "", 800, 600)
hgamma1_dR.Draw("hist")
hgamma1_dR.GetXaxis().SetTitle("#DeltaR")
c5.SaveAs(args.outputdir+ "/Gamma1_dR_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c5.SaveAs(args.outputdir+ "/Gamma1_dR_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c5.Draw()  

c8 = R.TCanvas("c8", "", 800, 600)
hgamma12_dEta.Draw("hist")
hgamma12_dEta.GetXaxis().SetTitle("#Delta#eta")
c8.SaveAs(args.outputdir+ "/Gamma1-2_dEta_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c8.SaveAs(args.outputdir+ "/Gamma1-2_dEta_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c8.Draw()  

c9 = R.TCanvas("c9", "", 800, 600)
hgamma12_dPhi.Draw("hist")
hgamma12_dPhi.GetXaxis().SetTitle("#Delta#phi")
c9.SaveAs(args.outputdir+ "/Gamma1-2_dPhi_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c9.SaveAs(args.outputdir+ "/Gamma1-2_dPhi_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c9.Draw()

c10 = R.TCanvas("c10", "", 800, 600)
hgamma12_dR.Draw("hist")
hgamma12_dR.GetXaxis().SetTitle("#DeltaR")
c10.SaveAs(args.outputdir+ "/Gamma1-2_dR_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c10.SaveAs(args.outputdir+ "/Gamma1-2_dR_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c10.Draw() 

c6 = R.TCanvas("c6", "", 1100, 800)
c6.SetRightMargin(0.15)
hscan_Egamma1.Draw("colz")
hscan_Egamma1.GetZaxis().SetRangeUser(-15,15)
hscan_Egamma1.SetTitle("En PF #gamma1 - E true #gamma1;E #gamma2;#DeltaR 1-2;#DeltaE (GeV)")
c6.SaveAs(args.outputdir+ "/Overlapscan_gamma1_dEta_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c6.SaveAs(args.outputdir+ "/Overlapscan_gamma1_dEta_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c6.Draw()  

c7 = R.TCanvas("c7", "", 1100, 800)
c7.SetRightMargin(0.15)
hscan_Egamma2.Draw("colz")
hscan_Egamma2.GetZaxis().SetRangeUser(-15,15)
hscan_Egamma2.SetTitle("En PF #gamma2 - E true #gamma2;E #gamma2;#DeltaR 1-2;#DeltaE (GeV)")
c7.SaveAs(args.outputdir+ "/Overlapscan_gamma2_dEta_E{:.1f}_eta{:.1f}.png".format(args.energy, args.eta))
c7.SaveAs(args.outputdir+ "/Overlapscan_gamma2_dEta_E{:.1f}_eta{:.1f}.C".format(args.energy, args.eta))
c7.Draw()   

# c8 = R.TCanvas("c8")
# hbadevent_dR.Draw("hist")
# c6.SaveAs(args.outputdir+ "/Overlapscan_gamma1_dEta_E{:1f}_eta{:1f}.png".format(args.energy, args.eta))
# c6.SaveAs(args.outputdir+ "/Overlapscan_gamma1_dEta_E{:1f}_eta{:1f}.C".format(args.energy, args.eta))
# c6.Draw()  
# c8.Draw()  

totbadevents = sum([v for v in nbadevents.values()])
print "Number of bad events: {} ({:.2f}%)".format(totbadevents, 100*totbadevents/totevents)
print "Number of bad events (no PfClusters): {} ({:.3f}%)".format(nbadevents["noPfClusters"], 
                                                        100*nbadevents["noPfClusters"]/totevents)
print "Number of bad events (no Gamma1 cluster): {} ({:.3f}%)".format(nbadevents["noGamma1Cluster"], 
                                                        100*nbadevents["noGamma1Cluster"]/totevents)
print "Number of bad events (no Gamma2 cluster): {} ({:.3f}%)".format(nbadevents["noGamma2Cluster"], 
                                                        100*nbadevents["noGamma2Cluster"]/totevents)
##################################################
# Snippets

# sum the cluster hit associated with the clusterassociated to the gammaX caloparticle
# assing also noise
# cluster_energies[gammaX]  =  sum( 
#         map(itemgetter(1), 
#             # filter cluster index
#             filter( lambda v: v[0] == gammaX_iclu, 
#                 chain.from_iterable( 
#                     # chain clusters hits
#                     map( lambda (k,clusters): clusters,
#                         # select gammaX
#                         filter( lambda (k, v): k[2] == gammaX, xtal_calo.items() )
#     ) ) ) ) 
#     )  + sum(
#             map(itemgetter(1), #summing the noise energy in that cluster 
#                 chain.from_iterable( 
#                     map(lambda (k, v): v,
#                         filter(lambda (k,v): k==gammaX_iclu, xtal_cluster_noise.items()) 
#     ) ) ) ) 