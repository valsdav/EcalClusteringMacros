import ROOT as R
import sys
from DataFormats.FWLite import Events, Handle
from collections import defaultdict
from itertools import islice
import numpy as np

events = Events (sys.argv[1])

handleEB  = Handle ('edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >' )
handleEE  = Handle ('edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >')

if len(sys.argv)>2:
    nevent = int(sys.argv[2])
    events = islice(events, nevent, nevent+1)

hs =[]

for iev, event in enumerate(events): 
    hb = R.TH2F("hitsEB_{}".format(iev), "hitsEB_{}".format(iev), 170, -85,+85,360,0, 360,)
    he = R.TH2F("hitsEE_{}".format(iev), "hitsEE_{}".format(iev),200, -100, 100, 100, 1, 100)
    hh = R.TH1F("hits_{}".format(iev), "hits_{}".format(iev), 300, 0, 3)
    htime = R.TH1F("hits_timing_{}".format(iev), "hits_timing_{}".format(iev),200, 0, 50) 	
    htime2 = R.TH2F("hits2_timing_{}".format(iev), "hits2_timing_{}".format(iev),200, 0, 50, 300, 0, 3) 	


    event.getByLabel( ("ecalRecHit" ,"EcalRecHitsEB"), handleEB)
    event.getByLabel( ("ecalRecHit" ,"EcalRecHitsEE"), handleEE)
    EBhits = handleEB.product()  
    EEhits = handleEE.product()

    xtals_energiesEB = defaultdict(float)
    xtals_energiesEE = defaultdict(float)
    xtals_times = defaultdict(list)

    for hit in EBhits:
        #print "ID: ", hit.id(), " | Track: ", hit.geantTrackId(), " | Energy: ", hit.energy()
        xtals_energiesEB[hit.id()] += hit.energy()
        xtals_times[hit.id()].append(hit.time())	
        htime.Fill(hit.time())

    for hit in EEhits:
        #print "ID: ", hit.id(), " | Track: ", hit.geantTrackId(), " | Energy: ", hit.energy()
        xtals_energiesEE[hit.id()] += hit.energy()    
        xtals_times[hit.id()].append(hit.time())	
        htime.Fill(hit.time())


    print "NXtals Barrel", len(xtals_energiesEB)
    for xtal, en in xtals_energiesEB.items():
        #print "XTAL: ", xtal, " | GeV: ", en
        detid =  R.EBDetId(xtal)
        hb.Fill( detid.ieta(), detid.iphi(), en)
        hh.Fill(en)


    print "NXtals Endcap ", len(xtals_energiesEE)
    for xtal, en in xtals_energiesEE.items():
        #print "XTAL: ", xtal, " | GeV: ", en
        detid =  R.EEDetId(xtal)
        he.Fill( detid.ix()*detid.zside(), detid.iy(), en)
        hh.Fill(en)

    for xtal, ts in xtals_times.items():
        if xtal in xtals_energiesEB:
            en = xtals_energiesEB[xtal] 
        else:
            en = xtals_energiesEE[xtal]
	    #print ts
        htime2.Fill(np.mean(ts), en)    


    print "Tot Energy: ", sum(xtals_energiesEB.values()) + sum(xtals_energiesEE.values())
    hs.append(hb)
    hs.append(he)
    hs.append(hh)
    hs.append(htime)    
    hs.append(htime2)

outfile = R.TFile("output_rechit_analysis.root", "recreate")
outfile.cd()
for h in hs:
    h.Write()
outfile.Close()
