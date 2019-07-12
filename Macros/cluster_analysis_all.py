import ROOT as R 
import sys 
from collections import defaultdict

R.TH1.SetDefaultSumw2()
R.gStyle.SetOptFit(1111)

def getSigma(inputfile)

    f = R.TFile(inputfile)
    tree = f.Get("recosimdumper/caloTree")
    Npfclusters = 0 

    cprof = R.TH2F("cprofile", "pfCluster profile",7,-3.5,3.5,7,-3.5,3.5);

    Npfclusters = 0

    for event in tree:
        hasPfCluster = False
        maxhit_ieta = -999
        maxhit_iphi = -999
        maxhit_energy = -1

        for ieta,iphi,clusterhit, rechit, rechit_ismatched in zip(event.simHit_ieta,
            event.simHit_iphi, event.pfClusterHit_energy, 
            event.recHit_energy, event.pfRecHit_isMatched):

            if (clusterhit>0 and rechit_ismatched):
                hasPfCluster = True
                if (rechit > maxhit_energy):
                    maxhit_ieta = ieta 
                    maxhit_iphi = iphi 
                    maxhit_energy = rechit
            
        if hasPfCluster: Npfclusters+=1
        
        for ieta,iphi,clusterhit, rechit, rechit_ismatched in zip(event.simHit_ieta,
            event.simHit_iphi, event.pfClusterHit_energy, 
            event.recHit_energy, event.pfRecHit_isMatched):
            if (clusterhit>0 and rechit_ismatched):
                cprof.Fill(ieta-maxhit_ieta, iphi-maxhit_iphi, rechit)

        # end of loop on enets

    print("Tot pfClusters: ", Npfclusters)
    #print(hits)
    # for (ieta,iphi), n in hits.items():
    #     bin = cprof.FindBin(ieta,iphi)
    #     print(ieta, iphi, cprof.GetBinContent(bin))
    #     cprof.SetBinContent(bin, cprof.GetBinContent(bin) / n )

    for ix in range(cprof.GetNbinsX()+1):
        for iy in range(cprof.GetNbinsY()+1):
            bin = cprof.GetBin(ix, iy)
            cprof.SetBinContent(bin, cprof.GetBinContent(bin)/ Npfclusters)
        

    print("Fitting")
    f2 = R.TF2("xygaus","xygaus",-3.5,3.5,-3.5,3.5);
    f2.SetParameter(3, 0.5)
    f2.SetParameter(5, 0.5)
    f2.SetParLimits(3, -1,1)
    f2.SetParLimits(5, -1,1)
    results = cprof.Fit(f2,"S");
    print(results.Chi2())

    return f2.GetParameter("SigmaX", "SigmaY")


if __name__ == "__main__":

    etas = [0.2,0.5, 1, 1.2, 1.8, 2.5]
    ens = [5,10,30,50,100]

   results= []

   inputdir = sys.argv[1]
   
   for eta in etas:
       for en in ens:
           file = inputdir+ "/cluster_en{:.1f}_eta{:.1f}.root".format(en,eta)
           results[(eta,en)] = getSigma(file)

    print(results)



