#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TProfile2D.h>
#include <TProfile.h>   
#include <TF2.h>
#include <TStyle.h>
#include "stdlib.h"


using  namespace std;


pair<float, float> getSigma(string inputfile){
    TFile * f = new TFile (inputfile.c_str());
    TTreeReader tree ("recosimdumper/caloTree", f);
    
    TTreeReaderValue<UInt_t> genParticle_id (tree,"genParticle_id");
    TTreeReaderValue<float> genParticle_energy(tree,"genParticle_energy");
    TTreeReaderValue<float> genParticle_pt(tree,"genParticle_pt");
    TTreeReaderValue<float> genParticle_eta(tree,"genParticle_eta");
    TTreeReaderValue<float> genParticle_phi(tree,"genParticle_phi");
    TTreeReaderValue<float> caloParticle_energy(tree,"caloParticle_energy");
    TTreeReaderValue<float> caloParticle_pt(tree,"caloParticle_pt");
    TTreeReaderValue<float> caloParticle_eta(tree,"caloParticle_eta");
    TTreeReaderValue<float> caloParticle_phi(tree,"caloParticle_phi");
    TTreeReaderValue<std::vector<float>> simHit_energy(tree,"simHit_energy");
    TTreeReaderValue<std::vector<float>> simHit_eta(tree,"simHit_eta");
    TTreeReaderValue<std::vector<float>> simHit_phi(tree,"simHit_phi");
    TTreeReaderValue<std::vector<int>> simHit_ieta(tree,"simHit_ieta");
    TTreeReaderValue<std::vector<int>> simHit_iphi(tree,"simHit_iphi");
    TTreeReaderValue<std::vector<int>> simHit_iz(tree,"simHit_iz");    
    TTreeReaderValue<std::vector<float>> recHit_energy(tree, "recHit_energy");
    TTreeReaderValue<std::vector<bool>> pfRecHit_isMatched(tree, "pfRecHit_isMatched");
    //TTreeReaderValue<std::vector<float>> pfRecHit_energy(tree, "pfRecHit_energy");
    TTreeReaderValue<std::vector<float>> pfClusterHit_energy(tree, "pfClusterHit_energy");
    TTreeReaderValue<std::vector<float>> pfCluster_energy(tree, "pfCluster_energy");
    TTreeReaderValue<std::vector<float>> pfCluster_eta(tree, "pfCluster_eta");
    TTreeReaderValue<std::vector<float>> pfCluster_phi(tree, "pfCluster_phi");
    // TTreeReaderValue<std::vector<float>> superClusterHit_energy(tree, "superClusterHit_energy");
    // TTreeReaderValue<std::vector<float>> superCluster_energy(tree, "superCluster_energy");
    // TTreeReaderValue<std::vector<float>> superCluster_eta(tree, "superCluster_eta");
    // TTreeReaderValue<std::vector<float>> superCluster_phi(tree, "superCluster_phi");
    TTreeReaderValue<std::map<int,int>>  map_simHit_pfCLuster(tree, "map_simHit_pfCluster"); 
    // TTreeReaderValue<std::map<int,int>>  map_simHit_superCLuster(tree, "map_simHit_superCLuster");  


    TH2F* cprof = new TH2F("cprofile", "pfCluster profile",19,-9.5,+9.5,19,-9.5,+9.5);

    vector<int> Nrechits; 

    int Npfclusters = 0;

    while(tree.Next()){
        cout << "========" <<endl;
        //Nrechits.push_back(count_if(recHit_energy->begin(), recHit_energy->end(), [](float f){return f>0;}));
        float maxhit_ieta = -1;
        float maxhit_iphi = -1;
        float maxhit_energy = -1;

        for (int is =0; is < simHit_ieta->size() ; is++){

            if (pfClusterHit_energy->at(is) >0){
                cout << pfClusterHit_energy->at(is) << endl;
                Npfclusters+=1;
                if (pfClusterHit_energy->at(is) > maxhit_energy){
                    maxhit_energy = pfClusterHit_energy->at(is);
                    maxhit_ieta = simHit_ieta->at(is);
                    maxhit_iphi = simHit_iphi->at(is);
                }
            }
        
        }

        //cout << maxhit_energy << " " <<maxhit_ieta << " " << maxhit_iphi <<endl;
        for (int is =0; is < simHit_ieta->size() ; is++){
            if (pfClusterHit_energy->at(is)>0){
                cprof->Fill(simHit_ieta->at(is)- maxhit_ieta, 
                            simHit_iphi->at(is)- maxhit_iphi, 
                            pfClusterHit_energy->at(is));
            }
        }



        
    }

    cout << "Tot pfClusters: "<<Npfclusters <<endl;
    for (int ix =0; ix < cprof->GetNbinsX()+1; ix++){
        for (int iy =0; iy < cprof->GetNbinsY()+1; iy++){
            int bin = cprof->GetBin(ix, iy);
            cprof->SetBinContent(bin, cprof->GetBinContent(bin)/ Npfclusters);
        }
    }

    cout << "Fitting"<<endl;
    auto f2 = new TF2("bigaus","bigaus",-9.5,9.5,-9.5,9.5);
    cprof->Fit("bigaus", "N");

    f->Close();
    return make_pair(f2->GetParameter("SigmaX"), f2->GetParameter("SigmaY"));

}



int cluster_profile_all(string inputdir){
    TH1::SetDefaultSumw2();
    gStyle->SetOptFit(1111);

    vector<string> etas = {"0.2", "0.5", "1", "1.2", "1.8", "2.5"}; 
    vector<string> energies = {"5", "10", "30", "50", "100"};

    map<pair<string,string>,pair<float,float>> sigmas;

    for (auto const & eta: etas){
        for (auto const & en: energies){
            char * file = "/cluster_en%s_eta%s.root";
            auto ss = getSigma(inputdir + to_string(sprintf(file, en.c_str(), eta.c_str() )));
            sigmas[make_pair(eta,en)] = ss;
            cout << "En: " << en <<" eta: "<< eta<< " Sigma_ieta: " << 
                ss.first << " Sigma_iphi: " << ss.second << endl;
        }
    }
  

    return 0;
}