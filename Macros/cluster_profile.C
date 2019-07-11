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


using  namespace std;

bool analyse_PCalohits = false;
bool useRechits_ = true;
bool usePFRechits_ = true;
bool usePFCluster_ = true;
bool useSuperCluster_ = false;

int cluster_profile(char* inputfile){
    TFile * f = new TFile (inputfile);
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

    
    TH2F * cluster_map = new TH2F("cluster_map", "cluster_energy", 170, -85,+85, 360, 0, 360);  

    TH1F* rechitEnergy = new TH1F("rechit_energy", "rechit_energy", 300, 0, 50);
    TH1F* pfRechitEnergy = new TH1F("pfRechit_energy", "pfRechit_energy", 300, 0, 50);


    TH2F* cprof = new TH2F("cprofile", "pfCluster profile",21,-10.5,+10.5,21,-10.5,+10.5);

    vector<int> Nrechits; 

    while(tree.Next()){

        //Nrechits.push_back(count_if(recHit_energy->begin(), recHit_energy->end(), [](float f){return f>0;}));
        float maxhit_ieta = -1;
        float maxhit_iphi = -1;
        float maxhit_energy = -1;

        for (int is =0; is < simHit_ieta->size() ; is++){

            rechitEnergy->Fill(recHit_energy->at(is));
            bool isMatched = false;
            if (pfRecHit_isMatched->at(is)){
                isMatched= true;
                pfRechitEnergy->Fill(recHit_energy->at(is));
            }

            if (pfClusterHit_energy->at(is) >0){
                cluster_map->Fill(simHit_ieta->at(is), simHit_iphi->at(is));

                if (pfClusterHit_energy->at(is) > maxhit_energy){
                    maxhit_energy = pfClusterHit_energy->at(is);
                    maxhit_ieta = simHit_ieta->at(is);
                    maxhit_iphi = simHit_iphi->at(is);
                }
            }
        }

        cout << maxhit_energy << " " <<maxhit_ieta << " " << maxhit_iphi <<endl;
        for (int is =0; is < simHit_ieta->size() ; is++){
            if (pfClusterHit_energy->at(is)>0){
                cprof->Fill(simHit_ieta->at(is)- maxhit_ieta, 
                            simHit_iphi->at(is)-maxhit_iphi, 
                            pfClusterHit_energy->at(is));
            }
        }

        
    }

    cout << "Fitting"<<endl;
    cprof->Fit("bigaus");
    //cout << "N_rechits mean:" << std::accumulate(Nrechits.begin(), Nrechits.end(), 0.0) / Nrechits.size() << endl;

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    cluster_map->Draw("COLZ");
    c1->Draw();

    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    rechitEnergy->Draw();
    pfRechitEnergy->Draw("same");
    pfRechitEnergy->SetLineColor(kRed);
    c2->SetLogy();
    c2->Draw();

    TCanvas *c3= new TCanvas("c3", "", 800,800);
    cprof->Draw("Lego");
    c3->Draw();

    return 0;

}
