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


using  namespace std;



int cluster_profile(char* inputfile){
    TH1::SetDefaultSumw2();
    gStyle->SetOptFit(1111);

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

    
    TH2F * cluster_map_EB = new TH2F("cluster_map_EB", "cluster_energy_EB", 170, -85,+85, 360, 0, 360); 
    TH2F * cluster_map_EE = new TH2F("cluster_map_EE", "cluster_energy_EE", 100, 0, 100, 100, 0, 100);   

    TH1F* rechitEnergy_EB = new TH1F("rechit_energyEB", "rechit_energy EB", 300, 0, 50);
    TH1F* rechitEnergy_EE = new TH1F("rechit_energyEE", "rechit_energy EE", 300, 0, 50);

    TH1F* pfRechitEnergy_EB = new TH1F("pfRechit_energyEB", "pfRechit_energy EB", 300, 0, 50);
    TH1F* pfRechitEnergy_EE = new TH1F("pfRechit_energyEE", "pfRechit_energy EE", 300, 0, 50);

    TH2F* cprof = new TH2F("cprofile", "pfCluster profile",9,-4.5,+4.5,9,-4.5,+4.5);

    vector<int> Nrechits; 

    int Npfclusters = 0;

    while(tree.Next()){
        //Nrechits.push_back(count_if(recHit_energy->begin(), recHit_energy->end(), [](float f){return f>0;}));
        float maxhit_ieta = -1;
        float maxhit_iphi = -1;
        float maxhit_energy = -1;

        bool hasPfCluster = false;

        for (int is =0; is < simHit_ieta->size() ; is++){
            
            if (recHit_energy->at(is)>0){
                if (simHit_iz->at(is) == 0){
                    rechitEnergy_EB->Fill(recHit_energy->at(is));
                    if (pfRecHit_isMatched->at(is)){
                        pfRechitEnergy_EB->Fill(recHit_energy->at(is));
                    }
                    if(pfClusterHit_energy->at(is)>0){
                        cluster_map_EB->Fill(simHit_ieta->at(is), simHit_iphi->at(is));
                    }
                }else{
                    rechitEnergy_EE->Fill(recHit_energy->at(is));
                    
                    if (pfRecHit_isMatched->at(is)){
                        pfRechitEnergy_EE->Fill(recHit_energy->at(is));
                    }
                    if(pfClusterHit_energy->at(is)>0){
                        cluster_map_EE->Fill(simHit_ieta->at(is), simHit_iphi->at(is));
                    }
                }              
            }

            if (pfClusterHit_energy->at(is) >0){
                hasPfCluster = true;
                if (pfRecHit_isMatched->at(is) && recHit_energy->at(is) > maxhit_energy){
                    maxhit_energy = recHit_energy->at(is);
                    maxhit_ieta = simHit_ieta->at(is);
                    maxhit_iphi = simHit_iphi->at(is);
                }
            }
        
        }

        //cout << maxhit_energy << " " <<maxhit_ieta << " " << maxhit_iphi <<endl;
        for (int is =0; is < simHit_ieta->size() ; is++){
            if (pfClusterHit_energy->at(is)>0 && pfRecHit_isMatched->at(is)){
                cprof->Fill(simHit_ieta->at(is)- maxhit_ieta, 
                            simHit_iphi->at(is)- maxhit_iphi, 
                            recHit_energy->at(is));
            }
        }

        if (hasPfCluster) Npfclusters+=1;

        
    }

    cout << "Tot pfClusters: "<<Npfclusters <<endl;
    for (int ix =0; ix < cprof->GetNbinsX()+1; ix++){
        for (int iy =0; iy < cprof->GetNbinsY()+1; iy++){
            int bin = cprof->GetBin(ix, iy);
            cprof->SetBinContent(bin, cprof->GetBinContent(bin)/ Npfclusters);
        }
    }

    cout << "Fitting"<<endl;
    auto f2 = new TF2("bigaus","bigaus",-4.5,4.5,-4.5,4.5);
    cprof->Fit("bigaus");

    TCanvas *c3= new TCanvas("c3", "", 800,800);
    cprof->Draw("LEGO");
    c3->Draw();
    //cout << "N_rechits mean:" << std::accumulate(Nrechits.begin(), Nrechits.end(), 0.0) / Nrechits.size() << endl;

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    cluster_map_EB->Draw("COLZ");
    c1->Draw();

    // TCanvas *c1b = new TCanvas("c1b", "", 800, 600);
    // cluster_map_EE->Draw("COLZ");
    // c1b->Draw();

    TCanvas *c2 = new TCanvas("c2", "", 800, 600);
    rechitEnergy_EB->Draw("hist");
    pfRechitEnergy_EB->Draw("hist same");
    pfRechitEnergy_EB->SetLineColor(kRed);
    cout << "EB rechits: " << rechitEnergy_EB->Integral() << " pf EB rechits: " << pfRechitEnergy_EB->Integral() << endl;
    c2->SetLogy();
    c2->Draw();

    TCanvas *c4 = new TCanvas("c4", "", 800, 600);
    rechitEnergy_EE->Draw("hist");
    pfRechitEnergy_EE->Draw("hist same");
    pfRechitEnergy_EE->SetLineColor(kRed);
    c4->SetLogy();
    c4->Draw();



    return 0;

}
