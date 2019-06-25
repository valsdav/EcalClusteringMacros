
// Macro for producing plots for investigating the weights

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
//#include <boost/tokenizer.hpp>
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
#include "DataFormats/Math/interface/deltaR.h"


using  namespace std;

bool analyse_PCalohits = true;
bool useRechits_ = true;
bool usePFRechits_ = true;
bool usePFCluster_ = true;
bool useSuperCluster_ = false;

int analyze_gen_reco(){
    TFile * f = new TFile ("RecoSimDumper.root");
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
    // TTreeReaderValue<std::vector<float>> pfRecHit_energy(tree, "pfRecHit_energy");
    TTreeReaderValue<std::vector<float>> pfClusterHit_energy(tree, "pfClusterHit_energy");
    TTreeReaderValue<std::vector<float>> pfCluster_energy(tree, "pfCluster_energy");
    TTreeReaderValue<std::vector<float>> pfCluster_eta(tree, "pfCluster_eta");
    TTreeReaderValue<std::vector<float>> pfCluster_phi(tree, "pfCluster_phi");
    // TTreeReaderValue<std::vector<float>> superClusterHit_energy(tree, "superClusterHit_energy");
    // TTreeReaderValue<std::vector<float>> superCluster_energy(tree, "superCluster_energy");
    // TTreeReaderValue<std::vector<float>> superCluster_eta(tree, "superCluster_eta");
    // TTreeReaderValue<std::vector<float>> superCluster_phi(tree, "superCluster_phi");
    // TTreeReaderValue<std::map<int,int>>  map_simHit_pfCLuster(tree, "map_simHit_pfCLuster"); 
    // TTreeReaderValue<std::map<int,int>>  map_simHit_superCLuster(tree, "map_simHit_superCLuster");  

    TTreeReaderValue<std::vector<float>> caloHit_energy(tree,"caloHit_energy");       
    TTreeReaderValue<std::vector<float>> caloHit_time(tree,"caloHit_time");
    

    /*  The first check if E_reco/ E gen true: 
    *   Take the pfCluster with R < 0.05 from the gen particle and calculate
    *   the ratio of the energies, 
    */
    //TCanvas * c1 = new TCanvas("c1");
    //TCanvas * c2 = new TCanvas("c2");
    TFile* outputfile = new TFile("output.root", "RECREATE");
    TH1F* hreco_gen = new TH1F("reco_gen", "", 50, 0.6,1.1);
    TH1F* hsimhits_gen = new TH1F("simhits_gen", "", 100, 0.6,1.1);


    TH1F* pcalohit_energy = new TH1F("pcalohit_energy", "", 300, 0 , 30);
    TH1F* pcalohit_energy_zoom = new TH1F("pcalohit_energy_zoom", "", 100, 0 , 0.5);
    TH1F* pcalohit_time = new TH1F("pcalohit_time", "", 500, 0, 100);

    // pcalohit over genparticle energy
    map<pair<float, float>, TH1F*> hcalo_gen_histos;
    hcalo_gen_histos[make_pair(0.000, 50)] = new TH1F("calo_gen_cut000MeV_50ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.000, 100)] = new TH1F("calo_gen_cut000MeV_100ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.020, 50)] = new TH1F("calo_gen_cut020MeV_50ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.030, 50)] = new TH1F("calo_gen_cut030MeV_50ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.050, 50)] = new TH1F("calo_gen_cut050MeV_50ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.100, 50)] = new TH1F("calo_gen_cut100MeV_50ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.200, 50)] = new TH1F("calo_gen_cut200MeV_50ns", "", 100, 0.6,1.1);
    // for timing cuts
    hcalo_gen_histos[make_pair(0.000, 30)] = new TH1F("calo_gen_cut000MeV_30ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.000, 20)] = new TH1F("calo_gen_cut000MeV_20ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.000, 15)] = new TH1F("calo_gen_cut000MeV_15ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.000, 10)] = new TH1F("calo_gen_cut000MeV_10ns", "", 100, 0.6,1.1);
    hcalo_gen_histos[make_pair(0.000, 5)]  = new TH1F("calo_gen_cut000MeV_5ns" , "", 100, 0.6,1.1);

    // pcalohit over simhit total
    map<pair<float, float>, TH1F*> hcalo_simhit_histos;
    hcalo_simhit_histos[make_pair(0.000, 50)] = new TH1F("calo_simhit_cut000MeV_50ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.000, 100)] = new TH1F("calo_simhit_cut000MeV_100ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.020, 50)] = new TH1F("calo_simhit_cut020MeV_50ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.030, 50)] = new TH1F("calo_simhit_cut030MeV_50ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.050, 50)] = new TH1F("calo_simhit_cut050MeV_50ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.100, 50)] = new TH1F("calo_simhit_cut100MeV_50ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.200, 50)] = new TH1F("calo_simhit_cut200MeV_50ns", "", 100, 0.6,1.1);
    // for timing cuts
    hcalo_simhit_histos[make_pair(0.000, 30)] = new TH1F("calo_simhit_cut000MeV_30ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.000, 20)] = new TH1F("calo_simhit_cut000MeV_20ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.000, 15)] = new TH1F("calo_simhit_cut000MeV_15ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.000, 10)] = new TH1F("calo_simhit_cut000MeV_10ns", "", 100, 0.6,1.1);
    hcalo_simhit_histos[make_pair(0.000, 5)]  = new TH1F("calo_simhit_cut000MeV_5ns" , "", 100, 0.6,1.1);


    while(tree.Next()){
        //cout << "Event"<<endl;
        // Reco- gentrue check
        for (int icl =0; icl < pfCluster_eta->size(); icl++){
            double dR = deltaR(*genParticle_eta, *genParticle_phi, pfCluster_eta->at(icl), pfCluster_phi->at(icl));
            //cout  << "dR " << dR << " energy "<< pfCluster_energy->at(icl)<< endl;
            if ( dR < 0.05 ){
                hreco_gen->Fill(pfCluster_energy->at(icl) /  *genParticle_energy);
            }
        }   

        float tot_simhit_energy = 0.;
        for (auto & se : *simHit_energy){
            tot_simhit_energy += se;
        }
        hsimhits_gen->Fill(tot_simhit_energy / *genParticle_energy);

        //Save pcalohit info
        for (int ich=0; ich< caloHit_energy->size(); ich++){
                pcalohit_energy->Fill(caloHit_energy->at(ich));
                pcalohit_energy_zoom->Fill(caloHit_energy->at(ich));
                pcalohit_time->Fill(caloHit_time->at(ich));
        }

        // PCalo hit  over genparticle energy
        for (auto const& [cut, histo] : hcalo_gen_histos){
            float tot_energy_cut = 0.;
            for (int ich=0; ich< caloHit_energy->size(); ich++){
                //cout << "caloHit energy: "<< caloHit_energy->at(ich) << " | time: "<< caloHit_time->at(ich) <<endl;
                if (caloHit_energy->at(ich) > cut.first && caloHit_time->at(ich) < cut.second){
                    tot_energy_cut+= caloHit_energy->at(ich);
                }
            }
            histo->Fill(tot_energy_cut / *genParticle_energy);
        }
        
        // PCalo hit over simhit total
        for (auto const& [cut, histo] : hcalo_simhit_histos){
            float tot_energy_cut = 0.;
            for (int ich=0; ich< caloHit_energy->size(); ich++){
                //cout << "caloHit energy: "<< caloHit_energy->at(ich) << " | time: "<< caloHit_time->at(ich) <<endl;
                if (caloHit_energy->at(ich) > cut.first && caloHit_time->at(ich) < cut.second){
                    tot_energy_cut+= caloHit_energy->at(ich);
                }
            }
            histo->Fill(tot_energy_cut / tot_simhit_energy);
        }
        
    }
    
    outputfile->Write();
    outputfile->Close();
       

    return 0;

}