#include "TCanvas.h"
#include "TF2.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <TChain.h>
#include <TList.h>
#include <TBranch.h>
#include <TFile.h>
#include <TSystem.h>
#include "TH1D.h"
#include <TStyle.h>
#include <TMath.h>

namespace fs = std::filesystem;

int crtana_time() {
	std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";
    
        // Iterate over files in the directory
        for (const auto& entry : fs::directory_iterator(directory)) {
        	std::cout << entry.path().filename().string() << std::endl;
        }
    
	
        std::vector<std::string> filenames = {
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13178_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13476_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13666_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13680_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13758_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13690_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13828_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13268_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13830_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13689_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13470_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_fixed_channel_map_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13281_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13268_crtana_all.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13470_run13666_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13320_crtana.root",
        "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13693_crtana.root"                                                       
        };
 	TChain chain("crtana/tree");

	for (const auto& filename : filenames) {
        	chain.Add(filename.c_str());
    	}	

    // Check if the TChain has entries
        Long64_t n_entries = chain.GetEntries();
        std::cout << "Number of entries in the TChain: " << n_entries << std::endl;
    
        if (n_entries == 0) {
        std::cerr << "No entries found in the TChain. Please check the filenames and tree names." << std::endl;
        }
    
	    if (!chain.GetListOfBranches()->FindObject("cl_sp_ts1")) {
        std::cerr << "Branch 'cl_sp_ts1' not found in the TChain. Please check the variable name." << std::endl;
        return 1;
    }
	TH1D* histogram_t_t = new TH1D("histogram2D", "Time", 300, 1520e3, 1540e3);
	
        std::vector<bool>* cl_has_sp = nullptr;
    	std::vector<double>* cl_sp_x = nullptr;
    	std::vector<double>* cl_sp_ts1 = nullptr;

    	chain.SetBranchAddress("cl_has_sp", &cl_has_sp);
    	chain.SetBranchAddress("cl_sp_x", &cl_sp_x);
    	chain.SetBranchAddress("cl_sp_ts1", &cl_sp_ts1);		

	for (int i = 0; i < n_entries; ++i) {
        	chain.GetEntry(i);
        	if ((i % 10000) == 0) {
            		std::cout << i << "k" << std::endl;
        	}
        	if (!cl_has_sp) {
            		continue;
        	}
        	for (size_t j = 0; j < cl_sp_x->size(); ++j) {
            		double t1 = cl_sp_ts1->at(j);
            		histogram_t_t->Fill(t1);
        	}
    	}

        gStyle->SetOptStat(0);
        double entries = histogram_t_t->GetEntries();
        std::cout << "Entries: " << entries << std::endl;

    	histogram_t_t->GetXaxis()->SetTitle("Time (ns)");
    	histogram_t_t->GetYaxis()->SetTitle("Number of hits");	
    	TCanvas* c_t_t = new TCanvas("c_t_t", "t1", 800, 600);
    	c_t_t->SetWindowPosition(0, 0);	
	histogram_t_t->Draw();
	c_t_t->Update();
    	c_t_t->SaveAs("t1.png");	
	//c_t_t->SaveAs("t1_all.png");	

	return 0;	
}

