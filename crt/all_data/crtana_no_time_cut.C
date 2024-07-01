//Histograms with no time cut
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
namespace fs = std::filesystem;

int crtana_time() {
	std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

	Iterate over files in the directory
	for (const auto& entry : fs::directory_iterator(directory)) {
    		std::cout << entry.path().filename().string() << std::endl;
	}

	std::vector<std::string> filenames = {
		"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13178_crtana.root",
		"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13476_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13666_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13680_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13758_crtana.root",
		"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13690_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_crtana_flat.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13828_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13268_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13830_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_13689_13690_13693_13758_crtana_flat_etrig.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13689_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13470_crtana.root",
    		"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_fixed_channel_map_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13281_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13268_crtana_all.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13470_run13666_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13320_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13693_crtana.root",
    		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_13689_13690_13693_13758_crtana_flat.root"
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
        	return 1;
    	}

	if (!chain.GetListOfBranches()->FindObject("cl_sp_ts1")) {
        	std::cerr << "Branch 'cl_sp_ts1' not found in the TChain. Please check the variable name." << std::endl;
        	return 1;
    	}

	TH1D* histogram1_f = new TH1D("histogram1Df_x", "Front face (x)", 50, -400, 400);
	TH1D* histogram2_f = new TH1D("histogram1Df_y", "Front face (y)", 50, -400, 400);
	TH1D* histogram3_f = new TH2D("histogram2Df", "Front face", 50, -400, 400, 50, -400, 400);
	TH1D* histogram1_b = new TH1D("histogram1Db_x", "Back face (x)", 50, -400, 400);
	TH1D* histogram2_b = new TH1D("histogram1Db_y", "Back face (y)", 50, -400, 400);
        TH1D* histogram3_b = new TH2D("histogram2Db", "back face", 50, -400, 400, 50, -400, 400);
	TH1D* histogram1_l = new TH1D("histogram1Dl_y", "Left face (y)", 50, -400, 400);
        TH1D* histogram2_l = new TH1D("histogram1Dl_z", "Left face (z)", 50, -200, 800);
        TH1D* histogram3_l = new TH2D("histogram2Dl", "Left face", 50, -400, 400, 50, -200, 800);
        TH1D* histogram1_r = new TH1D("histogram1Dr_y", "Right face (y)", 50, -400, 400);
        TH1D* histogram2_r = new TH1D("histogram1Dr_z", "Right face (z)", 50, -200, 800);
        TH1D* histogram3_r = new TH2D("histogram2Dr", "Right face", 50, -400, 400, 50, -200, 800);
        TH1D* histogram1_t = new TH1D("histogram1Dt_x", "Top face (x)", 50, -400, 400);
        TH1D* histogram2_t = new TH1D("histogram1Dt_z", "Top face (z)", 50, -200, 800);
        TH1D* histogram3_t = new TH2D("histogram2Dt", "Top face", 50, -400, 400, 50, -200, 800);
        TH1D* histogram1_d = new TH1D("histogram1Dd_x", "Bottom face (x)", 50, -400, 400);
        TH1D* histogram2_d = new TH1D("histogram1Dd_z", "Bottom face (z)", 50, -200, 800);
        TH1D* histogram3_d = new TH2D("histogram2Dd", "Bottom face", 50, -400, 400, 50, -200, 800);


    	std::vector<bool>* cl_has_sp = nullptr; 
    	std::vector<double>* cl_sp_x = nullptr;  
    	std::vector<double>* cl_sp_y = nullptr;
	std::vector<double>* cl_sp_z = nullptr;
	std::vector<double>* cl_sp_ts1 = nullptr;

    	chain.SetBranchAddress("cl_has_sp", &cl_has_sp);
    	chain.SetBranchAddress("cl_sp_x", &cl_sp_x);
    	chain.SetBranchAddress("cl_sp_y", &cl_sp_y);
	chain.SetBranchAddress("cl_sp_z", &cl_sp_z);
	chain.SetBranchAddress("cl_sp_ts1", &cl_sp_ts1);
	
    	for (Long64_t i = 0; i < n_entries; ++i) {
        	chain.GetEntry(i);
        	if ((i % 1000) == 0) {
            		std::cout << i << "k" << std::endl;
        	}
        	if (cl_has_sp == nullptr || !cl_has_sp->at(0)) {
            		continue;
        	}
        	for (size_t j = 0; j < cl_sp_x->size(); ++j) {
            		double t1 = cl_sp_ts1->at(j);
			double x = cl_sp_x->at(j);
            		double y = cl_sp_y->at(j);
			double z = cl_sp_z->at(j);
			if (-450 < y < -350) {				
			}
			else if (350 < y < 450) {
			}
			else if (

        	}
   	}







