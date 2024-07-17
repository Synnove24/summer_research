//import everything
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <iostream>
#include "TF1.h"
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "TPaveText.h"
#include "TF2.h"
#include <TStyle.h>
#include <TMath.h>

namespace fs = std::filesystem;


int crtana_all_hitstograms() {
    	std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

	//import files
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

        Long64_t n_entries = chain.GetEntries();

        std::cout << "Number of entries in the TChain: " << n_entries << std::endl;

        if (n_entries == 0) {
        std::cerr << "No entries found in the TChain. Please check the filenames and tree names." << std::endl;
        return 1;
        }

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






