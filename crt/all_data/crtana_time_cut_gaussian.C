#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <TChain.h>
#include <TF1.h>
#include <TF2.h>


namespace fs = std::filesystem;

int crtana_time_cut_gaussian() {
    std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

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

	TH1D* histogram1_f_t_f = new TH1D("histogram1Dftf_x", "Front face (x)", 10, -400, 400);
	TH1D* histogram2_f_t_f = new TH1D("histogram1Dftf_y", "Front face (y)", 10, -300, 400);
	TH2D* histogram3_f_t_f = new TH2D("histogram2Dftf", "Front face", 10, -400, 400, 10, -350, 400);
	TH1D* histogram1_b_t_f = new TH1D("histogram1Dbtf_x", "Back face (x)", 10, -400, 400);
	TH1D* histogram2_b_t_f = new TH1D("histogram1Dbtf_y", "Back face (y)", 10, -350, 400);
	TH2D* histogram3_b_t_f = new TH2D("histogram2Dbtf", "Back face", 10, -400, 400, 10, -350, 400);

	TCanvas* c1_f_t_f = new TCanvas("c1_f_t_f", "Front Face (x)", 800, 600);
	TCanvas* c2_f_t_f = new TCanvas("c2_f_t_f", "Front Face (y)", 800, 600);
	TCanvas* c3_f_t_f = new TCanvas("c3_f_t_f", "Front Face", 800, 600);
	TCanvas* c1_b_t_f = new TCanvas("c1_b_t_f", "Back Face (x)", 800, 600);
	TCanvas* c2_b_t_f = new TCanvas("c2_b_t_f", "Back Face (y)", 800, 600);
	TCanvas* c3_b_t_f = new TCanvas("c3_b_t_f", "Back Face", 800, 600);

	for (Long64_t i = 0; i < n_entries; ++i) {
        	chain.GetEntry(i);

        	if ((i % 10000) == 0) {
            		std::cout << i << "k" << std::endl;
       		}

        	if (cl_has_sp == nullptr || cl_has_sp->size() == 0 || !cl_has_sp->at(0)) {
            		continue;
        	}

        	for (size_t j = 0; j < cl_sp_x->size(); ++j) {
            		double t1 = cl_sp_ts1->at(j);
            		double x = cl_sp_x->at(j);
            		double y = cl_sp_y->at(j);
            		double z = cl_sp_z->at(j);
	    		if (1529.8e3 < t1 && t1 < 1532.8e3) {		
            			if (y > -350) {  //cut off feet of detector
	    				if (-250 < z && z < -150) {
                				histogram1_f_t_f->Fill(x);
                				histogram2_f_t_f->Fill(y);
                				histogram3_f_t_f->Fill(x, y);
            					}	 
					else if (750 < z && z < 850) {
   			             	histogram1_b_t_f->Fill(x);
                			histogram2_b_t_f->Fill(y);
                			histogram3_b_t_f->Fill(x, y);
            				}	
				}
			}
       	 	}
    	}

    // Fit functions


	TF1* fit1D_x_f = new TF1("fit1D_x_f", "gaus", -400, 400);
	TF1* fit1D_y_f = new TF1("fit1D_y_f", "gaus", -350, 400);
	TF1* fit1D_x_b = new TF1("fit1D_x_b", "gaus", -400, 400);  
	TF1* fit1D_y_b = new TF1("fit1D_y_b", "gaus", -350, 400);
	TF2* fit2D_f = new TF2("fit2D_f", "[0]*exp(-0.5*((x-[1])*(x-[1])/([2]*[2]) + (y-[3])*(y-[3])/([4]*[4])))", -400, 400, -200, 800);
    	TF2* fit2D_b = new TF2("fit2D_b", "[0]*exp(-0.5*((x-[1])*(x-[1])/([2]*[2]) + (y-[3])*(y-[3])/([4]*[4])))", -400, 400, -200, 800);

	fit2D_f->SetParameters(1,0,100,0,100);
	fit2D_b->SetParameters(1,0,100,0,100);

   	histogram1_f_t_f->Fit("fit1D_x_f");
   	histogram2_f_t_f->Fit("fit1D_y_f");
    	histogram3_f_t_f->Fit("fit2D_f");
  	histogram1_b_t_f->Fit("fit1D_x_b");
    	histogram2_b_t_f->Fit("fit1D_y_b");
 	histogram3_b_t_f->Fit("fit2D_b");

    	std::cout << "Fit parameters for Front face (2D): " << std::endl;
    	fit2D_f->Print();

    	std::cout << "Fit parameters for Back face (2D): " << std::endl;
    	fit2D_b->Print();

    	c1_f_t_f->cd();
    	histogram1_f_t_f->Draw();
	fit1D_x_f->Draw("same");
    	c2_f_t_f->cd();
    	histogram2_f_t_f->Draw();
	fit1D_y_f->Draw("same");
    	c3_f_t_f->cd();
	c3_f_t_f->SetLogz();
    	histogram3_f_t_f->Draw("COLZ");
	fit2D_f->Draw("same CONT3");

    	c1_b_t_f->cd();
    	histogram1_b_t_f->Draw();
	fit1D_x_b->Draw("same");
    	c2_b_t_f->cd();
    	histogram2_b_t_f->Draw();
	fit1D_y_b->Draw("same");
    	c3_b_t_f->cd();
    	c3_b_t_f->SetLogz();
	histogram3_b_t_f->Draw("COLZ");
	fit2D_b->Draw("same CONT3");

	c1_f_t_f->SaveAs("Front_face_x_t_f_all.png");
	c2_f_t_f->SaveAs("Front_face_y_t_f_all.png");
	c3_f_t_f->SaveAs("Front_face_t_f_all.png");

	c1_b_t_f->SaveAs("Back_face_x_t_f_all.png");
	c2_b_t_f->SaveAs("Back_face_y_t_f_all.png");
	c3_b_t_f->SaveAs("Back_face_t_f_all.png");

	return 0;
	}

