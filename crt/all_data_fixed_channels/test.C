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
#include <TLegend.h>
#include <TGaxis.h>
#include <TColor.h>

namespace fs = std::filesystem;


int test() {
    	std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

	//import files
	std::vector<std::string> filenames = {
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_crtana_22jul2024.root",
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13689_crtana_22jul2024.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13690_crtana_22jul2024.root",
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13693_crtana_22jul2024.root",
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13758_crtana_22jul2024.root",

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

	//Define histograms
	//Time
        TH1D* histogram_t_t = new TH1D("histogram2D", "Time", 300, 1520e3, 1540e3);
	
	//No time cut
	TH1D* histogram1_f = new TH1D("histogram1_f", "Front face (x)", 128, -400, 400);
	TH1D* histogram2_f = new TH1D("histogram2_f", "Front face (y)", 128, -400, 400);
	TH2D* histogram3_f = new TH2D("histogram3_f", "Front face", 128, -400, 400, 128, -400, 400);
	TH1D* histogram1_b = new TH1D("histogram1_b", "Back face (x)", 128, -400, 400);
	TH1D* histogram2_b = new TH1D("histogram2_b", "Back face (y)", 128, -400, 400);
	TH2D* histogram3_b = new TH2D("histogram3_b", "Back face", 128, -400, 400, 128, -400, 400);
	TH1D* histogram1_l = new TH1D("histogram1_l", "Left face (y)", 128, -400, 400);
	TH1D* histogram2_l = new TH1D("histogram2_l", "Left face (z)", 128, -200, 800);
	TH2D* histogram3_l = new TH2D("histogram3_l", "Left face", 128, -400, 400, 128, -200, 800);
	TH1D* histogram1_r = new TH1D("histogram1_r", "Right face (y)", 128, -400, 400);
	TH1D* histogram2_r = new TH1D("histogram2_r", "Right face (z)", 128, -200, 800);
	TH2D* histogram3_r = new TH2D("histogram3_r", "Right face", 128, -400, 400, 128, -200, 800);
	TH1D* histogram1_t = new TH1D("histogram1_t", "Top face (x)", 128, -400, 400);
	TH1D* histogram2_t = new TH1D("histogram2_t", "Top face (z)", 128, -200, 800);
	TH2D* histogram3_t = new TH2D("histogram3_t", "Top face", 128, -400, 400, 128, -200, 800);
	TH1D* histogram1_d = new TH1D("histogram1_d", "Bottom face (x)", 128, -400, 400);
	TH1D* histogram2_d = new TH1D("histogram2_d", "Bottom face (z)", 128, -200, 800);
	TH2D* histogram3_d = new TH2D("histogram3_d", "Bottom face", 128, -400, 400, 128, -200, 800);	

	//Time cut
        TH1D* histogram1_f_t = new TH1D("histogram1_f_t", "Front face (x)", 10, -360, 360);
        TH1D* histogram2_f_t = new TH1D("histogram2_f_t", "Front face (y)", 10, -360, 360);
        TH2D* histogram3_f_t = new TH2D("histogram3_f_t", "Front face", 10, -360, 360, 10, -360, 360);
	
	//Time cut gaussian	       
        TH1D* histogram1_f_t_g = new TH1D("histogram1_f_t_g", "Front face (x)", 10, -360, 360);
        TH1D* histogram2_f_t_g = new TH1D("histogram2_f_t_g", "Front face (y)", 10, -360, 360);
        TH2D* histogram3_f_t_g = new TH2D("histogram3_f_t_g", "Front face", 10, -360, 360, 10, -360, 360);
       
	//No cosmic rays
        TH1D* histogram1_f_n = new TH1D("histogram1_f_n", "Front face no beam (x)", 10, -360, 360);
        TH1D* histogram2_f_n = new TH1D("histogram2_f_n", "Front face no beam (y)", 10, -360, 360);
        TH2D* histogram3_f_n = new TH2D("histogram3_f_n", "Front face no beam", 10, -360, 360, 10, -360, 360);
     
	//Define times 
	double start_spill = 1529e3;
        double end_spill = 1533e3;
        double start_time = 0;
        double end_time = 1000e6;

	//number of entries
	double entries_number = 175000;


	//Outer for-loop (loop through entries)
        for (Long64_t i = 0; i < entries_number; ++i) {
                chain.GetEntry(i);
                if ((i % 10000) == 0) {
                        std::cout << i << "k" << std::endl;
                }
                if (cl_has_sp == nullptr || cl_has_sp->size() == 0 || !cl_has_sp->at(0)) {
                        continue;
                }
        	for (size_t j = 0; j < cl_sp_x->size(); ++j) { //begin inner for-loop (loop through events)
            		double t1 = cl_sp_ts1->at(j);
            		double x = cl_sp_x->at(j);
            		double y = cl_sp_y->at(j);
            		double z = cl_sp_z->at(j);
	                histogram_t_t->Fill(t1);           	 	
			if (-450 < y && y < -350) {
				histogram1_d->Fill(x);
				histogram2_d->Fill(z);
				histogram3_d->Fill(x, z);
			} 
			if (350 < y && y < 450) {
				histogram1_t->Fill(x);
				histogram2_t->Fill(z);
				histogram3_t->Fill(x, z);
			} 
			if (-450 < x && x < -350) {
				histogram1_r->Fill(y);
				histogram2_r->Fill(z);
				histogram3_r->Fill(y, z);
			} 	
			if (350 < x && x < 450) {
				histogram1_l->Fill(y);
				histogram2_l->Fill(z);
				histogram3_l->Fill(y, z);
			} 	 	
			if (750 < z && z < 850) {
				histogram1_b->Fill(x);
				histogram2_b->Fill(y);
				histogram3_b->Fill(x, y);
			}
			if (-250 < z && z < -150) {
				if (y > -360 && y < 360) {
				if (x > -360 && x < 360) {
					histogram1_f->Fill(x);
					histogram2_f->Fill(y);
					histogram3_f->Fill(x, y);
					if (start_spill < t1 && t1 < end_spill) {
						histogram1_f_t->Fill(x);
						histogram2_f_t->Fill(y);
						histogram3_f_t->Fill(x, y);
						histogram1_f_t_g->Fill(x);
						histogram2_f_t_g->Fill(y);
						histogram3_f_t_g->Fill(x, y);
					}
					else if ((start_time < t1 && t1 < start_spill) || (end_spill < t1 && t1 < end_time)) {
						histogram1_f_n->Fill(x);
						histogram2_f_n->Fill(y);
						histogram3_f_n->Fill(x, y);	
					}
				}
				}
			}
		}
    	} 
	
	//Calculate scaling
        double total_time = end_time - start_time;
        double spill_time = end_spill - start_spill;
        double no_spill = total_time - spill_time;
        double weight = spill_time / no_spill;

	//Scale the no beam histograms
        histogram1_f_n->Scale(weight);
        histogram2_f_n->Scale(weight);
        histogram3_f_n->Scale(weight);
       
	//Subtract off cosmic ray distribution
        TH2F *histogram1_f_neut = (TH2F*)histogram1_f_t->Clone("histogram1_f_neut");
        histogram1_f_neut->Add(histogram1_f_n, -1);
        TH2F *histogram2_f_neut = (TH2F*)histogram2_f_t->Clone("histogram2_f_neut");
        histogram2_f_neut->Add(histogram2_f_n, -1);
        TH2F *histogram3_f_neut = (TH2F*)histogram3_f_t->Clone("histogram3_f_neut");
        histogram3_f_neut->Add(histogram3_f_n, -1);
        TH2F *histogram1_f_neut_g = (TH2F*)histogram1_f_t->Clone("histogram1_f_neut_g");
        histogram1_f_neut_g->Add(histogram1_f_n, -1);
        TH2F *histogram2_f_neut_g = (TH2F*)histogram2_f_t->Clone("histogram2_f_neut_g");
        histogram2_f_neut_g->Add(histogram2_f_n, -1);
        TH2F *histogram3_f_neut_g = (TH2F*)histogram3_f_t->Clone("histogram3_f_neut_g");
        histogram3_f_neut_g->Add(histogram3_f_n, -1);       

	//Fit functions
        TF1* fit1D_x_f = new TF1("fit1D_x_f", "gaus", -400, 400);
        TF1* fit1D_y_f = new TF1("fit1D_y_f", "gaus", -350, 400);
        TF2* fit2D_f = new TF2("fit2D_f", "[0]*exp(-0.5*((x-[1])*(x-[1])/([2]*[2]) + (y-[3])*(y-[3])/([4]*[4])))", -400, 400, -400, 400);
        TF1* fit1D_x_f_neut = new TF1("fit1D_x_f_neut", "gaus", -360, 360);
        TF1* fit1D_y_f_neut = new TF1("fit1D_y_f_neut", "gaus", -360, 360);
        TF2* fit2D_f_neut = new TF2("fit2D_f_neut", "[0]*exp(-0.5*((x-[1])*(x-[1])/([2]*[2]) + (y-[3])*(y-[3])/([4]*[4])))", -400,400,-400,400);
      
	//Set parameters
        fit2D_f->SetParameters(1,0,100,0,100);
        fit2D_f_neut->SetParameters(1, 0, 100, 0, 100); 

	//Fit
        histogram1_f_t_g->Fit("fit1D_x_f");
        histogram2_f_t_g->Fit("fit1D_y_f");
        histogram3_f_t_g->Fit("fit2D_f", "N");
        histogram1_f_neut->Fit("fit1D_x_f_neut");
        histogram2_f_neut->Fit("fit1D_y_f_neut");
        histogram3_f_neut->Fit("fit2D_f_neut", "N");
	
	//Print Fits
        std::cout << "Fit parameters for Front face (2D) Time cut data: " << std::endl;
        fit2D_f->Print();
        std::cout << "Fit parameters for Front face (2D) No cosmic rays: " << std::endl;
        fit2D_f_neut->Print();

	//Get Fit info
        double chi2f_t = fit2D_f->GetChisquare();
        int ndff_t = fit2D_f->GetNDF();
        double chi_per_deg_f_t = chi2f_t / ndff_t;        
	double mean_x2D_f_t = fit2D_f->GetParameter(1);
        double stddev_x2D_f_t = fit2D_f->GetParameter(2);
        double mean_y2D_f_t = fit2D_f->GetParameter(3);
        double stddev_y2D_f_t = fit2D_f->GetParameter(4);

        double chi2f_neut = fit2D_f_neut->GetChisquare();
        int ndff_neut = fit2D_f_neut->GetNDF();
        double chi_per_deg_f_neut = chi2f_neut / ndff_neut;
        double mean_x2D_f_neut = fit2D_f_neut->GetParameter(1);
        double stddev_x2D_f_neut = fit2D_f_neut->GetParameter(2);
        double mean_y2D_f_neut = fit2D_f_neut->GetParameter(3);
        double stddev_y2D_f_neut = fit2D_f_neut->GetParameter(4);

	//Manually print Fit info
        std::cout << "2D Front X Mean (time cut): " << mean_x2D_f_t << ", 2D Front X StdDev (time cut): " << stddev_x2D_f_t << std::endl;
        std::cout << "2D Front Y Mean (time cut): " << mean_y2D_f_t << ", 2D Front Y StdDev (time cut): " << stddev_y2D_f_t << std::endl;
        std::cout << "Front Chi-Squared (time cut): " << chi2f_t << std::endl;
        std::cout << "Front Number Degrees Freedom (time cut): " << ndff_t << std::endl;
        std::cout << "Front Chi-Squared/Degrees Freedom (time cut): " << chi_per_deg_f_t << std::endl;

        std::cout << "Front X Mean (no cosmic rays): " << mean_x2D_f_neut << ", Front X StdDev (no cosmic rays): " << stddev_x2D_f_neut << std::endl;
        std::cout << "Front Y Mean (no cosmic rays): " << mean_y2D_f_neut << ", Front Y StdDev (no cosmic rays): " << stddev_y2D_f_neut << std::endl;
        std::cout << "Front Chi-Squared (no cosmic rays): " << chi2f_neut << std::endl;
        std::cout << "Front Number Degrees Freedom (no cosmic rays): " << ndff_neut << std::endl;
        std::cout << "Front Chi-Squared/Degrees Freedom (no cosmic rays): " << chi_per_deg_f_neut << std::endl;
       
	//Remove box and print entries
        gStyle->SetOptStat(0);
        double entries = histogram_t_t->GetEntries();
        std::cout << "Entries (time): " << entries << std::endl;
        double entries_f = histogram3_f->GetEntries();
        std::cout << "Front entries: " << entries_f << std::endl;
        double entries_b = histogram3_b->GetEntries();
        std::cout << "Back entries: " << entries_b << std::endl;
        double entries_l = histogram3_l->GetEntries();
        std::cout << "Left entries: " << entries_l << std::endl;
        double entries_r = histogram3_r->GetEntries();
        std::cout << "Right entries: " << entries_r << std::endl;
        double entries_t = histogram3_t->GetEntries();
        std::cout << "Top entries: " << entries_t << std::endl;
        double entries_f_t = histogram3_f_t->GetEntries();
        std::cout << "Front entries (time cut): " << entries_f_t << std::endl;
	double entries_d = histogram3_d->GetEntries();
        std::cout << "Bottom entries: " << entries_d << std::endl;
        double entries_f_neut = histogram3_f_neut->GetEntries();
        std::cout << "Front entries (no cosmic rays): " << entries_f_neut << std::endl;


	return 0;
}


