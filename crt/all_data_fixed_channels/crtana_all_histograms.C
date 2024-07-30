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


int crtana_all_histograms() {
    	std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

	//import files
	std::vector<std::string> filenames = {
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_crtana_22jul2024.root",
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13689_crtana_22jul2024.root",
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13690_crtana_22jul2024.root",
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
     
	//Create canvases
	//Time
	TCanvas* c_t_t = new TCanvas("c_t_t", "t1", 800, 600);

	//No time cut
	TCanvas* c1_f = new TCanvas("c1_f", "Front Face (x)", 800, 600);
	TCanvas* c2_f = new TCanvas("c2_f", "Front Face (y)", 800, 600);
	TCanvas* c3_f = new TCanvas("c3_f", "Front Face", 800, 600);
	TCanvas* c1_b = new TCanvas("c1_b", "Back Face (x)", 800, 600);
	TCanvas* c2_b = new TCanvas("c2_b", "Back Face (y)", 800, 600);
	TCanvas* c3_b = new TCanvas("c3_b", "Back Face", 800, 600);
	TCanvas* c1_l = new TCanvas("c1_l", "Left Face (y)", 800, 600);
	TCanvas* c2_l = new TCanvas("c2_l", "Left Face (z)", 800, 600);
	TCanvas* c3_l = new TCanvas("c3_l", "Left Face", 800, 600);
	TCanvas* c1_r = new TCanvas("c1_r", "Right Face (y)", 800, 600);
	TCanvas* c2_r = new TCanvas("c2_r", "Right Face (z)", 800, 600);
	TCanvas* c3_r = new TCanvas("c3_r", "Right Face", 800, 600);
	TCanvas* c1_t = new TCanvas("c1_t", "Top Face (x)", 800, 600);
	TCanvas* c2_t = new TCanvas("c2_t", "Top Face (z)", 800, 600);
	TCanvas* c3_t = new TCanvas("c3_t", "Top Face", 800, 600);
	TCanvas* c1_d = new TCanvas("c1_d", "Bottom Face (x)", 800, 600);
	TCanvas* c2_d = new TCanvas("c2_d", "Bottom Face (z)", 800, 600);
	TCanvas* c3_d = new TCanvas("c3_d", "Bottom Face", 800, 600);
	
	//Time cut
        TCanvas* c1_f_t = new TCanvas("c1_f_t", "Front Face (x)", 800, 600);
        TCanvas* c2_f_t = new TCanvas("c2_f_t", "Front Face (y)", 800, 600);
        TCanvas* c3_f_t = new TCanvas("c3_f_t", "Front Face", 800, 600);	

	//Time cut gaussian
        TCanvas* c1_f_t_g = new TCanvas("c1_f_t_g", "Front Face (x)", 800, 600);
        TCanvas* c2_f_t_g = new TCanvas("c2_f_t_g", "Front Face (y)", 800, 600);
        TCanvas* c3_f_t_g = new TCanvas("c3_f_t_g", "Front Face", 800, 600);

	//No cosmic rays
        TCanvas* c1_f_nocr = new TCanvas("c1_f_nocr", "Front Face (x)", 800, 600);
        TCanvas* c2_f_nocr = new TCanvas("c2_f_nocr", "Front Face (y)", 800, 600);
        TCanvas* c3_f_nocr = new TCanvas("c3_f_nocr", "Front Face", 800, 600);
        TCanvas* c1_f_n = new TCanvas("c1_f_n", "Front Face cosmic rays (x)", 800, 600);
        TCanvas* c2_f_n = new TCanvas("c2_f_n", "Front Face cosmic rays (y)", 800, 600); 
        TCanvas* c3_f_n = new TCanvas("c3_f_n", "Front Face cosmic rays", 800, 600);

	//No cosmic rays gaussian
        TCanvas* c1_f_nocr_g = new TCanvas("c1_f_nocr_g", "Front Face (x)", 800, 600);
        TCanvas* c2_f_nocr_g = new TCanvas("c2_f_nocr_g", "Front Face (y)", 800, 600);
        TCanvas* c3_f_nocr_g = new TCanvas("c3_f_nocr_g", "Front Face", 800, 600);

	//No beam and time cut
	TCanvas* c2_f_t_n = new TCanvas("c2_f_t_n", "Front Face Beam and No Beam (y)", 800, 600); 

	//Define times 
	double start_spill = 1529e3;
        double end_spill = 1533e3;
        double start_time = 0;
        double end_time = 1000e6;

	//number of entries
	//double entries_number = 25000;


	//Outer for-loop (loop through entries)
        for (Long64_t i = 0; i < n_entries; ++i) {
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


	//Set contours
        const int nContours = 1;
        double sigma1_f = fit2D_f->GetParameter(0) * exp(-0.5);
        double contours_f[nContours] = {sigma1_f};
        fit2D_f->SetContour(nContours, contours_f);

        double sigma1_f_neut = fit2D_f_neut->GetParameter(0) * exp(-0.5);
        double contours_f_neut[nContours] = {sigma1_f_neut};
        fit2D_f_neut->SetContour(nContours, contours_f_neut);

	//purple gradient	
    	//const int numLevels = 100;
    	//int colors[numLevels];
	//for (int i = 0; i < numLevels; ++i) {
    	//	double ratio = (double)i / (numLevels - 1);
    	//	int r = (int)(230 + ratio * (48 - 230));
    	//	int g = (int)(204 + ratio * (25 - 204));
    	//	int b = (int)(255 + ratio * (52 - 255));
    	//	colors[i] = TColor::GetColor(r, g, b);
	//}
	//gStyle->SetPalette(numLevels, colors); 
	// greyscale
	//const int numLevels = 100;
	//int colors[numLevels];
	//for (int i = 0; i < numLevels; ++i) {
	//    	colors[i] = TColor::GetColor((255 - (i * 255 / (numLevels - 1))),
	//				 (255 - (i * 255 / (numLevels - 1))),
	//				 (255 - (i * 255 / (numLevels - 1))));
	//}
	//gStyle->SetPalette(numLevels, colors);

	//histogram3_f->SetContour(numLevels);
	//histogram3_b->SetContour(numLevels);
	//histogram3_r->SetContour(numLevels);
	//histogram3_l->SetContour(numLevels);
	//histogram3_t->SetContour(numLevels);
	//histogram3_d->SetContour(numLevels);
	//histogram3_f_t->SetContour(numLevels);
	//histogram3_f_t_g->SetContour(numLevels);
	//histogram3_f_n->SetContour(numLevels);
	//histogram3_f_neut->SetContour(numLevels);
	//histogram3_f_neut_g->SetContour(numLevels);


	//Set axes and draw
	//Time
	c_t_t->cd();
        histogram_t_t->GetXaxis()->SetTitle("Time (ns)");
        histogram_t_t->GetYaxis()->SetTitle("Number of hits");
        histogram_t_t->Draw();

	//No time cut
        c1_f->cd();
        histogram1_f->GetXaxis()->SetTitle("X (cm)");
        histogram1_f->GetYaxis()->SetTitle("Number of Hits");
        histogram1_f->Draw();
        c2_f->cd();
        histogram2_f->GetXaxis()->SetTitle("Y (cm)");
        histogram2_f->GetYaxis()->SetTitle("Number of Hits");
        histogram2_f->Draw();
        c3_f->cd();
        c3_f->SetLogz();
        histogram3_f->GetXaxis()->SetTitle("X (cm)");
        histogram3_f->GetYaxis()->SetTitle("Y (cm)");
        histogram3_f->Draw("COLZ");
        c1_b->cd();
        histogram1_b->GetXaxis()->SetTitle("X (cm)");
        histogram1_b->GetYaxis()->SetTitle("Number of Hits");
        histogram1_b->Draw();
        c2_b->cd();
        histogram2_b->GetXaxis()->SetTitle("Y (cm)");
        histogram2_b->GetYaxis()->SetTitle("Number of Hits");
        histogram2_b->Draw();	
        c3_b->cd();
        c3_b->SetLogz();
        histogram3_b->GetXaxis()->SetTitle("X (cm)");
        histogram3_b->GetYaxis()->SetTitle("Y (cm)");
        histogram3_b->Draw("COLZ");
        c1_l->cd();
        histogram1_l->GetXaxis()->SetTitle("Y (cm)");
        histogram1_l->GetYaxis()->SetTitle("Number of Hits");
        histogram1_l->Draw();
        c2_l->cd();
        histogram2_l->GetXaxis()->SetTitle("Z (cm)");
        histogram2_l->GetYaxis()->SetTitle("Number of Hits");
        histogram2_l->Draw();
        c3_l->cd();
        c3_l->SetLogz();
        histogram3_l->GetXaxis()->SetTitle("Y (cm)");
        histogram3_l->GetYaxis()->SetTitle("Z (cm)");
        histogram3_l->Draw("COLZ");
        c1_r->cd();
        histogram1_r->GetXaxis()->SetTitle("Y (cm)");
        histogram1_r->GetYaxis()->SetTitle("Number of Hits");
        histogram1_r->Draw();
        c2_r->cd();
        histogram2_r->GetXaxis()->SetTitle("Z (cm)");
        histogram2_r->GetYaxis()->SetTitle("Number of Hits");
        histogram2_r->Draw();
        c3_r->cd();
        c3_r->SetLogz();
        histogram3_r->GetXaxis()->SetTitle("Y (cm)");
        histogram3_r->GetYaxis()->SetTitle("Z (cm)");
        histogram3_r->Draw("COLZ");
        c1_t->cd();
        histogram1_t->GetXaxis()->SetTitle("X (cm)");
        histogram1_t->GetYaxis()->SetTitle("Number of Hits");
        histogram1_t->Draw();
        c2_t->cd();
        histogram2_t->GetXaxis()->SetTitle("Z (cm)");
        histogram2_t->GetYaxis()->SetTitle("Number of Hits");
        histogram2_t->Draw();
        c3_t->cd();
        c3_t->SetLogz();
        histogram3_t->GetXaxis()->SetTitle("X (cm)");
        histogram3_t->GetYaxis()->SetTitle("Z (cm)");
        histogram3_t->Draw("COLZ");
        c1_d->cd();
        histogram1_d->GetXaxis()->SetTitle("X (cm)");
        histogram1_d->GetYaxis()->SetTitle("Number of Hits");
        histogram1_d->Draw();
        c2_d->cd();
        histogram2_d->GetXaxis()->SetTitle("Z (cm)");
        histogram2_d->GetYaxis()->SetTitle("Number of Hits");
        histogram2_d->Draw();
        c3_d->cd();
        c3_d->SetLogz();
        histogram3_d->GetXaxis()->SetTitle("X (cm)");
        histogram3_d->GetYaxis()->SetTitle("Z (cm)");
        histogram3_d->Draw("COLZ");

	//Time cut
        c1_f_t->cd();
        histogram1_f_t->GetXaxis()->SetTitle("X (cm)");
        histogram1_f_t->GetYaxis()->SetTitle("Number of Hits");
        histogram1_f_t->Draw();
        c2_f_t->cd();
        histogram2_f_t->GetXaxis()->SetTitle("Y (cm)");
        histogram2_f_t->GetYaxis()->SetTitle("Number of Hits");
        histogram2_f_t->Draw();
        c3_f_t->cd();
        histogram3_f_t->GetXaxis()->SetTitle("X (cm)");
        histogram3_f_t->GetYaxis()->SetTitle("Y (cm)");
        //c3_f_t->SetLogz();
        histogram3_f_t->Draw("COLZ");

	//Time cut gaussian
        c1_f_t_g->cd();
        histogram1_f_t_g->GetXaxis()->SetTitle("X (cm)");
        histogram1_f_t_g->GetYaxis()->SetTitle("Number of Hits");
        histogram1_f_t_g->Draw();
        fit1D_x_f->Draw("same");
        c2_f_t_g->cd();
        histogram2_f_t_g->GetXaxis()->SetTitle("Y (cm)");
        histogram2_f_t_g->GetYaxis()->SetTitle("Number of Hits");
        histogram2_f_t_g->Draw();
        fit1D_y_f->Draw("same");
        c3_f_t_g->cd();
        histogram3_f_t_g->GetXaxis()->SetTitle("X (cm)");
        histogram3_f_t_g->GetYaxis()->SetTitle("Y (cm)");
        histogram3_f_t_g->Draw("COLZ");
        fit2D_f->Draw("same CONT3"); 

	//Cosmic rays only
	double inverse_weight = 1 / weight;
        histogram3_f_n->Scale(inverse_weight);
        c1_f_n->cd();
        histogram1_f_n->GetXaxis()->SetTitle("X (cm)");
        histogram1_f_n->GetYaxis()->SetTitle("Number of Hits");
        histogram1_f_n->Draw("HIST");
        c2_f_n->cd();
        histogram2_f_n->GetXaxis()->SetTitle("Y (cm)");
        histogram2_f_n->GetYaxis()->SetTitle("Number of Hits");
        histogram2_f_n->Draw("HIST");
        c3_f_n->cd();
        histogram3_f_n->GetXaxis()->SetTitle("X (cm)");
        histogram3_f_n->GetYaxis()->SetTitle("Y (cm)");
        histogram3_f_n->Draw("COLZ");

	//No cosmic rays
        c1_f_nocr->cd();
        histogram1_f_neut->GetXaxis()->SetTitle("X (cm)");
        histogram1_f_neut->GetYaxis()->SetTitle("Number of Hits");
        histogram1_f_neut->Draw("HIST");
        c2_f_nocr->cd();
        histogram2_f_neut->GetXaxis()->SetTitle("Y (cm)");
        histogram2_f_neut->GetYaxis()->SetTitle("Number of Hits");
        histogram2_f_neut->Draw("HIST");
        c3_f_nocr->cd();
        //c3_f_nocr->SetLogz();
        histogram3_f_neut->GetXaxis()->SetTitle("X (cm)");
        histogram3_f_neut->GetYaxis()->SetTitle("Y (cm)");
        histogram3_f_neut->Draw("COLZ");
	//histogram3_f_neut->Draw("LEGO2 same");

	//No cosmic rays gaussian
        c1_f_nocr_g->cd();
        histogram1_f_neut_g->GetXaxis()->SetTitle("X (cm)");
        histogram1_f_neut_g->GetYaxis()->SetTitle("Number of Hits");
        histogram1_f_neut_g->Draw("HIST");
        fit1D_x_f_neut->Draw("same");
        c2_f_nocr_g->cd();
        histogram2_f_neut_g->GetXaxis()->SetTitle("Y (cm)");
        histogram2_f_neut_g->GetYaxis()->SetTitle("Number of Hits");
        histogram2_f_neut_g->Draw("HIST");
        fit1D_y_f_neut->Draw("same");
        c3_f_nocr_g->cd();
        histogram3_f_neut_g->GetXaxis()->SetTitle("X (cm)");
        histogram3_f_neut_g->GetYaxis()->SetTitle("Y (cm)");
        histogram3_f_neut_g->Draw("COLZ");
	fit2D_f_neut->Draw("CONT3 SAME");

	//Time cut and cosmic ray
	int purple = TColor::GetColor(128, 0, 128);
	double scale = histogram2_f_t->GetMaximum() / histogram2_f_n->GetMaximum();
	histogram2_f_n->Scale(scale);
	c2_f_t_n->cd();
	histogram2_f_t->SetLineColor(purple);
	histogram2_f_t->Draw("HIST");
	histogram2_f_n->SetLineColor(purple);
	histogram2_f_n->SetLineStyle(2);
	histogram2_f_n->Draw("HIST SAME");
	auto legend = new TLegend(0.1, 0.7, 0.3, 0.9);
	legend->AddEntry(histogram2_f_t, "Time Cut Data", "l");
	legend->AddEntry(histogram2_f_n, "No Beam Cosmic Ray Data", "l");
	legend->Draw();
	double rightAxisMin = histogram2_f_n->GetMinimum() / scale;
	double rightAxisMax = histogram2_f_n->GetMaximum() / scale;
	TGaxis *axis = new TGaxis(gPad->GetUxmax(), histogram2_f_t->GetMinimum(), 
				  gPad->GetUxmax(), histogram2_f_t->GetMaximum(), 
				  rightAxisMin, rightAxisMax, 510, "+L");
	axis->SetLabelSize(0.035);
	axis->SetTitleSize(0.035);
	axis->SetLabelFont(histogram2_f_t->GetYaxis()->GetLabelFont());
	axis->SetTitleFont(histogram2_f_t->GetYaxis()->GetTitleFont());
	axis->SetTitle("Number of Hits (Cosmic Rays)");	
 	axis->Draw();
	c2_f_t_n->Update();

	//Save as pngs
	//Time
        c_t_t->SaveAs("t1_fc_test.png");

	//No time cut
        c1_f->SaveAs("Front_face_x_no_time_cut_fc_test.png");
        c2_f->SaveAs("Front_face_y_no_time_cut_fc_test.png");
        c3_f->SaveAs("Front_face_no_time_cut_fc_test.png");
        c1_b->SaveAs("Back_face_x_no_time_cut_fc_test.png");
        c2_b->SaveAs("Back_face_y_no_time_cut_fc_test.png");
        c3_b->SaveAs("Back_face_no_time_cut_fc_test.png");
        c1_l->SaveAs("Left_face_y_no_time_cut_fc_test.png");
        c2_l->SaveAs("Left_face_z_no_time_cut_fc_test.png");
        c3_l->SaveAs("Left_face_no_time_cut_fc_test.png");
        c1_r->SaveAs("Right_face_y_no_time_cut_fc_test.png");
        c2_r->SaveAs("Right_face_z_no_time_cut_fc_test.png");
        c3_r->SaveAs("Right_face_no_time_cut_fc_test.png");
        c1_t->SaveAs("Top_face_x_no_time_cut_fc_test.png");
        c2_t->SaveAs("Top_face_z_no_time_cut_fc_test.png");
        c3_t->SaveAs("Top_face_no_time_cut_fc_test.png");
        c1_d->SaveAs("Bottom_face_x_no_time_cut_fc_test.png");
        c2_d->SaveAs("Bottom_face_z_no_time_cut_fc_test.png");
        c3_d->SaveAs("Bottom_face_no_time_cut_fc_test.png");	

	//Time cut
        c1_f_t->SaveAs("Front_face_x_time_cut_fc_test.png");
        c2_f_t->SaveAs("Front_face_y_time_cut_fc_test.png");
        c3_f_t->SaveAs("Front_face_time_cut_fc_test.png");

        //Time cut gaussian
        c1_f_t_g->SaveAs("Front_face_x_time_cut_gaussian_fc_test.png");
        c2_f_t_g->SaveAs("Front_face_y_time_cut_gaussian_fc_test.png");
        c3_f_t_g->SaveAs("Front_face_time_cut_gaussian_fc_test.png");	
	
	//Cosmic rays only
        c1_f_n->SaveAs("Front_face_x_cosmic_rays_fc_test.png");
        c2_f_n->SaveAs("Front_face_y_cosmic_rays_fc_test.png");
        c3_f_n->SaveAs("Front_face_cosmic_rays_fc_test.png");
       
	//No cosmic rays
        c1_f_nocr->SaveAs("Front_face_x_nocr_fc_test.png");
        c2_f_nocr->SaveAs("Front_face_y_nocr_fc_test.png");
        c3_f_nocr->SaveAs("Front_face_nocr_fc_all_purple.png");	

	//No cosmic rays gaussian
        c1_f_nocr_g->SaveAs("Front_face_x_nocr_gaussian_fc_test.png");
        c2_f_nocr_g->SaveAs("Front_face_y_nocr_gaussian_fc_test.png");
        c3_f_nocr_g->SaveAs("Front_face_nocr_gaussian_fc_test.png");

	//Time cut and cosmic ray distribution
	c2_f_t_n->SaveAs("Front_face_y_time_cut_cosmic_rays_fc_purple_scaled_test.png");

	return 0;
}


