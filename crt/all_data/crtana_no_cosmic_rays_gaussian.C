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


int crtana_no_cosmic_rays_gaussian() {
    std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

        std::vector<std::string> filenames = {       
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13178_crtana.root",
		//"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13476_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13666_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13680_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13758_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13690_crtana.root", 
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13828_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13268_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13830_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13689_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13470_crtana.root",
                "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13688_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_fixed_channel_map_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13281_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13268_crtana_all.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13470_run13666_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13320_crtana.root",
                //"/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13693_crtana.root"
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

	TH1D* histogram1_f_t = new TH1D("histogram1Dft_x", "Front face (x)", 10, -400, 400);
	TH1D* histogram2_f_t = new TH1D("histogram1Dft_y", "Front face (y)", 10, -400, 400);
	TH2D* histogram3_f_t = new TH2D("histogram2Dft", "Front face", 10, -400, 400, 10, -400, 400);
	TH1D* histogram1_b_t = new TH1D("histogram1Dbt_x", "Back face (x)", 10, -400, 400);
	TH1D* histogram2_b_t = new TH1D("histogram1Dbt_y", "Back face (y)", 10, -400, 400);
	TH2D* histogram3_b_t = new TH2D("histogram2Dbt", "Back face", 10, -400, 400, 10, -400, 400);

        TH1D* histogram1_f_n = new TH1D("histogram1Dfn_x", "Front face no beam (x)", 10, -400, 400);
        TH1D* histogram2_f_n = new TH1D("histogram1Dfn_y", "Front face no beam (y)", 10, -400, 400);
        TH2D* histogram3_f_n = new TH2D("histogram2Dfn", "Front face no beam", 10, -400, 400, 10, -400, 400);
        TH1D* histogram1_b_n = new TH1D("histogram1Dbn_x", "Back face no beam (x)", 10, -400, 400);
        TH1D* histogram2_b_n = new TH1D("histogram1Dbn_y", "Back face no beam (y)", 10, -400, 400);
        TH2D* histogram3_b_n = new TH2D("histogram2Dbn", "Back face no beam", 10, -400, 400, 10, -400, 400);

	TCanvas* c1_f_nocr = new TCanvas("c1_f_nocr", "Front Face (x)", 800, 600);
	TCanvas* c2_f_nocr = new TCanvas("c2_f_nocr", "Front Face (y)", 800, 600);
	TCanvas* c3_f_nocr = new TCanvas("c3_f_nocr", "Front Face", 800, 600);
	TCanvas* c1_b_nocr = new TCanvas("c1_b_nocr", "Back Face (x)", 800, 600);
	TCanvas* c2_b_nocr = new TCanvas("c2_b_nocr", "Back Face (y)", 800, 600);
	TCanvas* c3_b_nocr = new TCanvas("c3_b_nocr", "Back Face", 800, 600);
	TCanvas* c3_f_n = new TCanvas("c2_f_n", "Front Face cosmic rays", 800, 600);
        TCanvas* c3_b_n = new TCanvas("c3_f_n", "Back Face cosmic rays", 800, 600);

	double start_spill = 1529e3;
	double end_spill = 1533e3;
	double start_time = 0;
	double end_time = 1000e6;
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
	    		if (start_spill < t1 && t1 < end_spill) {		
            			if (y > -350) {  //cut off feet of detector
	    				if (-250 < z && z < -150) {
                				histogram1_f_t->Fill(x);
                				histogram2_f_t->Fill(y);
                				histogram3_f_t->Fill(x, y);
            				}	 
					else if (750 < z && z < 850) {
   			             		histogram1_b_t->Fill(x);
                				histogram2_b_t->Fill(y);
                				histogram3_b_t->Fill(x, y);
            				}	
				}
			}
			else if ((start_time < t1 && t1 < start_spill) || (end_spill < t1 && t1 < end_time)) {
				if (y > -350) {
					if (-250 < z && z < -150) {
                                                histogram1_f_n->Fill(x);
                                                histogram2_f_n->Fill(y);
                                                histogram3_f_n->Fill(x, y);		
					}
					else if (750 < z && z < 850) {
                                                histogram1_b_n->Fill(x);
                                                histogram2_b_n->Fill(y);
                                                histogram3_b_n->Fill(x, y);			
					}
				}
			}
       	 	}
    	}



//scale no beam areas
	double total_time = end_time - start_time;
	double spill_time = end_spill - start_spill;
	double no_spill = total_time - spill_time;
	double weight = spill_time / no_spill;	

	histogram1_f_n->Scale(weight);
	histogram2_f_n->Scale(weight);
	histogram3_f_n->Scale(weight);
	histogram1_b_n->Scale(weight);
	histogram2_b_n->Scale(weight);
	histogram3_b_n->Scale(weight);

    	TH2F *histogram1_f_neut = (TH2F*)histogram1_f_t->Clone("histogram1_fn");
    	histogram1_f_neut->Add(histogram1_f_n, -1);
        TH2F *histogram2_f_neut = (TH2F*)histogram2_f_t->Clone("histogram2_fn");
        histogram2_f_neut->Add(histogram2_f_n, -1);
        TH2F *histogram3_f_neut = (TH2F*)histogram3_f_t->Clone("histogram3_fn");
        histogram3_f_neut->Add(histogram3_f_n, -1);
        TH2F *histogram1_b_neut = (TH2F*)histogram1_b_t->Clone("histogram1_bn");
        histogram1_b_neut->Add(histogram1_b_n, -1);
        TH2F *histogram2_b_neut = (TH2F*)histogram2_b_t->Clone("histogram2_bn");
        histogram2_b_neut->Add(histogram2_b_n, -1);
        TH2F *histogram3_b_neut = (TH2F*)histogram3_b_t->Clone("histogram3_bn");
        histogram3_b_neut->Add(histogram3_b_n, -1);


        TF1* fit1D_x_f_neut = new TF1("fit1D_x_f_neut", "gaus", -400, 400);
        TF1* fit1D_y_f_neut = new TF1("fit1D_y_f_neut", "gaus", -350, 400);
        TF1* fit1D_x_b_neut = new TF1("fit1D_x_b_neut", "gaus", -400, 400);
        TF1* fit1D_y_b_neut = new TF1("fit1D_y_b_neut", "gaus", -350, 400);
        TF2* fit2D_f_neut = new TF2("fit2D_f_neut", "[0]*exp(-0.5*((x-[1])*(x-[1])/([2]*[2]) + (y-[3])*(y-[3])/([4]*[4])))", -400,400,-400,400);
        TF2* fit2D_b_neut = new TF2("fit2D_b_neut", "[0]*exp(-0.5*((x-[1])*(x-[1])/([2]*[2]) + (y-[3])*(y-[3])/([4]*[4])))", -400,400,-400,400);

        fit2D_f_neut->SetParameters(1, 0, 100, 0, 100);
        fit2D_b_neut->SetParameters(1, 0, 100, 0, 100);

        histogram1_f_neut->Fit("fit1D_x_f_neut");
        histogram2_f_neut->Fit("fit1D_y_f_neut");
        histogram3_f_neut->Fit("fit2D_f_neut");
        histogram1_b_neut->Fit("fit1D_x_b_neut");
        histogram2_b_neut->Fit("fit1D_y_b_neut");
        histogram3_b_neut->Fit("fit2D_b_neut");

        std::cout << "Fit parameters for Front face (2D): " << std::endl;
        fit2D_f_neut->Print();

        std::cout << "Fit parameters for Back face (2D): " << std::endl;
        fit2D_b_neut->Print();

        double chi2f = fit2D_f_neut->GetChisquare();
        int ndff = fit2D_f_neut->GetNDF();
        double chi2b = fit2D_b_neut->GetChisquare();
        int ndfb = fit2D_b_neut->GetNDF();
        double chi_per_deg_f = chi2f / ndff;
        double chi_per_deg_b = chi2b / ndfb;

        std::cout << "Front Chi-Squared: " << chi2f << std::endl;
        std::cout << "Front Number Degrees Freedom: " << ndff << std::endl;
        std::cout << "Front Chi-Squared/Degrees Freedom: " << chi_per_deg_f << std::endl;
        std::cout << "Back Chi-Squared: " << chi2b << std::endl;
        std::cout << "Back Number Degrees Freedom: " << ndfb << std::endl;
        std::cout << "Back Chi-Squared/Degrees Freedom: " << chi_per_deg_b << std::endl;


        //TPaveText *ptf = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
        //ptf->AddText(Form("Chi2: %.2f", chi2f));
        //ptf->AddText(Form("NDF: %d", ndff));
        //ptf->AddText(Form("Chi2/NDF: %.2f", chi2f/ndff));
        //ptf->SetFillColor(0);
        //ptf->SetTextColor(kRed);

        //TPaveText *ptb = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");
        //ptb->AddText(Form("Chi2: %.2f", chi2b));
        //ptb->AddText(Form("NDF: %d", ndfb));
        //ptb->AddText(Form("Chi2/NDF: %.2f", chi2b/ndfb));
        //ptb->SetFillColor(0);
        //ptb->SetTextColor(kRed);
	
        gStyle->SetOptStat(0);
        double entriesf = histogram3_f_neut->GetEntries();
        std::cout << "Front entries: " << entriesf << std::endl;
        double entriesb = histogram3_b_neut->GetEntries();
        std::cout << "Back entries: " << entriesb << std::endl;

	c1_f_nocr->cd();
	histogram1_f_neut->GetXaxis()->SetTitle("X (cm)");
	histogram1_f_neut->GetYaxis()->SetTitle("Number of Hits");
	histogram1_f_neut->Draw("HIST");
	fit1D_x_f_neut->Draw("same");

	c2_f_nocr->cd();
	histogram2_f_neut->GetXaxis()->SetTitle("Y (cm)");
	histogram2_f_neut->GetYaxis()->SetTitle("Number of Hits");
	histogram2_f_neut->Draw("HIST");
	fit1D_y_f_neut->Draw("same");

	c3_f_nocr->cd();
	//c3_f_nocr->SetLogz();
	histogram3_f_neut->GetXaxis()->SetTitle("X (cm)");
	histogram3_f_neut->GetYaxis()->SetTitle("Y (cm)");
	histogram3_f_neut->Draw("COLZ");
	fit2D_f_neut->Draw("CONT3 SAME");
	//ptf->Draw();

	c1_b_nocr->cd();
	histogram1_b_neut->GetXaxis()->SetTitle("X (cm)");
	histogram1_b_neut->GetYaxis()->SetTitle("Number of Hits");
	histogram1_b_neut->Draw("HIST");
	fit1D_x_b_neut->Draw("same");

	c2_b_nocr->cd();
	histogram2_b_neut->GetXaxis()->SetTitle("Y (cm)");
	histogram2_b_neut->GetYaxis()->SetTitle("Number of Hits");
	histogram2_b_neut->Draw("HIST");
	fit1D_y_b_neut->Draw("same");

	c3_b_nocr->cd();
	//c3_b_nocr->SetLogz();
	histogram3_b_neut->GetXaxis()->SetTitle("X (cm)");
	histogram3_b_neut->GetYaxis()->SetTitle("Y (cm)");
	histogram3_b_neut->Draw("COLZ");
	fit2D_b_neut->Draw("CONT3 SAME");
	//ptb->Draw();



	c3_f_n->cd();
        histogram3_f_n->GetXaxis()->SetTitle("X (cm)");
        histogram3_f_n->GetYaxis()->SetTitle("Y (cm)");
	histogram3_f_n->Draw("COLZ");
	c3_b_n->cd();
        histogram3_b_n->GetXaxis()->SetTitle("X (cm)");
        histogram3_b_n->GetYaxis()->SetTitle("Y (cm)");
	histogram3_b_n->Draw("COLZ");

	c3_f_n->SaveAs("Front_face_cosmic_rays_fitted");
	c3_b_n->SaveAs("Back_face_cosmic_rays_fitted");
	c1_f_nocr->SaveAs("Front_face_x_nocr_fitted.png");
	c2_f_nocr->SaveAs("Front_face_y_nocr_fitted.png");
	c3_f_nocr->SaveAs("Front_face_nocr_fitted.png");

	c1_b_nocr->SaveAs("Back_face_x_nocr_fitted.png");
	c2_b_nocr->SaveAs("Back_face_y_nocr_fitted.png");
	c3_b_nocr->SaveAs("Back_face_nocr_fitted.png");

	return 0;
	}

