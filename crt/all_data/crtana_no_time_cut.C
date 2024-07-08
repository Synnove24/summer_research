#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <TChain.h>

namespace fs = std::filesystem;

int crtana_no_time_cut() {
    std::string directory = "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/";

        std::vector<std::string> filenames = {                                                                                          "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13178_crtana.root",                                           "/pnfs/sbnd/persistent/users/hlay/crt_comm_summer_2024/run13466_crtana.root",
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

    TH1D* histogram1_f = new TH1D("histogram1Df_x", "Front face (x)", 50, -400, 400);
    TH1D* histogram2_f = new TH1D("histogram1Df_y", "Front face (y)", 50, -400, 400);
    TH2D* histogram3_f = new TH2D("histogram2Df", "Front face", 50, -400, 400, 50, -400, 400);
    TH1D* histogram1_b = new TH1D("histogram1Db_x", "Back face (x)", 50, -400, 400);
    TH1D* histogram2_b = new TH1D("histogram1Db_y", "Back face (y)", 50, -400, 400);
    TH2D* histogram3_b = new TH2D("histogram2Db", "Back face", 50, -400, 400, 50, -400, 400);
    TH1D* histogram1_l = new TH1D("histogram1Dl_y", "Left face (y)", 50, -400, 400);
    TH1D* histogram2_l = new TH1D("histogram1Dl_z", "Left face (z)", 50, -200, 800);
    TH2D* histogram3_l = new TH2D("histogram2Dl", "Left face", 50, -400, 400, 50, -200, 800);
    TH1D* histogram1_r = new TH1D("histogram1Dr_y", "Right face (y)", 50, -400, 400);
    TH1D* histogram2_r = new TH1D("histogram1Dr_z", "Right face (z)", 50, -200, 800);
    TH2D* histogram3_r = new TH2D("histogram2Dr", "Right face", 50, -400, 400, 50, -200, 800);
    TH1D* histogram1_t = new TH1D("histogram1Dt_x", "Top face (x)", 50, -400, 400);
    TH1D* histogram2_t = new TH1D("histogram1Dt_z", "Top face (z)", 50, -200, 800);
    TH2D* histogram3_t = new TH2D("histogram2Dt", "Top face", 50, -400, 400, 50, -200, 800);
    TH1D* histogram1_d = new TH1D("histogram1Dd_x", "Bottom face (x)", 50, -400, 400);
    TH1D* histogram2_d = new TH1D("histogram1Dd_z", "Bottom face (z)", 50, -200, 800);
    TH2D* histogram3_d = new TH2D("histogram2Dd", "Bottom face", 50, -400, 400, 50, -200, 800);

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

    for (Long64_t i = 0; i < n_entries; ++i) {
        chain.GetEntry(i);

        if ((i % 1000) == 0) {
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

            if (-450 < y && y < -350) {
                histogram1_d->Fill(x);
                histogram2_d->Fill(z);
                histogram3_d->Fill(x, z);
            } else if (350 < y && y < 450) {
                histogram1_t->Fill(x);
                histogram2_t->Fill(z);
                histogram3_t->Fill(x, z);
            } else if (-450 < x && x < -350) {
                histogram1_r->Fill(y);
                histogram2_r->Fill(z);
                histogram3_r->Fill(y, z);
            } else if (350 < x && x < 450) {
                histogram1_l->Fill(y);
                histogram2_l->Fill(z);
                histogram3_l->Fill(y, z);
            } else if (-250 < z && z < -150) {
                histogram1_f->Fill(x);
                histogram2_f->Fill(y);
                histogram3_f->Fill(x, y);
            } else if (750 < z && z < 850) {
                histogram1_b->Fill(x);
                histogram2_b->Fill(y);
                histogram3_b->Fill(x, y);
            }
        }
    }

    	c1_f->cd();
    	histogram1_f->Draw();
    	c2_f->cd();
    	histogram2_f->Draw();
    	c3_f->cd();
	c3_f->SetLogz();
    	histogram3_f->Draw("COLZ");

    	c1_b->cd();
    	histogram1_b->Draw();
    	c2_b->cd();
    	histogram2_b->Draw();
    	c3_b->cd();
    	c3_b->SetLogz();
	histogram3_b->Draw("COLZ");

    	c1_l->cd();
   	histogram1_l->Draw();
    	c2_l->cd();
   	histogram2_l->Draw();
    	c3_l->cd();
    	c3_l->SetLogz();
	histogram3_l->Draw("COLZ");

	c1_r->cd();
	histogram1_r->Draw();
	c2_r->cd();
	histogram2_r->Draw();
	c3_r->cd();
	c3_r->SetLogz();
	histogram3_r->Draw("COLZ");

	c1_t->cd();
	histogram1_t->Draw();
	c2_t->cd();
	histogram2_t->Draw();
	c3_t->cd();
	c3_t->SetLogz();
	histogram3_t->Draw("COLZ");

	c1_d->cd();
	histogram1_d->Draw();
	c2_d->cd();
	histogram2_d->Draw();
	c3_d->cd();
	c3_d->SetLogz();
	histogram3_d->Draw("COLZ");


	c1_f->SaveAs("Front_face_x.png");
	c2_f->SaveAs("Front_face_y.png");
	c3_f->SaveAs("Front_face.png");

	c1_b->SaveAs("Back_face_x.png");
	c2_b->SaveAs("Back_face_y.png");
	c3_b->SaveAs("Back_face.png");

	c1_l->SaveAs("Left_face_y.png");
	c2_l->SaveAs("Left_face_z.png");
	c3_l->SaveAs("Left_face.png");

	c1_r->SaveAs("Right_face_y.png");
	c2_r->SaveAs("Right_face_z.png");
	c3_r->SaveAs("Right_face.png");

	c1_t->SaveAs("Top_face_x.png");
	c2_t->SaveAs("Top_face_y.png");
	c3_t->SaveAs("Top_face.png");

	c1_d->SaveAs("Bottom_face_x.png");
	c2_d->SaveAs("Bottom_face_y.png");
	c3_d->SaveAs("Bottom_face.png");
	
	return 0;
	}

