///////////////////////////////////////////////////////////////////////////////
///      plotter.c
///----------------------------------------------------------------------------
///
///     AUTHOR: Jacob Miner
///
///  Saves the histograms from the meer.c or meer.cc output as PNG files
///    unfortunately, it isn't intelligent about the plots 
///    (ie. it doesn't look for every histogram in the root file)
///
///////////////////////////////////////////////////////////////////////////////


#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <TH1.h>
#include <TH2.h>
#include <TString.h>


// In order to run this code as a "macro" in Root, the  main function must
// have the same name as the source file

int HepMCReader_test(TString DataFileName)
{
	gROOT->Reset();

	//
	// Initialize input filestream and output file pointers 
	//
	
	TFile *file = new TFile(DataFileName+".root"); //load the input file
	TCanvas *c = new TCanvas; // make a canvas


	//get the histograms, plot to canvas and save to file
	TH1F *h = (TH1F*)file->Get("Ttbar_0_evtnr");
	h->Draw();
	c->Print("Ttbar_0_evtnr.png");
	c->Clear();  //clear the canvas for a new histogram
/*
	TH1F *h1 = (TH1F*)file->Get("L2_z_theta");
	h1->Draw();
	c->Print("L2_z_theta.png");
	c->Clear();

	TH1F *h2 = (TH1F*)file->Get("M1_M2");
	h2->Draw();
	c->Print("M1_M2.png");
	c->Clear();

	TH1F *h3 = (TH1F*)file->Get("V1_2j");
	h3->Draw();
	c->Print("V1_2j.png");
	c->Clear();

	TH1F *h4 = (TH1F*)file->Get("V1_Md1");
	h4->Draw();
	c->Print("V1_Md1.png");
	c->Clear();

	TH1F *h5 = (TH1F*)file->Get("V1_Md2");
	h5->Draw();
	c->Print("V1_Md2.png");
	c->Clear();

	TH1F *h6 = (TH1F*)file->Get("V1_Mj");
	h6->Draw();
	c->Print("V1_Mj.png");
	c->Clear();

	TH1F *h7 = (TH1F*)file->Get("V1_V2_H");
	h7->Draw();
	c->Print("V1_V2_H.png");
	c->Clear();

	TH1F *h8 = (TH1F*)file->Get("V1_V2_L");
	h8->Draw();
	c->Print("V1_V2_L.png");
	c->Clear();

	TH1F *h9 = (TH1F*)file->Get("V1_V2_all");
	h9->Draw();
	c->Print("V1_V2_all.png");
	c->Clear();

	TH1F *h10 = (TH1F*)file->Get("V1_pt");
	h10->Draw();
	c->Print("V1_pt.png");
	c->Clear();

	TH1F *h11 = (TH1F*)file->Get("V1inv");
	h11->Draw();
	c->Print("V1inv.png");
	c->Clear();

	TH1F *h12 = (TH1F*)file->Get("V1inv_full");
	h12->Draw();
	c->Print("V1inv_full.png");
	c->Clear();

	TH1F *h13 = (TH1F*)file->Get("pt12_m1_m2");
	h13->Draw();
	c->Print("pt12_m1_m2.png");
	c->Clear();

	TH1F *h14 = (TH1F*)file->Get("pt1_pt2");
	h14->Draw();
	c->Print("pt_1pt_2.png");
	c->Clear();
	
	TH1F *h15 = (TH1F*)file->Get("z_theta");
	h15->Draw();
	c->Print("z_theta.png");
	c->Clear();

	TH1F *h16 = (TH1F*)file->Get("z_theta_100");
	h16->Draw();
	c->Print("z_theta_100.png");
	c->Clear();

	TH1F *h17 = (TH1F*)file->Get("z_theta_300");
	h17->Draw();
	c->Print("z_theta_300.png");
	c->Clear();

	TH1F *h18 = (TH1F*)file->Get("z_theta_400");
	h18->Draw();
	c->Print("z_theta_400.png");
	c->Clear();
*/S
	return 0;
}
