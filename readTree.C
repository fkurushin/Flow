#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TStopwatch.h>

// Define max number of particles in the event
#define MAX_TRACKS 15000
#define HARMONICNUMBER 2
#define WEIGHT 1

Float_t Qy(Float_t phi)
{
	return WEIGHT * sin(HARMONICNUMBER * phi);
}

Float_t Qx(Float_t phi)
{
	return WEIGHT * cos(HARMONICNUMBER * phi);
}

Float_t PHI(Float_t momy, Float_t momx)
{
	return TMath::ATan2(momy, momx);
}

Float_t PSI_EP(Float_t Qx, Float_t Qy)
{
	return TMath::ATan2(Qy, Qx) / HARMONICNUMBER;
}

Float_t vnobs(std::vector<float> angles, float ep_angle)
{
	TH1F *htmp = new TH1F("h1", "h1 title", 100, -10, 10);
	Float_t tmp;
	Float_t mean = 0.0;

	for (Int_t i = 0; i < angles.size(); i++) {
		tmp = cos(HARMONICNUMBER * (angles[i] - ep_angle));
		htmp->Fill(tmp);
	}
	mean = htmp->GetMean();
	delete htmp;
	return mean;
}


void readTree(TString inListName, TString outFileName)
{
	// Start timer
	TStopwatch timer;
	timer.Start();

	// Read input from the file list
	std::ifstream file(inListName.Data());
	std::string line;
	TChain *chain = new TChain("mctree");
	Int_t nFiles = 0;
	while (std::getline(file, line))
	{
		chain->Add(line.c_str());
		nFiles++;
	}

		// Init variables to read to
		Float_t         bimp; // - impact parameter
  	Float_t         phi2;
  	Float_t         phi3;
  	Float_t         ecc2;
  	Float_t         ecc3;
  	Int_t           npart;
  	Int_t           nh; // - number of particles in the event
  	Float_t         momx[MAX_TRACKS];   //[nh] - px
  	Float_t         momy[MAX_TRACKS];   //[nh] - py
  	Float_t         momz[MAX_TRACKS];   //[nh] - pz
  	Float_t         ene[MAX_TRACKS];   //[nh] - Energy
  	Int_t           hid[MAX_TRACKS];   //[nh]
  	Int_t           pdg[MAX_TRACKS];   //[nh] - PDG code (211 - pion, 321 - kaon, 2212 - proton)
  	Short_t         charge[MAX_TRACKS];   //[nh]	- charge

		chain->SetBranchAddress("bimp", &bimp);
  	chain->SetBranchAddress("phi2", &phi2);
  	chain->SetBranchAddress("phi3", &phi3);
  	chain->SetBranchAddress("ecc2", &ecc2);
  	chain->SetBranchAddress("ecc3", &ecc3);
  	chain->SetBranchAddress("npart", &npart);
  	chain->SetBranchAddress("nh", &nh);
  	chain->SetBranchAddress("momx", momx);
  	chain->SetBranchAddress("momy", momy);
  	chain->SetBranchAddress("momz", momz);
 		chain->SetBranchAddress("ene", ene);
  	chain->SetBranchAddress("hid", hid);
  	chain->SetBranchAddress("pdg", pdg);
  	chain->SetBranchAddress("charge", charge);

	 // Manage output
	 TFile *fo = new TFile(outFileName.Data(), "recreate");

	/*Init my own variables*/
	float bimpValues[11] = {
       	 0.00, 4.17,
       	 6.03, 7.40,
       	 8.56, 9.63,
       	 10.59, 11.46,
       	 12.29, 13.30,
       	 14.94
   	 };


    	Float_t Qx_neg = 0.0, Qx_pos = 0.0;
    	Float_t Qy_neg = 0.0, Qy_pos = 0.0;
    	Float_t phi = 0.0;
    	Float_t qW = 0.0;
    	Float_t Psi_EP_pos = 0.0;
    	Float_t Psi_EP_neg = 0.0;

    	Float_t error = 0.0;
    	Float_t etta = 0.0;
    	Float_t momModule = 0.0;
    	Float_t momTransverse = 0.0;
    	Float_t pre_Resolution = 0.0;
			Float_t vnobs_pos = 0.0;
			Float_t vnobs_neg = 0.0;

    	Bool_t flag = false;

			char name[20];
    	char title[20];

			std::vector<float> resolution;

			std::vector<std::vector<Float_t>> negative_phi(10);
			std::vector<std::vector<Float_t>> positive_phi(10);
			std::vector<Float_t> Psi_EP_posVec(10);
			std::vector<Float_t> Psi_EP_negVec(10);


			TH1F *h1st = new TH1F("h1temp", "Flow, #phi > 0 and #Psi_{EP} < 0", 10, 0.0, 10);
    	h1st->GetXaxis()->SetTitle("Centrality bin");
   		h1st->GetYaxis()->SetTitle("v_{n}");

			TH1F *h2st = new TH1F("h2temp", "Flow, #phi < 0 and #Psi_{EP} > 0", 10, 0.0, 10);
    	h2st->GetXaxis()->SetTitle("Centrality bin");
   		h2st->GetYaxis()->SetTitle("v_{n}");

			TH1F *hist[10];

	/* End of init */

	// Init output histograms, profile, etc.

	// Reading events

	Long64_t nEvents = chain->GetEntries();
	std::cout << "Read " << nFiles << " files. " << nEvents << " events total." << std::endl;


	for (Int_t i = 0; i < 10; i++) {
		sprintf(name, "%.2f < bimp < %.2f", bimpValues[i], bimpValues[i + 1]);
		sprintf(title, "Psi_Ep, %d", i + 1);
		hist[i] = new TH1F(name, title, 200, -1.0, 1.0);
	}

// Collisions
	for (Long64_t ievent = 0; ievent < nEvents; ievent++)
	{
		if (ievent % 1000 == 0) std::cout << "Event [" << ievent << "/" << nEvents << "]" << std::endl;
		chain->GetEntry(ievent);

		// Do event-wise stuff (fill histograms, etc.)

		for (Int_t i = 0; i < 10; i++) {

			if(bimpValues[i] < bimp &&  bimp < bimpValues[i + 1]) {

				// Reading particles from one collision
				for (int itrack = 0; itrack < nh; itrack++)
				{
					// Do particle-wise stuff (fill histograms, etc.)

					/* Conditions */
					if (charge[itrack] != 0) {

						    momModule = sqrt(pow(momx[itrack], 2) + pow(momy[itrack], 2) + pow(momz[itrack], 2));
						    etta = 0.5 * log((momModule + momz[itrack]) / (momModule - momz[itrack]));
						    momTransverse = sqrt(pow(momx[itrack], 2) + pow(momy[itrack], 2));

						    if (- 1.5 < etta && etta < -0.1) {
							if (0.15 < momTransverse &&  momTransverse < 2.0) {
							    phi = PHI(momy[itrack], momx[itrack]);
							    Qx_neg += Qx(phi);
							    Qy_neg += Qy(phi);
							    qW += WEIGHT;
							    flag = true;
									negative_phi[i].push_back(phi);
							}
						    }

						    if (0.1 < etta && etta < 1.5) {
							if (0.15 < momTransverse &&  momTransverse < 2.0) {
							    phi = PHI(momy[itrack], momx[itrack]);
							    Qx_pos += Qx(phi);
							    Qy_pos += Qy(phi);
							    qW += WEIGHT;
							    flag = true;
									positive_phi[i].push_back(phi);
							}
						    }

						}
						etta = 0.0;
						momModule = 0.0;
						momTransverse = 0.0;
						phi = 0.0;
						/* End of conditions */
					}

					Psi_EP_neg = PSI_EP(Qy_neg, Qx_neg);
					Psi_EP_pos = PSI_EP(Qy_pos, Qx_pos);
					
					Psi_EP_negVec.append(Psi_EP_neg);
					Psi_EP_posVec.append(Psi_EP_pos);
					if (flag) {
						if (qW != 0) {

						    pre_Resolution = cos(HARMONICNUMBER * (Psi_EP_neg - Psi_EP_pos));
						    //pre_Resolution /= qW;
						    hist[i]->Fill(pre_Resolution);

				} else { std::cout << "Warning! qW variable has zero quantity! The program will stop immediately!" << std::endl; break; }
			    }

			    error = 0.0;
			    flag = false;
			    qW = 0.0;
			    Qx_neg = 0.0;
			    Qy_neg = 0.0;
			    Qx_pos = 0.0;
			    Qy_pos = 0.0;
			    Psi_EP_neg = 0.0;
			    Psi_EP_pos = 0.0;
			    pre_Resolution = 0.0;
			   //if (Cut(ientry) < 0) continue;
			}
		}
	}
	// Fill histograms here:


	for (Int_t i = 0; i < 10; i++) {
		vnobs_neg = vnobs(positive_phi[i], Psi_EP_negVec[i]);
		vnobs_pos = vnobs(negative_phi[i], Psi_EP_posVec[i]);
		//TMath::sqrt(hist[i]->GetMean()); Resolution
		h1st->SetBinContent(i + 1, vnobs_neg / TMath::sqrt(hist[i]->GetMean()));
		h2st->SetBinContent(i + 1, vnobs_pos / TMath::sqrt(hist[i]->GetMean()));
		vnobs_neg = 0.0;
		vnobs_pos = 0.0;
	}



	// Save output to the file
	std::cout << "Save output information to the file: " << outFileName.Data() << std::endl;
	fo->cd();

	// Save histograms, etc. like histogram->Write();
	h1st->GetYaxis()->SetRangeUser(-0.5, 0.5);
	h2st->GetYaxis()->SetRangeUser(-0.5, 0.5);
	h1st->Write();
	h2st->Write();

	fo->Close();

	std::cout << "Program is finished successfully!" << std::endl;

	// Print out timer info
	timer.Stop();
	timer.Print();


}

